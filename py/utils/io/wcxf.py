import wcxf, yaml, json5
import numpy as np
from ..core import eval_wc, create_tasks
from ..io.yaml import complex_representer


def matchete_to_wcxf_name(full_name: str) -> str:
    """returns the WCxf counterpart of a Wilson coefficient defined in the Matchete EFT file"""
    
    name_idx = full_name.split("_")
    name = name_idx[0]

    if name in blv_ops["matchete"].keys():
        wcxf_name = blv_ops["matchete"][name]
    else:
        wcxf_name = name.replace("t", "tilde").replace("H", "phi")[1:]

    if len(name_idx) == 1:
        return wcxf_name
    else:
        return wcxf_name + "_" + name_idx[1]


def wcxf_to_matchete_name(name: str) -> str:
    """returns the counterpart in Matchete format of a Wilson coefficient defined in the WCxf convention"""
   
    if name in blv_ops["wcxf"].keys():
        return blv_ops["wcxf"][name]
    else:
        return "c" + name.replace("tilde", "t").replace("phi", "H") 


def permute_chars(seq: str, ordering: list[int]) -> str:
    """returns a string produced by permuting the characters in an input string based on the specified ordering"""
    new_seq = ""
    for i in ordering:
        new_seq += seq[i-1]

    return new_seq


def collect_permutations(arr1: np.ndarray, arr2: np.ndarray) -> np.ndarray:
    """returns an array containing all non-null arrays provided in the 2 input arrays"""

    assert(len(arr1) and len(arr2)) 

    if np.ndim(arr1) == np.ndim(arr2):
        return np.array([arr1, arr2])
    elif np.ndim(arr1) > np.ndim(arr2):
        return np.concatenate((arr1, [arr2]))
    else:
        return np.concatenate(([arr1], arr2))    


def generate_evaluation_list(wcs_wcxf: list[str], wc_info: dict) -> tuple[dict[str, int], list[str]]:
    """Given a list of Wilson coeffcient names (in the WCxf convention) and permutation symmetry information about 
    each coefficient, returns a tuple containing the list of unique index permutations for each case and a list of 
    Wilson coefficient names in Matchete's convention"""

    unique_ops = {}
    wcs_matchete = []
    # duplicates = {}

    for wc in wcs_wcxf:
        name_idx = wc.split("_")
        m_name = wcxf_to_matchete_name(name_idx[0])
        
        if len(name_idx) == 1: # i.e., no flavour index (purely bosonic operators)
            wcs_matchete.append(m_name)
            unique_ops[wc] = 1
            # duplicates[wc] = [m_name]

        else:
            idx_seq = name_idx[1]
            wcs_matchete.append(m_name + "_" + idx_seq)
            
            conj = wc_info[m_name].get("conj")
            symm = wc_info[m_name].get("symm")
            
            match (conj, symm):
                case (None, None): 
                    perms = [] 
                case (None, x) | (x, None): 
                    perms = [x]
                case (x, y): 
                    perms = collect_permutations(x, y)
                
            idx_perms = [permute_chars(idx_seq, perm) for perm in perms]
            idx_perms = idx_perms + [permute_chars(idx_perm, perm) for perm in perms for idx_perm in idx_perms]
            
            idx_perms.append(idx_seq)
            unique_perms = list(set(idx_perms))
            unique_ops[wc] = len(unique_perms)
            # duplicates[wc] = [m_name + "_" + unique_perm for unique_perm in unique_perms]

    return (unique_ops, wcs_matchete)


def write_to_wcxf(filename: str, model, eft_info: dict[str, str], wc_names: list[str] | None = None, opt: str = "par") -> None:
    """
    Writes the values of evaluated Wilson coefficients to a .yaml file in the WCxf convention

    Parameters
    ----------
        filename : str
            A string specifying the path of the output .yaml file
        model
            Instance of the UV model class defined in the match_to_py module
        eft_info : dict
            A dictionary specifying the "eft" and the "basis".
        wc_names : list[str] | None, optional
            A list of Wilson coefficient names in the WCxf convention. The default option is 
            to evaluate the entire list of coefficients in the specified EFT-basis. 
        opt : str, "seq" | "par"
            Option for specifying whether the operation should proceed sequentially, or 
            in parallel with OpenMP at the C++ backend 
        kw : dict
            Fixed global parameters such as "mubarsq" for renormalization scale and
            "hbar" for the matching order, specified using a dictionary
    
    """
    match opt:
        case "seq":
            write_to_wcxf_seq(filename, model, eft_info, wc_names)
        case "par":
            write_to_wcxf_par(filename, model, eft_info, wc_names)


def write_to_wcxf_seq(filename: str, model, eft_info: dict[str, str], wc_names: list[str] | None = None) -> None:
    """
    Writes the values of evaluated Wilson coefficients to a .yaml file in the WCxf convention, follows sequential evaluation.
    """
    output_dict = {
        "eft": eft_info["eft"], 
        "basis": eft_info["basis"], 
        "scale": model.getScale()
    }

    eft_basis = wcxf.Basis[eft_info["eft"], eft_info["basis"]]
    wcs_wcxf = eft_basis.all_wcs

    wc_info_path = "EFT_db" + "/" + eft_info["eft"] + "/" + eft_info["basis"] + "/wc_info.json"
    with open(wc_info_path, "r") as f:
        wc_info = json5.load(f)

    if not wc_names:
        wc_names = wcs_wcxf
    
    unique_ops, wcs_matchete = generate_evaluation_list(wc_names, wc_info["permutations"])
    values_1_wc = [eval_wc(model, wc) for wc in wcs_matchete]

    values = [unique_ops[wc_names[i]] * values_1_wc[i] for i in range(len(wc_names))]
    values_dict = dict(zip(wc_names, values))
    
    output_dict["values"] = values_dict

    yaml.add_representer(complex, complex_representer)
    with open(filename,'w') as out_file:
        yaml.dump(output_dict, out_file, sort_keys=False)


def write_to_wcxf_par(filename: str, model, eft_info: dict[str, str], wc_names: list[str] | None = None) -> None:
    """
    Writes the values of evaluated Wilson coefficients to a .yaml file in the WCxf convention, utilizes parallelisation
    through OpenMP and batch processing at the C++ end
    """
    output_dict = {
        "eft": eft_info["eft"], 
        "basis": eft_info["basis"], 
        "scale": model.getScale()
    }

    eft_basis = wcxf.Basis[eft_info["eft"], eft_info["basis"]]
    wcs_wcxf = eft_basis.all_wcs

    wc_info_path = "EFT_db" + "/" + eft_info["eft"] + "/" + eft_info["basis"] + "/wc_info.json"
    with open(wc_info_path, "r") as f:
        wc_info = json5.load(f)

    if not wc_names:
        wc_names = wcs_wcxf
    
    unique_ops, wcs_matchete = generate_evaluation_list(wc_names, wc_info["permutations"])
    tasks = create_tasks(model, wcs_matchete)
    res_dict = model.batch_eval(tasks)

    values_1_wc = list(res_dict.values())
    keys_wc = [matchete_to_wcxf_name(k) for k in list(res_dict.keys())]

    values = [unique_ops[keys_wc[i]] * values_1_wc[i] for i in range(len(wc_names))]
    values_dict = dict(zip(keys_wc, values))
    
    output_dict["values"] = values_dict

    yaml.add_representer(complex, complex_representer)
    with open(filename,'w') as out_file:
        yaml.dump(output_dict, out_file, sort_keys=False)


blv_ops = {
    "matchete": {"cduq": "duql", "cqqu": "qque", "cqqq": "qqql", "cduu": "duue"},
    "wcxf": {"duql": "cduq", "qque": "cqqu", "qqql": "cqqq", "duue": "cduu"},
}
