import wcxf, yaml
import numpy as np
from math import sqrt
from ..core import eval_wc
from ..io.yaml import float_representer
from ..io.yaml import complex_representer

# globally loaded list of all Wilson coefficient names included in the wcxf SMEFT-Warsaw basis.
smeft_warsaw = wcxf.Basis["SMEFT","Warsaw"]
wcs_wcxf = smeft_warsaw.all_wcs


def matchete_to_wcxf_name(name: str) -> str:
    """returns the wcxf counterpart of a Wilson coefficient defined in the Matchete EFT file"""
    
    if name in blv_ops["matchete"].keys():
        return blv_ops["matchete"][name]
    else:
        return name.replace("t", "tilde").replace("H", "phi")[1:]


def wcxf_to_matchete_name(name: str) -> str:
    """returns the counterpart in Matchete format of a Wilson coefficient defined in the wcxf convention"""
   
    if name in blv_ops["wcxf"].keys():
        return blv_ops["matchete"][name]
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
    

def permute_and_compare(name: str, idx_seq: str, perms: list[list[int]]) -> str | None:
    """
    Generates permutations of the flavour indices for a specified Wilson coefficient and
    returns the one available in the wcxf SMEFT-Warsaw list
    """
    
    global wcs_wcxf
    try:
        for idx in idx_seq:
            assert (int(idx)>0 and int(idx)<4)
    except AssertionError:
        print(f"Invalid flavour index in {name}_{idx_seq}, skipping")
        return
        
    if (name + "_" + idx_seq) in wcs_wcxf:
        # if the given sequence itself is present in the list then permutations are not needed
        return idx_seq 

    all_seq = [idx_seq]
    while True:
        new_perm = False
        for perm in perms:
            new_seq = permute_chars(all_seq[-1], perm)
            if (name + "_" + new_seq) in wcs_wcxf:
                return new_seq
            if new_seq != all_seq[-1]:
                new_perm = True
                all_seq.append(new_seq)
    
        if not new_perm:
            return 


def wcxf_name_val(model, wc_name: str, convention: str, **kw: dict) -> tuple[str, float | complex] | None:
    """
    evaluates a Wilson coefficient given in Matchete or wcxf convention and returns a tuple containing
    the coefficient name in wcxf format and the value 
    """
    name_idx = wc_name.split("_")
    name = name_idx[0]

    if convention == "matchete":
        m_name = name
        wcxf_name = matchete_to_wcxf_name(name)
    elif convention == "wcxf":
        m_name = wcxf_to_matchete_name(name)
        wcxf_name = name
    else:
        print("Invalid naming convention for Wilson coefficients")
        return
    
    if len(name_idx) == 1: # i.e., no flavour index (purely bosonic operators)
        return (wcxf_name, eval_wc(model,m_name,**kw))
    else:
        idx_seq = name_idx[1]
        conj = wc_info[m_name].get("conj")
        symm = wc_info[m_name].get("symm")
        
        match (conj, symm):
            case (None, None): 
                return (
                    wcxf_name + "_" + idx_seq, 
                    eval_wc(model, m_name + "_" + idx_seq, **kw)
                ) 
            case (None, x) | (x, None): 
                perms = [x]
            case (x, y): 
                perms = collect_permutations(x, y)
            
        updated_seq = permute_and_compare(wcxf_name, idx_seq, perms)
        if not updated_seq:
            return
        else:
            return (
                wcxf_name + "_" + updated_seq, 
                eval_wc(model, m_name + "_" + updated_seq, **kw) 
            )


def create_wcxf_dict(model, wc_names: list[str], convention: str, **kw: dict) -> dict[str, float | complex]:
    """
    creates a dictionary of Wilson coefficient name-value pairs for a specified model and a list of coefficient names
    """
    wc_dict = dict()
    for wc in wc_names:
        k, v = wcxf_name_val(model,wc,convention,**kw)
        wc_dict[k] = v

    return wc_dict


def write_to_wcxf(filename: str, model, wc_names: list[str], convention: str, **kw:dict) -> None:
    """
    Writes the values of evaluated Wilson coefficients to a .yaml file in the wcxf convention

    Parameters
    ----------
        filename : str
            A string specifying the path of the output .yaml file
        model
            Instance of the UV model class defined in the match_to_py module
        wc_names : list[str]
            A list of Wilson coefficient names in Matchete or wcxf convention
        convention : str, "matchete" | "wcxf"
            A flag to specify if the the entries of wc_names is provided in Matchete 
            or wcxf convention
        kw : dict
            Fixed global parameters such as "mubarsq" for renormalization scale and
            "hbar" for the matching order, specified using a dictionary
    
    """
    output_dict = {
        "eft": "SMEFT", 
        "basis": "Warsaw", 
        "scale": sqrt(kw["mubarsq"])*1000
        # convert scale from TeV to GeV
    }

    yaml.add_representer(float, float_representer)
    yaml.add_representer(complex, complex_representer)

    values_dict = create_wcxf_dict(model,wc_names,convention,**kw)
    output_dict["values"] = values_dict

    with open(filename,'w') as out_file:
        yaml.dump(output_dict, out_file, sort_keys=False)


wc_info = {
    "cuG": {},  
    "cuW": {},  
    "cuB": {},  
    "cdG": {},  
    "cdW": {},  
    "cdB": {},  
    "ceW": {},  
    "ceB": {},  
    "cHl1": { "conj": [2,1] }, 
    "cHl3": { "conj": [2,1] }, 
    "cHe": { "conj": [2,1] }, 
    "cHq1": { "conj": [2,1] }, 
    "cHq3": { "conj": [2,1] }, 
    "cHu": { "conj": [2,1] },  
    "cHd": { "conj": [2,1] }, 
    "cHud": {}, 
    "cll": { "conj": [2,1,4,3], "symm": [3,4,1,2] }, 
    "cqq1": { "conj": [2,1,4,3], "symm": [3,4,1,2] }, 
    "cqq3": { "conj": [2,1,4,3], "symm": [3,4,1,2] }, 
    "clq1": { "conj": [2,1,4,3] }, 
    "clq3": { "conj": [2,1,4,3] },  
    "cee": { "conj": [2,1,4,3], "symm": [[3,4,1,2], [1,4,3,2]] }, 
    "cuu": { "conj": [2,1,4,3], "symm": [3,4,1,2] },   
    "cdd": { "conj": [2,1,4,3], "symm": [3,4,1,2] }, 
    "ceu": { "conj": [2,1,4,3] }, 
    "ced": { "conj": [2,1,4,3] }, 
    "cud1": { "conj": [2,1,4,3] }, 
    "cud8": { "conj": [2,1,4,3] }, 
    "cle": { "conj": [2,1,4,3] },   
    "clu": { "conj": [2,1,4,3] },   
    "cld": { "conj": [2,1,4,3] },   
    "cqe": { "conj": [2,1,4,3] },   
    "cqu1": { "conj": [2,1,4,3] },   
    "cqu8": { "conj": [2,1,4,3] },   
    "cqd1": { "conj": [2,1,4,3] },   
    "cqd8": { "conj": [2,1,4,3] },   
    "cledq": {},   
    "cquqd1": {},   
    "cquqd8": {},   
    "clequ1": {},   
    "clequ3": {},   
    "cduq": {},   
    "cqqu": { "symm": [2,1,3,4] },   
    "cqqq": {},   
    "cduu": {},   
}

blv_ops = {
    "matchete": {"cduq": "duql", "cqqu": "qque", "cqqq": "qqql", "cduu": "duue"},
    "wcxf": {"duql": "cduq", "qque": "cqqu", "qqql": "cqqq", "duue": "cduu"},
}
