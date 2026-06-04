from utils.io import write_to_wcxf, write_to_yaml
from utils.io.wcxf import wcxf_to_matchete_name
from time import perf_counter
import wcxf


def wcxf_vs_yaml_all_seq(eft_info, model, param_dict) -> None:
    """
        Compares the execution times for evaluating and writing to file all Wilson coefficients
        defined within a specifc EFT-basis in (i) native yaml (ii) WCxf format
    """
    
    # the WCxf writer
    f_name = "./benchmark-points/all_wcxf.yaml"
    t_i = perf_counter() 
    write_to_wcxf(f_name,model,eft_info,opt="seq") 
    t_f = perf_counter()
    print(f"File written in WCxf format (using sequential process) in {t_f - t_i:.2f}s.")

    # the native yaml writer for all operators defined in the WCxf basis
    eft_basis = wcxf.Basis[eft_info["eft"], eft_info["basis"]]
    wcs_wcxf = eft_basis.all_wcs
    wc_names: list[str] = list()
    for wc in wcs_wcxf:
        name = wc.split("_")
        new_name = wcxf_to_matchete_name(name[0])
        if len(name) == 1:
            wc_names.append(new_name)
        else:
            wc_names.append(new_name + "_" + name[1])   

    f_name = "./benchmark-points/all_native.yaml" 
    t_i = perf_counter()
    keys = [[],wc_names]
    write_to_yaml(f_name,model,param_dict,keys,"seq") 
    t_f = perf_counter()
    print(f"File written in native format (using sequential process) in {t_f - t_i:.2f}s.")


def wcxf_vs_yaml_SMEFTd6_59_seq(model, param_dict) -> None:
    """
        Compares the execution times for evaluating and writing to file the 59 Wilson coefficients
        of (single fermion generation) dimension 6 SMEFT Warsaw basis (i) native yaml (ii) WCxf format
    """
    
    # 59 Wilson coefficients in WCxf format
    wc_names_wcxf = ["G", "W", "phi", "phiBox", "phiD", "phiG", "phiB", "phiW", "phiWB", "ephi_33", "dphi_33", "uphi_33", "eB_33", "eW_33", "dB_33", "dW_33", "dG_33", "uB_33", "uW_33", "uG_33", "phie_33", "phiu_33", "phid_33", "phiud_33","phil1_33", "phil3_33", "phiq1_33", "phiq3_33", "ll_3333", "qq1_3333", "qq3_3333", "lq1_3333", "lq3_3333", "ee_3333", "uu_3333", "dd_3333", "eu_3333", "ed_3333", "ld_3333", "lu_3333", "qe_3333", "qu1_3333", "qu8_3333", "qd1_3333", "qd8_3333", "ud1_3333", "ud8_3333", "quqd1_3333", "quqd8_3333", "ledq_3333", "lequ1_3333", "lequ3_3333"]

    
    eft_info = {"eft": "SMEFT", "basis": "Warsaw"}
    
    f_name = "./benchmark-points/59_wcxf.yaml"
    t_i = perf_counter() 
    write_to_wcxf(f_name,model,eft_info,wc_names_wcxf,"seq")
    t_f = perf_counter()
    print(f"File written in WCxf format in {t_f - t_i:.2f}s.")

    # 59 Wilson coefficients in Matchete format
    wc_names_matchete = ["cG", "cW", "cHBox", "cH", "cHD", "cHG", "cHB", "cHW", "cHWB", "ceH_33", "cdH_33", "cuH_33", "ceB_33", "ceW_33", "cdB_33", "cdW_33", "cdG_33", "cuB_33", "cuW_33", "cuG_33", "cHe_33", "cHu_33", "cHd_33", "cHud_33","cHl1_33", "cHl3_33", "cHq1_33", "cHq3_33", "cll_3333", "cqq1_3333", "cqq3_3333", "clq1_3333", "clq3_3333", "cee_3333", "cuu_3333", "cdd_3333", "ceu_3333", "ced_3333", "cld_3333", "clu_3333", "cqe_3333", "cqu1_3333", "cqu8_3333", "cqd1_3333", "cqd8_3333", "cud1_3333", "cud8_3333", "cquqd1_3333", "cquqd8_3333", "cledq_3333", "clequ1_3333", "clequ3_3333"] 

    f_name = "./benchmark-points/59_native.yaml" 
    t_i = perf_counter()
    keys = [[],wc_names_matchete]
    write_to_yaml(f_name,model,param_dict,keys,"seq") 
    t_f = perf_counter()
    print(f"File written in native format in {t_f - t_i:.2f}s.")


def wcxf_vs_yaml_SMEFTd6_59_par(model, param_dict) -> None:
    """
        Compares the execution times for evaluating and writing to file the 59 Wilson coefficients
        of (single fermion generation) dimension 6 SMEFT Warsaw basis (i) native yaml (ii) WCxf format
    """
    
    # 59 Wilson coefficients in WCxf format
    wc_names_wcxf = ["G", "W", "phi", "phiBox", "phiD", "phiG", "phiB", "phiW", "phiWB", "ephi_33", "dphi_33", "uphi_33", "eB_33", "eW_33", "dB_33", "dW_33", "dG_33", "uB_33", "uW_33", "uG_33", "phie_33", "phiu_33", "phid_33", "phiud_33","phil1_33", "phil3_33", "phiq1_33", "phiq3_33", "ll_3333", "qq1_3333", "qq3_3333", "lq1_3333", "lq3_3333", "ee_3333", "uu_3333", "dd_3333", "eu_3333", "ed_3333", "ld_3333", "lu_3333", "qe_3333", "qu1_3333", "qu8_3333", "qd1_3333", "qd8_3333", "ud1_3333", "ud8_3333", "quqd1_3333", "quqd8_3333", "ledq_3333", "lequ1_3333", "lequ3_3333"]

    
    eft_info = {"eft": "SMEFT", "basis": "Warsaw"}
    
    f_name = "./benchmark-points/59_wcxf_par.yaml"
    t_i = perf_counter() 
    write_to_wcxf(f_name,model,eft_info,wc_names_wcxf)
    t_f = perf_counter()
    print(f"File written in WCxf format in {t_f - t_i:.2f}s.")

    # 59 Wilson coefficients in Matchete format
    wc_names_matchete = ["cG", "cW", "cHBox", "cH", "cHD", "cHG", "cHB", "cHW", "cHWB", "ceH_33", "cdH_33", "cuH_33", "ceB_33", "ceW_33", "cdB_33", "cdW_33", "cdG_33", "cuB_33", "cuW_33", "cuG_33", "cHe_33", "cHu_33", "cHd_33", "cHud_33","cHl1_33", "cHl3_33", "cHq1_33", "cHq3_33", "cll_3333", "cqq1_3333", "cqq3_3333", "clq1_3333", "clq3_3333", "cee_3333", "cuu_3333", "cdd_3333", "ceu_3333", "ced_3333", "cld_3333", "clu_3333", "cqe_3333", "cqu1_3333", "cqu8_3333", "cqd1_3333", "cqd8_3333", "cud1_3333", "cud8_3333", "cquqd1_3333", "cquqd8_3333", "cledq_3333", "clequ1_3333", "clequ3_3333"] 

    f_name = "./benchmark-points/59_native_par.yaml" 
    t_i = perf_counter()
    keys = [[],wc_names_matchete]
    write_to_yaml(f_name,model,param_dict,keys,"par") 
    t_f = perf_counter()
    print(f"File written in native format in {t_f - t_i:.2f}s.")


def wcxf_seq_vs_par_all(eft_info, model) -> None:
    """
        Compares the execution times for evaluating and writing to file all Wilson coefficients
        defined within a specifc EFT-basis in the native yaml format for (i) sequential and 
        (ii) parallel execution 
    """
    
    f_name = "./benchmark-points/all_wcxf.yaml"
    t_i = perf_counter() 
    write_to_wcxf(f_name,model,eft_info,opt="seq") 
    t_f = perf_counter()
    print(f"File written in WCxf format (using sequential process) in {t_f - t_i:.2f}s.")

    f_name = "./benchmark-points/all_wcxf_par.yaml"
    t_i = perf_counter() 
    write_to_wcxf(f_name,model,eft_info) 
    t_f = perf_counter()
    print(f"File written in WCxf format (using parallel process) in {t_f - t_i:.2f}s.")


def yaml_seq_vs_par_all(eft_info, model, param_dict) -> None:
    """
        Compares the execution times for evaluating and writing to file all Wilson coefficients
        defined within a specifc EFT-basis in the native yaml format for (i) sequential and 
        (ii) parallel execution
    """
    
    eft_basis = wcxf.Basis[eft_info["eft"], eft_info["basis"]]
    wcs_wcxf = eft_basis.all_wcs
    wc_names: list[str] = list()
    for wc in wcs_wcxf:
        name = wc.split("_")
        new_name = wcxf_to_matchete_name(name[0])
        if len(name) == 1:
            wc_names.append(new_name)
        else:
            wc_names.append(new_name + "_" + name[1])   

    f_name = "./benchmark-points/all_native.yaml" 
    t_i = perf_counter()
    keys = [[],wc_names]
    write_to_yaml(f_name,model,param_dict,keys,"seq") 
    t_f = perf_counter()
    print(f"File written in native format (using sequential process) in {t_f - t_i:.2f}s.")

    f_name = "./benchmark-points/all_native_par.yaml" 
    t_i = perf_counter()
    keys = [[],wc_names]
    write_to_yaml(f_name,model,param_dict,keys,"par") 
    t_f = perf_counter()
    print(f"File written in native format (using parallel C++ process) in {t_f - t_i:.2f}s.")
