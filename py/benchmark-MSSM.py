from operator import call
from match_to_py import MSSM
from utils.io import *
from utils.io.wcxf import wcxf_to_matchete_name
from utils.core import exec
from time import perf_counter
import wcxf

###########################
#   parameter definitions
###########################

# SM parameters (common to all benchmark points)
param_dict_sm = {
    "g1": 0.354474, "g2": 0.639861, "g3": 1.06197,
    "Ye11": 2.93503e-6, "Ye22": 0.00060688, "Ye33": 0.0102062,
    "Yd11": 0.0000269954, "Yd22": 0.000537037, "Yd33": 0.0240259,
    "Yu11": 0.0000120923, "Yu12": -0.00164521, "Yu13": 0.00683461 + 0.00295002j, 
    "Yu21": 2.79157e-6, "Yu22": 0.00712664, "Yu23": -0.0362955, 
    "Yu31": 1.9044e-8 + 4.21698e-8j, "Yu32": 0.000305777,  "Yu33": 0.867899,
}

print("Loading data from file.")
t1 = perf_counter()

param_dicts = list()
for i in range(1,6):
    filename = "./benchmark-points/bp_" + str(i) + ".yaml"
    par_dict_mssm, _ = read_param_values(filename)
    par_dict = param_dict_sm | par_dict_mssm
    param_dicts.append(par_dict)

t2 = perf_counter()
print(f"Data for {len(param_dicts)} benchmark points loaded in {t2 - t1:.2f}s.")


###########################
#   Model Initialization 
###########################

print("Initializing MSSM model instances.")
t3 = perf_counter()
models = list()
for p_dict in param_dicts:
    models.append(MSSM(p_dict))

t4 = perf_counter()
print(f"Initialized {len(models)} models in {t4 - t3:.2f}s.")

scale = {"mubarsq": 1e6, "hbar": 0.006332574}
eft_info = {"eft": "SMEFT", "basis": "Warsaw"}

###############################################################
#     Total time for writing all operators, WCxf vs native 
###############################################################

# the WCxf writer
# f_name = "./benchmark-points/all_wcxf.yaml"
# t_i = perf_counter() 
# write_to_wcxf(f_name,models[0],eft_info,**scale) 
# t_f = perf_counter()
# print(f"File written in WCxf format in {t_f - t_i:.2f}s.")

# the native yaml writer for all operators defined in the WCxf basis
# smeft_warsaw = wcxf.Basis["SMEFT","Warsaw"]
# wcs_wcxf = smeft_warsaw.all_wcs
# wc_names: list[str] = list()
# for wc in wcs_wcxf:
#     name = wc.split("_")
#     new_name = wcxf_to_matchete_name(name[0])
#     if len(name) == 1:
#         wc_names.append(new_name)
#     else:
#         wc_names.append(new_name + "_" + name[1])   

# f_name = "./benchmark-points/all_native.yaml" 
# t_i = perf_counter()
# keys = [[],wc_names]
# write_to_yaml(f_name,models[0],param_dicts[0],keys,"seq",**scale) 
# t_f = perf_counter()
# print(f"File written in native format (using sequential process) in {t_f - t_i:.2f}s.")


########################################################################
#     Total time for writing all operators, parallel (C++) native 
########################################################################

# f_name = "./benchmark-points/all_native_par_c.yaml" 
# t_i = perf_counter()
# keys = [[],wc_names]
# write_to_yaml(f_name,models[0],param_dicts[0],keys,"par",**scale) 
# t_f = perf_counter()
# print(f"File written in native format (using parallel C++ process) in {t_f - t_i:.2f}s.")


# ###############################################################
#     Total time for writing 59 operators, WCxf vs native 
# ###############################################################

# # 59 Wilson coefficients in WCxf format
# wc_names_wcxf = ["G", "W", "phi", "phiBox", "phiD", "phiG", "phiB", "phiW", "phiWB", "ephi_33", "dphi_33", "uphi_33", "eB_33", "eW_33", "dB_33", "dW_33", "dG_33", "uB_33", "uW_33", "uG_33", "phie_33", "phiu_33", "phid_33", "phiud_33","phil1_33", "phil3_33", "phiq1_33", "phiq3_33", "ll_3333", "qq1_3333", "qq3_3333", "lq1_3333", "lq3_3333", "ee_3333", "uu_3333", "dd_3333", "eu_3333", "ed_3333", "ld_3333", "lu_3333", "qe_3333", "qu1_3333", "qu8_3333", "qd1_3333", "qd8_3333", "ud1_3333", "ud8_3333", "quqd1_3333", "quqd8_3333", "ledq_3333", "lequ1_3333", "lequ3_3333"]

# f_name = "./benchmark-points/59_wcxf.yaml"
# t_i = perf_counter() 
# write_to_wcxf(f_name,models[0],eft_info,wc_names_wcxf,**scale)
# t_f = perf_counter()
# print(f"File written in WCxf format in {t_f - t_i:.2f}s.")

# # 59 Wilson coefficients in Matchete format
# wc_names_matchete = ["cG", "cW", "cHBox", "cH", "cHD", "cHG", "cHB", "cHW", "cHWB", "ceH_33", "cdH_33", "cuH_33", "ceB_33", "ceW_33", "cdB_33", "cdW_33", "cdG_33", "cuB_33", "cuW_33", "cuG_33", "cHe_33", "cHu_33", "cHd_33", "cHud_33","cHl1_33", "cHl3_33", "cHq1_33", "cHq3_33", "cll_3333", "cqq1_3333", "cqq3_3333", "clq1_3333", "clq3_3333", "cee_3333", "cuu_3333", "cdd_3333", "ceu_3333", "ced_3333", "cld_3333", "clu_3333", "cqe_3333", "cqu1_3333", "cqu8_3333", "cqd1_3333", "cqd8_3333", "cud1_3333", "cud8_3333", "cquqd1_3333", "cquqd8_3333", "cledq_3333", "clequ1_3333", "clequ3_3333"] 

# f_name = "./benchmark-points/59_native.yaml" 
# t_i = perf_counter()
# keys = [[],wc_names_matchete]
# write_to_yaml(f_name,models[0],param_dicts[0],keys,"seq",**scale) 
# t_f = perf_counter()
# print(f"File written in native format in {t_f - t_i:.2f}s.")


#############################################
#    eval_wc() vs direct method call
#############################################

# exec(models[0], "cH", **scale)
# exec(models[0], "cuH_33", **scale)
# exec(models[0], "cqd1_2213", **scale)

# t_i = perf_counter()
# x = models[0].cH(**scale)
# t_f = perf_counter()
# print(f"{x}, {t_f-t_i}s")

# t_i = perf_counter()
# x = models[0].cuH(2,2,**scale)
# t_f = perf_counter()
# print(f"{x}, {t_f-t_i}s")

# t_i = perf_counter()
# x = models[0].cqd1(1,1,0,2,**scale)
# t_f = perf_counter()
# print(f"{x}, {t_f-t_i}s")


####################
#   File Writing 
####################

# print("Writing to files.")
# for i in range(5):
#     t_i = perf_counter()
#     f_name = "./benchmark-points/wcxf_all_" + str(i+1) + ".yaml" 
#     write_to_wcxf_all(f_name,models[i],**scale) 
#     t_f = perf_counter()
#     print(f"File {i+1} written in {t_f - t_i:.2f}s.")