from match_to_py import MSSM
from utils.io import read_param_values, write_to_wcxf
from utils.core import exec
from utils.benchmarks import *


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


#################
#   Benchmarks 
#################

# Writing all operators, WCxf (seq) vs native yaml (seq) 
# wcxf_vs_yaml_all_seq(eft_info, models[0], param_dicts[0], **scale)

# Writing all operators, native yaml, seq vs par
# yaml_seq_vs_par_all(eft_info, models[0], param_dicts[0], **scale) 

# Writing 59 operators of the SMEFT-d6 Warsaw basis
# wcxf_vs_yaml_SMEFTd6_59_seq(models[0], param_dicts[0], **scale)

# Single coefficient execution times
# exec(models[0], "cH", **scale)
# exec(models[0], "cuH_33", **scale)
# exec(models[0], "cqd1_2213", **scale)


####################
#   File Writing 
####################

print("Writing to files.")
for i in range(5):
    t_i = perf_counter()
    f_name = "./benchmark-points/wcxf_all_" + str(i+1) + ".yaml" 
    write_to_wcxf(f_name,models[i],eft_info,**scale) 
    t_f = perf_counter()
    print(f"File {i+1} written in {t_f - t_i:.2f}s.")