from match_to_py import SeesawTypeI
from utils.core import exec
from utils.benchmarks import *


###########################
#   parameter definitions
###########################

param_dict = {
    "gY": 0.36, "gL": 0.63, "lmbd": 0.085, 
    "Yu11": 7e-6, "Yu22": 3.3e-3, "Yu33": 0.86,
    "Yd11": 1.5e-5, "Yd22": 3e-4, "Yd33": 0.015,
    "Ye11": 2.9e-6, "Ye22": 6e-4, "Ye33": 0.01,
    "Yn11": 1.0, "Yn22": 2.0, "Yn33": 3.0,
    "MNR1": 1e8, "MNR2": 2e8, "MNR3": 3e8
}

###########################
#   Model Initialization
###########################

model1 = SeesawTypeI(param_dict, 1e3, True)
eft_info = {"eft": "SMEFT", "basis": "Warsaw"}


#################
#   Benchmarks
#################

# Writing all operators, WCxf (seq) vs native yaml (seq)
wcxf_vs_yaml_all_seq(eft_info, model1, param_dict)

# Writing all operators, WCxf seq vs par
wcxf_seq_vs_par_all(eft_info, model1)

# Writing all operators, native yaml, seq vs par
yaml_seq_vs_par_all(eft_info, model1, param_dict)

# Writing 59 operators of the SMEFT-d6 Warsaw basis, sequentially
wcxf_vs_yaml_SMEFTd6_59_seq(model1, param_dict)

# Writing 59 operators of the SMEFT-d6 Warsaw basis, in parallel
wcxf_vs_yaml_SMEFTd6_59_par(model1, param_dict)


# Single coefficient execution times
exec(model1, "cH")
exec(model1, "cuH_33")
exec(model1, "cqu8_3333")
