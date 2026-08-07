from match_to_py import MSSM
from utils.core import exec
from utils.benchmarks import *


###########################
#   parameter definitions
###########################

param_dict1 = {
    "g1": 0.37, "g3": 1.1, "cgamma": 0.01, "Yu11": 0.00001, "Yu22": 0.007, "Yu33": 0.9,
    "met1": 1e9, "met2": 1e9, "met3": 1e9, "mlt1": 1e9, "mlt2": 1e9, "mlt3": 1e9,
    "mqt1": 1e9, "mqt2": 1e9, "mqt3": 1e9, "mut1": 1e9, "mut2": 1e9, "mut3": 2e3,
    "mdt1": 1e9, "mdt2": 1e9, "mdt3": 1e9, "m1": 1.2e3, "m2": 1e9, "m3": 1e9,
    "mPhi": 1e9, "muTilde": 1e9
}


###########################
#   Model Initialization
###########################

model1 = MSSM(param_dict1, 1e3, True)
eft_info = {"eft": "SMEFT", "basis": "Warsaw"}


#################
#   Benchmarks
#################

# Writing all operators, WCxf (seq) vs native yaml (seq)
wcxf_vs_yaml_all_seq(eft_info, model1, param_dict1)

# Writing all operators, WCxf seq vs par
wcxf_seq_vs_par_all(eft_info, model1)

# Writing all operators, native yaml, seq vs par
yaml_seq_vs_par_all(eft_info, model1, param_dict1)

# Writing 59 operators of the SMEFT-d6 Warsaw basis, sequentially
wcxf_vs_yaml_SMEFTd6_59_seq(model1, param_dict1)

# Writing 59 operators of the SMEFT-d6 Warsaw basis, in parallel
wcxf_vs_yaml_SMEFTd6_59_par(model1, param_dict1)


# Single coefficient execution times
exec(model1, "cH")
exec(model1, "cuH_33")
exec(model1, "cqd1_2213")