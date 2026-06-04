from match_to_py import MSSM
from utils.io import write_to_wcxf

# define a parameter dictionary
param_dict = {
    "g1": 0.37, "g3": 1.1, "cgamma": 0.01,
    "Yu11": 0.00001, "Yu22": 0.007, "Yu33": 0.9,
    "mPhi": 1001000, "muTilde": 1000000,
    "m1": 1.2, "m2": 1002000, "m3": 1003000,
    "met1": 1004000, "met2": 1005000, "met3": 1006000,
    "mlt1": 1007000, "mlt2": 1008000, "mlt3": 1009000,
    "mqt1": 1010000, "mqt2": 1011000, "mqt3": 1012000,
    "mut1": 1013000, "mut2": 1014000, "mut3": 2.0,
    "mdt1": 1016000, "mdt2": 1017000, "mdt3": 1018000,
}

# create an instance of the MSSM class
model1 = MSSM(param_dict, 1.0, True)

# evaluate Wilson coefficients as method calls
print(model1.cG())
print(model1.cuB(2,2))

# evaluate coefficients and write to a WCxf file

wcs = ["uu_1331", "G", "phiG", "phiD", "uG_33"]
eft_info = { "eft": "SMEFT", "basis": "Warsaw" }
write_to_wcxf("example_wcxf.yaml",model1,eft_info,wcs,opt="seq")
write_to_wcxf("example_wcxf_all.yaml",model1,eft_info)
