from match_to_py import MSSM
from utils.io import write_to_wcxf

# define a parameter dictionary

param_dict = {
    "g1": 0.37, "g3": 1.1, "cgamma": 0.01, "Yu11": 0.00001, "Yu22": 0.007, "Yu33": 0.9,
    "met1": 1e9, "met2": 1e9, "met3": 1e9, "mlt1": 1e9, "mlt2": 1e9, "mlt3": 1e9,
    "mqt1": 1e9, "mqt2": 1e9, "mqt3": 1e9, "mut1": 1e9, "mut2": 1e9, "mut3": 2e3,
    "mdt1": 1e9, "mdt2": 1e9, "mdt3": 1e9, "m1": 1.2e3, "m2": 1e9, "m3": 1e9,
    "mPhi": 1e9, "muTilde": 1e9
}

# create an instance of the MSSM class
model1 = MSSM(param_dict, 1e3, True)

# evaluate Wilson coefficients as method calls
print(model1.cG())
print(model1.cuB(2,2))

# evaluate coefficients and write to a WCxf file

wcs = ["uu_1331", "G", "phiG", "phiD", "uG_33"]
eft_info = { "eft": "SMEFT", "basis": "Warsaw" }
write_to_wcxf("example_wcxf.yaml",model1,eft_info,wcs,opt="seq")
write_to_wcxf("example_wcxf_all.yaml",model1,eft_info)
