from match_to_py import SeesawTypeI
from utils.io import write_to_wcxf

# define a parameter dictionary

param_dict = {
    "gY": 0.36, "gL": 0.63, "lmbd": 0.085, 
    "Yu11": 7e-6, "Yu22": 3.3e-3, "Yu33": 0.86,
    "Yd11": 1.5e-5, "Yd22": 3e-4, "Yd33": 0.015,
    "Ye11": 2.9e-6, "Ye22": 6e-4, "Ye33": 0.01,
    "Yn11": 1.0, "Yn22": 2.0, "Yn33": 3.0,
    "MNR1": 1e8, "MNR2": 2e8, "MNR3": 3e8
}

# create an instance of the SeesawTypeI class
model1 = SeesawTypeI(param_dict, 1e3, True)

# evaluate Wilson coefficients as method calls
print(model1.cH())
print(model1.cuH(2,2))

# evaluate coefficients and write to a WCxf file

wcs = ["uu_1331", "G", "phiG", "phiD", "uG_33"]
eft_info = { "eft": "SMEFT", "basis": "Warsaw" }
write_to_wcxf("example_wcxf_all.yaml",model1,eft_info)
