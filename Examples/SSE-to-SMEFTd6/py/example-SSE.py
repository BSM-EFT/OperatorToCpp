from match_to_py import SingletScalarExtension
from utils.io import write_to_wcxf

# define a parameter dictionary

param_dict = {
    "gY": 0.36, "gL": 0.63, "lmbd": 0.085, 
    "Yu11": 7e-6, "Yu22": 3.3e-3, "Yu33": 0.86,
    "Yd11": 1.5e-5, "Yd22": 3e-4, "Yd33": 0.015,
    "Ye11": 2.9e-6, "Ye22": 6e-4, "Ye33": 0.01,
    "M": 1.5e3, "kappa": 0.2, "A": 0.1, "mu": 0.65, "lmbdPhi": 0.3
}

# initialize an instance of the SingletScalarExtension model with the parameter dictionary,
# renormalization scale set to 1000 GeV and loop contrbutions turned on
model1 = SingletScalarExtension(param_dict, 1e3, True)

# evaluate Wilson coefficients as method calls
print(model1.cH())
print(model1.cuH(2,2))

# evaluate coefficients and write to a WCxf file

wcs = ["uu_1331", "G", "phiG", "phiD", "uG_33"]
eft_info = { "eft": "SMEFT", "basis": "Warsaw" }
write_to_wcxf("example_wcxf_all.yaml",model1,eft_info)
