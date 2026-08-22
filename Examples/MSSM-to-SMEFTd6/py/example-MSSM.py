from match_to_py import MSSM
from utils.io import write_to_wcxf
from utils.core import eval_wc
from matplotlib import pyplot as plt
import numpy as np

# define a parameter dictionary

param_dict = {
    "g1": 0.37, "g3": 1.1, "cgamma": 0.01, "Yu11": 0.00001, "Yu22": 0.007, "Yu33": 0.9,
    "met1": 1e9, "met2": 1e9, "met3": 1e9, "mlt1": 1e9, "mlt2": 1e9, "mlt3": 1e9,
    "mqt1": 1e9, "mqt2": 1e9, "mqt3": 1e9, "mut1": 1e9, "mut2": 1e9, "mut3": 2e3,
    "mdt1": 1e9, "mdt2": 1e9, "mdt3": 1e9, "m1": 1.2e3, "m2": 1e9, "m3": 1e9,
    "mPhi": 1e9, "muTilde": 1e9
}

# initialize an instance of the MSSM model with the parameter dictionary,
# renormalization scale set to 1000 GeV and loop contrbutions turned on
model1 = MSSM(param_dict, 1e3, True)

# evaluate Wilson coefficients as method calls
print(model1.cG())
print(model1.cuB(2,2))

# evaluate coefficients and write to a WCxf file

wcs = ["uu_1331", "G", "phiG", "phiD", "uG_33"]
eft_info = { "eft": "SMEFT", "basis": "Warsaw" }
write_to_wcxf("example_wcxf.yaml",model1,eft_info,wcs,opt="seq")
write_to_wcxf("example_wcxf_all.yaml",model1,eft_info)

# obtain arrays of Wilson coefficients for varying m1, keeping other parameters fixe

m1_range = np.linspace(1200,2700,16)
cqq1_vals, cuG_vals, cHq1_vals = [], [], []
for m in m1_range:
    d = {"m1": m}
    model1.updateParams(d)
    cqq1_vals.append(eval_wc(model1, "cqq1_3333"))
    cuG_vals.append(eval_wc(model1, "cuG_33"))
    cHq1_vals.append(eval_wc(model1, "cHq1_33"))

# Set plot attributes
plt.rcParams['axes.labelpad'] = 12
plt.rcParams['xtick.labelsize'] = 20
plt.rcParams['ytick.labelsize'] = 20
plt.rcParams['axes.labelsize'] = 20
plt.rcParams['legend.fontsize'] = 16
plt.rcParams["text.usetex"] = True
plt.rcParams.update({"savefig.dpi" : 300})
plt.rcParams['text.latex.preamble'] = r'\usepackage{amssymb}'

# Create the plot 

plt.figure(figsize=(8,8))
plt.plot(m1_range/1e3, np.abs(cqq1_vals)*1e12, color="red", label=r"$|C_{qq}^{(1),3333}|$")
plt.plot(m1_range/1e3, np.abs(cuG_vals)*1e12, color="blue", label=r"$|C_{uG}^{33}|$")
plt.plot(m1_range/1e3, np.abs(cHq1_vals)*1e12, color="green", label=r"$|C_{Hq}^{(1),33}|$")
plt.xticks([1.2,1.5,1.8,2.1,2.4,2.7])
plt.yticks([1,2,3,4,5,6])
plt.xlim(1.2,2.7)
plt.ylim(1,6)
plt.xlabel(r"$m_1\,\,[{\rm TeV}]$")
plt.ylabel(r"$10^{-12}\,\,[{\rm GeV}^{-2}]$")
plt.legend()
plt.savefig("lineplot.pdf")
