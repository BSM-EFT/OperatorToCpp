from match_to_py import MSSM
from utils.io import read_param_values, write_to_wcxf_all
from time import perf_counter

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
print(f"Data for {len(param_dicts)} benchmark points loaded in {t2 - t1}s.")

print("Initializing MSSM model instances.")
t3 = perf_counter()
models = list()
for p_dict in param_dicts:
    models.append(MSSM(p_dict))

t4 = perf_counter()
print(f"Initialized {len(param_dicts)} models in {t4 - t3}s.")

scale = {"mubarsq": 1e6, "hbar": 0.006332574}

print("Writing to files.")

tf0 = perf_counter()
write_to_wcxf_all("./benchmark-points/wcxf_all_1.yaml",models[0],**scale)
tf1 = perf_counter()
print(f"File 1 written in {tf1 - tf0}s.")

write_to_wcxf_all("./benchmark-points/wcxf_all_2.yaml",models[1],**scale)
tf2 = perf_counter()
print(f"File 2 written in {tf2 - tf1}s.")

write_to_wcxf_all("./benchmark-points/wcxf_all_3.yaml",models[2],**scale)
tf3 = perf_counter()
print(f"File 3 written in {tf3 - tf2}s.")

write_to_wcxf_all("./benchmark-points/wcxf_all_4.yaml",models[3],**scale)
tf4 = perf_counter()
print(f"File 4 written in {tf4 - tf3}s.")

write_to_wcxf_all("./benchmark-points/wcxf_all_5.yaml",models[4],**scale)
tf5 = perf_counter()
print(f"File 5 written in {tf5 - tf4}s.")

print(f"Total time taken - {tf5 - t1}s")
# ideas for improving the speed for write_all --> maybe the list comprehension where we call eval_wc multiple times is slowing us down? perhaps we should create numpy arrays and store them for reuse, instead of generating the dynamically