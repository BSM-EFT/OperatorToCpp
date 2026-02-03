from operator import call
from time import perf_counter
import numpy as np
import pandas as pd
import yaml
import csv

#====================================================================
# Functions for extracting a callable from a Wilson coefficient name
#====================================================================

def split_wc_name(full_name):
    if "_" in full_name:
        name, idxs = full_name.split("_")
        kw = dict()
        
        for i in range(len(idxs)):
            kw["i" + str(i+1)] = int(list(idxs)[i]) - 1
        
        return (name, kw)

    else:
        return (full_name,)


def eval_wc(model, wc_name, **kw):
    wc_tuple = split_wc_name(wc_name)
    wc = getattr(model,wc_tuple[0])

    if len(wc_tuple) == 2:
        kw.update(wc_tuple[1])
    
    res = call(wc,**kw)
    if abs(res.imag) < 1e-18:
        return res.real
    else:
        return res

#==============================
# Functions to enable File IO
#==============================

# read param information from .yaml files
def read_param_values(filename): 
    with open(filename,'r') as params_file:
        param_dict = yaml.safe_load(params_file)
        ranges_dict = dict()
        keys_to_pop = list()

    for (key, val) in param_dict.items():
        if type(val) == dict:
            param_dict[key] = complex(val["real"], val["imag"])
        elif type(val) == list:
            ranges_dict[key] = np.linspace(val[0],val[1],val[2])
            keys_to_pop.append(key)
    
    for k in keys_to_pop:
        param_dict.pop(k)
    
    return (param_dict, ranges_dict)


# read wilson coefficient names from .txt files
def read_wc_names(filename):
    with open(filename, 'r') as wc_file:
        wc_names = []
        for line in wc_file:
            wc_names.append(line.strip())

    return wc_names


# compute all specified WCs for a fixed benchmark point
def create_wc_dict(model,wc_names,**kw):
    result_dict = dict()
    
    for wc in wc_names:
        result_dict[wc] = eval_wc(model,wc,**kw)

    return result_dict    


# write WC values along with specific parameter values to a .yaml file
def write_to_yaml(filename,model,param_dict,keys,**kw):
    output_dict = dict()
    p_keys = keys[0]
    wc_names = keys[1]

    for k in p_keys:
        output_dict[k] = param_dict[k]

    output_dict.update(create_wc_dict(model,wc_names,**kw))
    
    yaml.add_representer(float,float_representer)
    yaml.add_representer(complex, complex_representer)

    with open(filename, 'w') as out_file:
        yaml.dump(output_dict,out_file,sort_keys=False)


# create a grid based on combinations of parameters that vary between min, max values
def create_combs(arrays, grid=True):
    if grid:
        return np.array(np.meshgrid(*arrays)).T.reshape(-1,len(arrays))
    else:
        combs = []
        for i in range(len(arrays[0])):
            comb = []
            for array in arrays:
                comb.append(array[i])
            combs.append(comb)

        return np.array(combs)

# write WC values for each point on a parameter grid to a .csv file
def write_to_csv(filename,model,ranges_dict,wc_names,grid=True,**kw):
    assert(len(ranges_dict)!=0)
    
    param_keys = list(ranges_dict.keys())
    param_combs = create_combs(list(ranges_dict.values()),grid)

    colNames = list(param_keys) + wc_names

    with open(filename, 'w') as out_file:
        writer = csv.DictWriter(out_file,fieldnames=colNames)
        writer.writeheader()
        
        output_dict = dict()
        for param_comb in param_combs:
            for i in range(len(param_keys)):
                output_dict[param_keys[i]] = param_comb[i]
            
            model.updateParams(output_dict)
            output_dict.update(create_wc_dict(model,wc_names,**kw))
            formatted_output = {
                k: (f"{v.real:.1e}" if isinstance(v, complex) else f"{v:.1e}" if isinstance(v, float) else v) for k, v in output_dict.items()
            }
            writer.writerow(formatted_output)

# create a dataframe containing WC values for each point on a parameter grid
def create_dataframe(model,ranges_dict,wc_names,grid=True,**kw):
    assert(len(ranges_dict)!=0)
    
    param_keys = list(ranges_dict.keys())
    param_combs = create_combs(list(ranges_dict.values()),grid)

    colNames = list(param_keys) + wc_names
    df = pd.DataFrame(columns=colNames)

    output_dict = dict()
    for param_comb in param_combs:
        for i in range(len(param_keys)):
            output_dict[param_keys[i]] = param_comb[i]
        
        model.updateParams(output_dict)
        output_dict.update(create_wc_dict(model,wc_names,**kw))

        df = pd.concat([df, pd.DataFrame(data=[output_dict.values()], columns=list(output_dict.keys()))], ignore_index=True)

    return df 

#========================================================
# Functions for benchmarking the execution time
#========================================================

# decorator for recording the value returned and the execution time of any function that returns
def exec_time(f):
    def wrapper(*args, **kwargs):
        start = perf_counter()
        result = f(*args, **kwargs)
        stop = perf_counter()
        t = stop - start
        return (result, t)

    return wrapper


# decorator for printing the value and execution time of specific functions that admit a model and wc_name as arguments
def print_exec_time(f):
    def wrapper(model, wc_name, **kwargs):
        val, t = f(model, wc_name, **kwargs)
        print(f"{wc_name} = {val:.2e}, execution time = {t:.6} secs.")
    
    return wrapper


@print_exec_time
@exec_time
def exec(model, wc_name, **kwargs):
    return eval_wc(model, wc_name, **kwargs)

#==========================================================================
# Functions for defining float and complex value formatting in yaml output
#==========================================================================

def float_representer(dumper, value):
    expr = f"{value:.1e}"
    return dumper.represent_scalar('tag:yaml.org,2002:float',expr)

def complex_representer(dumper, value):
    return dumper.represent_mapping(
        'tag:yaml.org,2002:map',
        {
            'real': value.real,
            'imag': value.imag
        }
    )