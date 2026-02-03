from operator import call
from time import perf_counter
from typing import Callable, Iterable
import numpy as np
import pandas as pd
import yaml
import csv

#====================================================================
# Functions for extracting a callable from a Wilson coefficient name
#====================================================================

def split_wc_name(full_name: str) -> tuple[str,dict[str, int]]:
    """
    Splits the name and index information from a Wilson coefficient string
    
    Parameters
    ----------
        full_name : str
            A string containing the coefficient name in "cXYZ_abcd" format 
    
    Returns
    -------
        tuple[str, dict]
            A tuple containing the coefficient name as a string "cXYZ" and
            index information as a dictionary { "i1": a, "i2": b, ... }
    """
    if "_" in full_name:
        name, idxs = full_name.split("_")
        kw = dict()
        
        for i in range(len(idxs)):
            kw["i" + str(i+1)] = int(list(idxs)[i]) - 1
        
        return (name, kw)

    else:
        return (full_name,)


def eval_wc(model, wc_name: str, **kw: dict) -> float | complex:
    """
    Evaluates the specified Wilson coeffcient method for a given model
    
    Parameters
    ----------
        model : 
            Instance of the UV model class defined in the match_to_py module
        wc_name : str
            A string containing the coefficient name in "cXYZ_abcd" format
        kw : dict
            Additional parameters such as "mubarsq" for renormalization scale and
            "hbar" for the matching order, specified using a dictionary

    Returns
    -------
        float | complex
            Result of evaluating the Wilson coefficient method. Generally complex 
            but only the real part is returned if the imaginary part is 0. 

    """
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
def read_param_values(filename: str) -> tuple[dict[str, float | complex], dict[str, Iterable]]:
    """
    Reads the parameter values from a .yaml file and stores them into 
    relevant dictionaries

    Parameters
    ----------
        filename : str
            A string specifying the path of the input .yaml file

    Returns
    --------
        tuple[dict, dict]
            A tuple containing (i) a dictionary of (name, value) pairs
            for fixed-value parameters and (ii) a dictionary of (name, range)
            pairs for parameters that vary over a range specified using
            [min, max, step] values in the input file.     
    """
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
def read_wc_names(filename: str) -> list[str]:
    """
    Reads the Wilson coefficient names from a .txt file and returns 
    a list of strings
    """
    with open(filename, 'r') as wc_file:
        wc_names = []
        for line in wc_file:
            wc_names.append(line.strip())

    return wc_names


# compute all specified WCs for a fixed benchmark point
def create_wc_dict(model, wc_names: list[str], **kw: dict) -> dict[str, float | complex]:
    """
    Creates a dictionary to store (name, value) pairs for the specified
    list of Wilson coefficients for a given model
    """

    result_dict = dict()
    
    for wc in wc_names:
        result_dict[wc] = eval_wc(model,wc,**kw)

    return result_dict    


# write WC values along with specific parameter values to a .yaml file
def write_to_yaml(filename:str, model, param_dict: dict[str, float | complex], keys: list[list[str], list[str]], **kw: dict) -> None:
    """
    Writes the values of specified parameters and evaluated Wilson coefficients
    to a .yaml file

    Parameters
    ----------
        filename : str
            A string specifying the path of the output .yaml file
        model
            Instance of the UV model class defined in the match_to_py module
        param_dict : dict
            Dictionary containing (name, value) pairs of model parameters
        keys : list[list, list]
            List containing a list of parameter names and a list of coefficient names
        kw : dict
            Additional parameters such as "mubarsq" for renormalization scale and
            "hbar" for the matching order, specified using a dictionary
    
    """
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
def create_combs(arrays: list[list[float | complex]], grid=True) -> np.ndarray[np.ndarray]:
    """
    Creates combinations of parameter values from individual lists

    Parameters
    ----------
        arrays : list
            A list of lists containing values of individual parameters
        grid : bool, optional
            A flag used to specify whether the combinations should be created
            as Cartesian products (the default option) or by simply combining 
            the corresponding entries of each array within the arrays list 

    Returns
    --------
        A 2d numpy array containing the generated parameter combinations    
    """
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
def write_to_csv(filename: str, model, ranges_dict: dict[str,Iterable], wc_names: list[str], grid=True, **kw: dict) -> None:
    """
    Writes the values of parameters that vary over a range and the evaluated 
    Wilson coefficients (for each parameter combination) to a .csv file

    Parameters
    ----------
        filename : str
            A string specifying the path of the output .yaml file
        model
            Instance of the UV model class defined in the match_to_py module
        ranges_dict : dict
            Dictionary containing (name, range) pairs for model parameters that
            vary over specific ranges
        wc_names : list
            List of strings describing Wilson coefficient names
        grid : bool, optional
            A flag used to specify whether the combinations should be created
            as Cartesian products (the default option) or by simply combining 
            the corresponding entries of each array within the arrays list 
        kw : dict
            Additional parameters such as "mubarsq" for renormalization scale and
            "hbar" for the matching order, specified using a dictionary
    """
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
def create_dataframe(model, ranges_dict: dict[str,Iterable], wc_names: list[str], grid=True, **kw: dict) -> pd.DataFrame:
    """
    Generates combinations for parameters that vary over specified ranges, evaluates the
    corresponding Wilson coefficents (for the UV model) and returns a pandas dataframe 
    object containing the parameter and Wilson coefficient names as columns.    
    
    Parameters
    ----------
        model
            Instance of the UV model class defined in the match_to_py module
        ranges_dict : dict
            Dictionary containing (name, range) pairs for model parameters that
            vary over specific ranges
        wc_names : list
            List of strings describing Wilson coefficient names
        grid : bool, optional
            A flag used to specify whether the combinations should be created
            as Cartesian products (the default option) or by simply combining 
            the corresponding entries of each array within the arrays list 
        kw : dict
            Additional parameters such as "mubarsq" for renormalization scale and
            "hbar" for the matching order, specified using a dictionary

    Returns
    --------
        A pandas.Dataframe object containing parameter combinations and the corresponding
        Wilson coefficient values 
    
    """
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

def exec_time(f: Callable) -> Callable:
    """
    A decorator for recording the value and execution time for a function
    """
    def wrapper(*args, **kwargs):
        start = perf_counter()
        result = f(*args, **kwargs)
        stop = perf_counter()
        t = stop - start
        return (result, t)

    return wrapper


def print_exec_time(f: Callable) -> Callable:
    """
    A decorator for printing the value and execution time of functions 
    that admit a model and wc_name as arguments
    """
    def wrapper(model, wc_name, **kwargs):
        val, t = f(model, wc_name, **kwargs)
        print(f"{wc_name} = {val:.2e}, execution time = {t:.6} secs.")
    
    return wrapper


@print_exec_time
@exec_time
def exec(model, wc_name: str, **kwargs: dict) -> float | complex:
    """
    Decorated version of eval_wc() that prints the output and the execution time 
    """
    return eval_wc(model, wc_name, **kwargs)

#==========================================================================
# Functions for defining float and complex value formatting in yaml output
#==========================================================================

def float_representer(dumper, value: float):
    """Format specifier for floating-point output when written to a .yaml file"""

    expr = f"{value:.1e}"
    return dumper.represent_scalar('tag:yaml.org,2002:float',expr)

def complex_representer(dumper, value: complex):
    """Format specifier for complex output when written to a .yaml file"""

    return dumper.represent_mapping(
        'tag:yaml.org,2002:map',
        {
            'real': value.real,
            'imag': value.imag
        }
    )