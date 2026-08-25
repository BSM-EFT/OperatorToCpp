from ..core import create_wc_dict, create_tasks
import numpy as np
import yaml
from typing import Iterable


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


def write_to_yaml(filename:str, model, param_dict: dict[str, float | complex], keys: list[list[str], list[str]], opt: str = "par") -> None:
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
        opt : str, "seq" | "par"
            Option for specifying whether the operation should proceed sequentially, or 
            in parallel with OpenMP at the C++ backend 
        kw : dict
            Fixed global parameters such as "mubarsq" for renormalization scale and
            "hbar" for the matching order, specified using a dictionary
    
    """
    match opt:
        case "seq":
            write_to_yaml_seq(filename, model, param_dict, keys)
        case "par":
            write_to_yaml_omp(filename, model, param_dict, keys)


def write_to_yaml_seq(filename:str, model, param_dict: dict[str, float | complex], keys: list[list[str], list[str]]) -> None:
    """
    Writes the values of specified parameters and evaluated Wilson coefficients
    to a .yaml file, follows sequential execution during dictionary building
    """
    output_dict = dict()
    p_keys = keys[0]
    wc_names = keys[1]

    for k in p_keys:
        output_dict[k] = param_dict[k]

    output_dict.update(create_wc_dict(model,wc_names))

    zero_cutoff = 1e-20
    new_dict = {k: v for k, v in output_dict.items() if abs(v) > zero_cutoff}    
    
    yaml.add_representer(complex, complex_representer)
    with open(filename, 'w') as out_file:
        yaml.dump(new_dict,out_file,sort_keys=False)


def write_to_yaml_omp(filename:str, model, param_dict: dict[str, float | complex], keys: list[list[str], list[str]]) -> None:
    """
    Writes the values of specified parameters and evaluated Wilson coefficients
    to a .yaml file, utilizes parallelisation through OpenMP and batch processing at
    the C++ end
    """
    
    output_dict = dict()
    p_keys = keys[0]
    wc_names = keys[1]

    for k in p_keys:
        output_dict[k] = param_dict[k]

    tasks = create_tasks(model, wc_names)
    output_dict.update(model.batch_eval(tasks))

    zero_cutoff = 1e-20
    new_dict = {k: v for k, v in output_dict.items() if abs(v) > zero_cutoff}    
    
    yaml.add_representer(complex, complex_representer)
    with open(filename, 'w') as out_file:
        yaml.dump(new_dict,out_file,sort_keys=False)


def float_representer(dumper, value: float):
    """Format specifier for floating-point output when written to a .yaml file"""

    expr = f"{value:.1e}"
    return dumper.represent_scalar('tag:yaml.org,2002:float',expr)


def complex_representer(dumper, value: complex):
    """Format specifier for complex output when written to a .yaml file"""
    
    return dumper.represent_mapping(
        'tag:yaml.org,2002:map',
        {                
            'Re': value.real,
            'Im': value.imag
        }
    )
