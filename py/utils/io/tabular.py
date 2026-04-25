from ..core import create_combs, create_wc_dict, nworkers
from typing import Iterable
import pandas as pd
import csv
from concurrent.futures import ThreadPoolExecutor

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
            Fixed global parameters such as "mubarsq" for renormalization scale and
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



def create_row(p_comb, p_keys: list[str], wc_names, model, fixed_pars, **kw) -> dict[str, float | complex]:
    """Helper function for creating a single row with parameter and Wilson coefficient values"""

    row_model = model(fixed_pars)
    row_dict = dict()
    for i in range(len(p_keys)):
        row_dict[p_keys[i]] = p_comb[i]

    row_model.updateParams(row_dict)
    row_dict.update(create_wc_dict(row_model, wc_names, **kw))

    return row_dict
    

def create_dataframe(model, fixed_pars: dict[str, float|complex], ranges_dict: dict[str,Iterable], wc_names: list[str], grid:bool = True, max_workers:int = nworkers, **kw: dict):
    """
    Generates combinations for parameters that vary over specified ranges, evaluates the
    corresponding Wilson coefficents (for the UV model) and returns a pandas dataframe 
    object containing the parameter and Wilson coefficient names as columns.
    
    Parameters
    ----------
        model
            The UV model class defined in the match_to_py module
        fixed_pars : dict
            Dictionary containing (name, value) pairs for fixed model parameters
        ranges_dict : dict
            Dictionary containing (name, range) pairs for model parameters that
            vary over specific ranges
        wc_names : list
            List of strings describing Wilson coefficient names
        grid : bool, optional
            A flag used to specify whether the combinations should be created
            as Cartesian products (the default option) or by simply combining 
            the corresponding entries of each array within the arrays list
        max_workers : int, optional
            Number of cores over which the multi-threaded task will be distributed
        kw : dict
            Fixed global parameters such as "mubarsq" for renormalization scale and
            "hbar" for the matching order, specified using a dictionary

    Returns
    --------
        A pandas.Dataframe object containing parameter combinations and the corresponding
        Wilson coefficient values 
    
    """
    
    assert(len(ranges_dict)!=0)
    
    param_keys = list(ranges_dict.keys())
    param_combs = create_combs(list(ranges_dict.values()),grid)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(
            create_row, p_comb, param_keys, wc_names, model, fixed_pars, **kw
            ) for p_comb in param_combs 
        ]

        results = [f.result() for f in futures]

    return pd.DataFrame(results)
