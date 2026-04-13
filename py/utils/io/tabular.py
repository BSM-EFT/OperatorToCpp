from ..core import create_combs, create_wc_dict, generate_call_plan, build_dict, nworkers
from typing import Iterable
import pandas as pd
import csv

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


def create_dataframe_par(model, ranges_dict: dict[str,Iterable], wc_names: list[str], grid=True, **kw: dict) -> pd.DataFrame:
    """Same as create_dataframe() but with multithreading"""

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

        call_plan = generate_call_plan(model, wc_names)
        output_dict.update(build_dict(call_plan,nworkers,**kw))
   
        df = pd.concat([df, pd.DataFrame(data=[output_dict.values()], columns=list(output_dict.keys()))], ignore_index=True)

    return df