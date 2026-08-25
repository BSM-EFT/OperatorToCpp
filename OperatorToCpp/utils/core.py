from operator import call
from time import perf_counter
from typing import Callable
import numpy as np
from concurrent.futures import ThreadPoolExecutor

import subprocess
import platform

if platform.system() == "Darwin":
    cmd = ['sysctl', '-n', 'hw.perflevel0.logicalcpu']
    nworkers = int(subprocess.check_output(cmd))
elif platform.system() == "Linux":
    cmd = ['nproc']
    nworkers = int(subprocess.check_output(cmd))

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


def eval_wc(model, wc_name: str) -> complex:
    """
    Evaluates the specified Wilson coeffcient method for a given model
    
    Parameters
    ----------
        model : 
            Instance of the UV model class defined in the match_to_py module
        wc_name : str
            A string containing the coefficient name in "cXYZ_abcd" format

    Returns
    -------
        complex
            Result of evaluating the Wilson coefficient method. 

    """
    wc_tuple = split_wc_name(wc_name)
    wc = getattr(model,wc_tuple[0])

    kw = dict()
    if len(wc_tuple) == 2:
        kw = {**kw,**wc_tuple[1]}
    
    res = call(wc,**kw)
    return res


def create_wc_dict(model, wc_names: list[str]) -> dict[str, float | complex]:
    """
    Creates a dictionary to store (name, value) pairs for the specified
    list of Wilson coefficients for a given model
    """

    result_dict = dict()
    
    for wc in wc_names:
        result_dict[wc] = eval_wc(model,wc)

    return result_dict    


def generate_call_plan(model, wc_names: list[str]) -> list[dict]:
    plan = []
    for wc_name in wc_names:
        parsed = split_wc_name(wc_name)
        method_name = parsed[0]
        method = getattr(model, method_name)

        indices = parsed[1] if len(parsed) > 1 else {}
        weight = 3 if not indices else (2 if len(indices) == 2 else 1)

        plan.append({
            "weight": weight,
            "name": wc_name,
            "method": method,
            "kwargs": indices
        })

    plan.sort(key=lambda x: x['weight'], reverse=True)
    return plan


def build_dict(plan, num_workers) -> dict[str, float | complex]:
    result_dict: dict[str, float | complex] = {}

    def worker(task) -> tuple[str, float | complex]:
        return task['name'], task['method'](**task['kwargs'])

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        for name, val in executor.map(worker, plan, chunksize=20):
            result_dict[name] = val

    return result_dict


def create_tasks(model, names):
    plan = []
    for name in names:
        wc_tuple = split_wc_name(name)
        method_name = wc_tuple[0]

        if len(wc_tuple) == 1:
            task = model.wrap_0f(name, method_name)
        elif len(wc_tuple) == 2:
            if len(wc_tuple[1]) == 2:
                task = model.wrap_2f(name, method_name, **wc_tuple[1])
            elif len(wc_tuple[1]) == 4:
                task = model.wrap_4f(name, method_name, **wc_tuple[1])

        plan.append(task)
    return plan


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
