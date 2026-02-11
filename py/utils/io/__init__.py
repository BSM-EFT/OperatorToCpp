from .tabular import write_to_csv, create_dataframe
from .yaml import read_param_values, write_to_yaml


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
