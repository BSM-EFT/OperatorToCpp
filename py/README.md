This subdirectory provides a Python-based frontend for importing and calling the UV model class and the Wilson coefficient methods, within Python based HEP workflows.

It contains the following scripts and modules:

1. match_to_py.<...>.so - This is a shared object (.so) file not available in the repository by default but generated after compiling and building the match_to_py Python module by running the [run.sh](../run.sh) script. This module contains the UV model class and Wilson coefficient methods. 

2. [match_to_py.pyi](./match_to_py.pyi) -  This file is generated from the Mathematica frontend [OpExp_MSSM.nb](../mathematica/OpExp_MSSM.nb) and it stores Python-declarations for the UV model class, the constructor, the updater and the Wilson coefficient methods.

3. [utils](./utils/) - This is a Python module designed to provide file input-output facilities and to interface with popular formats such as csv, yaml, wcxf etc.

4. [EFT_db](./EFT_db/) - This subdirectory functions as a database for the various operator bases corresponding to different EFTs. The information about specific Wilson coefficients is stored in WCInfo.json files, e.g. [SMEFT/Warsaw/WCInfo.json](./EFT_db/SMEFT/Warsaw/WCInfo.json). The JSON file is generated from the Mathematica frontend [OpExp_MSSM.nb](../mathematica/OpExp_MSSM.nb). 

4. [example.py](./example.py) - This is a simple script tht demonstrates the basic functionality of the UV model class and Wilson coefficients provided by the match_to_py module.

5. [benchmark-MSSM.py](./benchmark-MSSM.py) - This script contains a variety of tasks to record and compare the execution speeds of various functions, e.g. (i) the native yaml writer vs the wcxf writer (ii) execution times of individual Wilson coefficient methods (iii) parallel vs sequential execution of large numbers of coefficients. This is also useful to study how the execution speeds differ across different compuatational platforms.  

6. [plots_mssm.ipynb](./plots_mssm.ipynb) - This jupyter notebook contains the code used to create the 2d and bar plots in Figures 2 and 3 of the publication [JHEP 04 (2026) 028](https://link.springer.com/article/10.1007/JHEP04(2026)028).
