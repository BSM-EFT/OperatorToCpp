This is the dedicated working directory for the example case of MSSM-to-SMEFT matching. It comes with the following user-interface files.

1. [OpExp_MSSM.nb](./OpExp_MSSM.nb) - a Mathematica notebook that enables the following steps:
   - Read the EFT information from [SMEFT.m](./input/SMEFT.m).
   - Read the matched expressions as input from [matching-conditions.m](./input/matching-conditions.m).
   - Generate the header and source files `LF.h` and `LF.cpp` that declare and define loop-function.
   - Generate header and source files `MSSM.h` and source `MSSM.cpp` that define the MSSM class (along with a constructor and an updater), incorporate the model parameters as member variables and declare SMEFT Wilson coefficients as methods of the class. The definitions of the SMEFT coefficients (as C++ functions) are stored in separate cXX.cpp files within the [lib](./lib) folder, further segregated based on the number of fermions.
   - Generate the `pyBindings.cpp` file that contains the definitions required to build a Python module named match_to_py using the pybind11 library.
   - Generate the `match_to_py.pyi` file that stores Python-declarations for the MSSM class, the constructor, the updater and the Wilson coefficient methods.

2. [src/main.cpp](./src/main.cpp) - a C++ source file that demonstrates in a simple way, how the model class and its methods can be imported and called within C++ code. 

3. [py/example-MSSM.py](./py/example-MSSM.py) - a simple Python script that demonstrates the basic functionality of the MSSM model class and Wilson coefficients.

5. [py/benchmark-MSSM.py](./py/benchmark-MSSM.py) - a script containing a variety of tasks to record and compare the execution speeds of various functions, e.g. (i) the native yaml writer vs the wcxf writer (ii) execution times of individual Wilson coefficient methods (iii) parallel vs sequential execution of large numbers of coefficients.   

6. [py/plots-MSSM.ipynb](./py/plots-MSSM.ipynb) - a jupyter notebook that reproduces the 2d and bar plots in Figures 2 and 3 of the publication [JHEP 04 (2026) 028](https://link.springer.com/article/10.1007/JHEP04(2026)028).

