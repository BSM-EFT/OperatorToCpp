# OperatorToC++

[![DOI](https://zenodo.org/badge/948956654.svg)](https://doi.org/10.5281/zenodo.15599937)

OperatorToC++ is an extensible, hybrid Mathematica and C++ based tool that facilitates next steps beyond the matching of parameters between a UV model and an Effective Field Theory. It efficiently remedies the complexities within the analytical matched expressions such as those due to loop-functions and lengthy sums and products involving tensor objects. It then translates and bundles the results into C++ classes and methods which provide a convenient platform for further numerical analyses. 

Additionally, it also provides a convenient Python front-end from where the transpiled classes and functions can be called and incorporated within Python-based HEP workflows. 

It admits [Matchete](https://gitlab.com/matchete/matchete/) (>= v0.3.0) output, stored in a <matching-results>.m file in a key-value pair format, where the keys correspond to Standard Model EFT coefficients and the values are the respective matching expressions, i.e. functions of the parameters of the specific UV model.

The backbone of OperatorToC++ is constituted of 
1.  [OperatorExport.m](./mathematica/OperatorExport.m) - a Mathematica package that provides functions to systematically translate Matchete output to C++ classes and methods and generates model-specific .h and .cpp files, as well as a .pyi file containing the corresponding Python declarations of the model class and its methods.
2.  [CollectLoopFunctions.wl](./mathematica/CollectLoopFunctions.wl) - a Mathematica script that extracts all unique loop-functions (defined in Matchete's convention), their analytical forms, and degenerate limits and stores them into source (LF.cpp) and header (LF.h) C++ files.
3.  [OperatorImport.cpp](./lib/OperatorImport.cpp) - a suite of C++ helper functions that enable computations such as loop-function evaluation, Einstein-summation et cetera.

Additionally, the source file [pyBindings.cpp](./src/pyBindings.cpp) enables the creation of a Python module with the C++ class and methods. We have also added IO facilities at the Python-end by bundling some helper functions within [utilities.py](./py/utilities.py).

### Prerequisites and Dependecies

OperatorToC++ requires a compiler that abides by the C++17 standard. The [Matchete](https://gitlab.com/matchete/matchete/) package is needed for the extraction of loop functions. It also needs OpenMP to parallelise some operations defined within OperatorImport.cpp and [pybind11](https://pybind11.readthedocs.io/en/stable/installing.html) to generate the Python module. To simplify the build process of the C++ code and to generate the Python module [cmake](https://cmake.org/) is required, We provide a [CMakeLists.txt](./CMakeLists.txt) and an [install.sh](./install.sh) file to automate the build process.

_The CMakeLists.txt file provides a manual fix for locating the OpenMP files on MacOS. This may need to be updated based on the specific libomp version installed._

The user can interact with the code at 3 distinct levels as described in the next section.

## An Example - Exporting the matching results for the MSSM

The basic functionality and ease of use of OperatorToC++ can be demonstrated using the highly non-trivial case of matching the general MSSM with the dimension-6 SMEFT coefficients. 

We organize the Mathematica and Python codes in separate folders: [mathematica](./mathematica/) and [py](./py/). 

There are 3 levels at which the user can interact with the code and the matching results.

1. The Mathematica frontend: The notebook [OpExp_MSSM.nb](./mathematica/OpExp_MSSM.nb) enables the following steps:
   - Read the matched expressions as input from [MSSM-matching-conditions.m](./matching/MSSM-matching-conditions.m).
   - Execute the [CollectLoopFunctions.wl](./mathematica/CollectLoopFunctions.wl) script and generate the files [LF.h](./include/LF.h) and [LF.cpp](./lib/LF.cpp).
   - Apply necessary modifications to the matching output (e.g. replacing Greek letters by a more C++ friendly form).
   - Generate header [MSSM.h](./include/MSSM.h) and source [MSSM.cpp](./lib/MSSM.cpp) files that define the MSSM class (along with a constructor and an updater), incorporate the model parameters as member variables and declare SMEFT Wilson coefficients as methods of the class. The definitions of the SMEFT coefficients (as C++ functions) are stored in separate cXX.cpp files within the [lib](./lib/) folder, further segregated based on the number of fermions.
   - Generate the [match_to_py.pyi](./py/match_to_py.pyi) file that stores Python-declarations for the MSSM class, the constructor, the updater and the Wilson coefficient methods.
 
2. The C++ frontend and the build system:
   - Using a file like [main.cpp](./src/main.cpp), one can create instances of the MSSM model with specific parameter values and the corresponding Wilson coefficient functions can be evaluated. Currently, [main.cpp](./src/main.cpp) contains some simple code that benchmarks the speed of execution for the cH() function (which was found to be the most time-consuming among all Wilson coefficients).
   - Additionally, the source file [pyBindings.cpp](./src/pyBindings.cpp) contains the necessary definitions to generate a Python module named _match_to_py_ using the pybind11 library.
   - To be able to compile and link main.cpp (as well as pyBindings.cpp) with all the MSSM files, LF.cpp and OperatorImport.cpp, the [CMakeLists.txt](./CMakeLists.txt) file can be used. However, the build process is easily simplified through the [install.sh](./install.sh) script
  
4. The Python frontend: After running install.sh, a Python shared object (match_to_py_<...>.so) file is created within the [py](./py/) subdirectory. Python functions that enable reading of parameter values from .yaml files and writing evaluated Wilson coefficient values to .yaml or .csv files are provided in [utilities.py](./py/utilities.py). Using a notebook such as [examples.ipynb](./py/examples.ipynb), one can create instances of the MSSM class, evaluate Wilson coefficients for fixed as well ranges of parameter values, create plots or incorporate the MSSM class and its methods into further numerical analyses.

By default all model parameters are internally defined as complex-valued. However, by declaring which parameters are actually complex before generating the C++ files in [OpExp_MSSM.nb](./mathematica/OpExp_MSSM.nb), we tell the code about parameters for which distinct complex conjugates should be defined. 

All Wilson coefficient functions return complex-valued output and the real and imaginary parts can be extracted suitably in C++ as well as Python.

This repository makes all C++ files generated after evaluating [OpExp_MSSM.nb](./mathematica/OpExp_MSSM.nb) already available. 

## Instructions for use beyond the MSSM matching

While originally developed to help generate numerical output from the MSSM to SMEFT matching results, over time OperatorToC++ has been made entirely model-independent. It is easy to use the code to generate C++ (and Python) classes and methods for any UV model matched onto the dimension-6 SMEFT (in the Warsaw basis).

One simply needs to replace the [MSSM-matching-conditions.m](./matching/MSSM-matching-conditions.m) file by the particular <...>-matching-conditions.m file and make sure to update the model-specific steps in a OpExp_<...>.nb file, e.g. the renaming of variables, identification of complex parameters and declaring the model name. This should enable one to create <...>.h and <...>.cpp files and update the definitions of the Wilson coefficient (cXX.cpp) and the match_to_py.pyi files, as well as the contents of LF.cpp. 

The only hardcoded part in pyBindings.cpp and CMakeLists.txt is the name of the Python module _match_to_py_, the name of the class, i.e. MSSM or <...> is introduced internally in such a way that one wouldn't have to copy-paste the "Model-Name" in multiple places throughout the code.

We recommend maintaining the same directory structure as in this repository. However, if one updates it or adds new source files to the [src](./src/) directory then the CMakeLists.txt should be suitably updated.

## Citation

If you use this code, please cite:

> Suraj Prakash, *OperatorToC++*. Zenodo, 2025. DOI: [10.5281/zenodo.15599937](https://doi.org/10.5281/zenodo.15599937)

