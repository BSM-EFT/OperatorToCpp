This subdirectory contains C++ source files for specific tasks:

1. [pyBindings.cpp](./pyBindings.cpp) contains the necessary definitions to generate a Python module named _match_to_py_ using the pybind11 library. The _match_to_py_ module provides access to the model class and its methods. 

2. [main.cpp](./src/main.cpp) demonstrates how the model class and its methods can be imported and called within C++ source code. One can create instances of the MSSM model with specific parameter values and the corresponding Wilson coefficient functions can be evaluated. 
