This is the dedicated working directory for the simple case of the matching between a Type I Seesaw Model and SMEFT. It comes with the following user-interface files.

1. [OpExp_SeesawTypeI.nb](OpExp_SeesawTypeI.nb) - a Mathematica notebook that enables the following steps:
   - Read the EFT information from [SMEFT.m](./input/SMEFT.m).
   - Read the matched expressions as input from [matching-conditions.m](./input/matching-conditions.m).
   - Generate the header and source files `LF.h` and `LF.cpp` that declare and define loop-function.
   - Generate header and source files `SeesawTypeI.h` and source `SeesawTypeI.cpp` that define the SeesawTypeI class (along with a constructor and an updater), incorporate the model parameters as member variables and declare SMEFT Wilson coefficients as methods of the class. The definitions of the SMEFT coefficients (as C++ functions) are stored in separate cXX.cpp files within the [lib](./lib) folder, further segregated based on the number of fermions.
   - Generate the `pyBindings.cpp` file that contains the definitions required to build a Python module named match_to_py using the pybind11 library.
   - Generate the `match_to_py.pyi` file that stores Python-declarations for the SeesawTypeI class, the constructor, the updater and the Wilson coefficient methods.

2. [src/main.cpp](./src/main.cpp) - a C++ source file that demonstrates in a simple way, how the model class and its methods can be imported and called within C++ code. 

3. [py/example-SeesawTypeI.py](./py/example-SeesawTypeI.py) - a simple Python script that demonstrates the basic functionality of the SeesawTypeI model class and Wilson coefficients.

5. [py/benchmark-SeesawTypeI.py](./py/benchmark-SeesawTypeI.py) - a script containing a variety of tasks to record and compare the execution speeds of various functions, e.g. (i) the native yaml writer vs the wcxf writer (ii) execution times of individual Wilson coefficient methods (iii) parallel vs sequential execution of large numbers of coefficients.  

We also provide a Mathematica notebook [SeesawTypeI-matching.nb](./input/SeesawTypeI-matching.nb), with the Matchete code to obtain the matching relations, in the input subfolder.
