This directory contains the core of OperatorToC++, which is comprised of the following files.

1. [OperatorExport.m](./OperatorExport.m) - a Mathematica package that provides functions to systematically translate Matchete output to C++ classes and methods and generates model-specific .h and .cpp files, as well as a .pyi file containing the corresponding Python declarations of the model class and its methods.

2. [EFTParser.m](./EFTParser.m) - a Mathematica package that provides simple functions to read EFT metadata from either the Matchete models-database or a user-provided file, and store them in JSON format.

3. [CollectLoopFunctions.wl](./CollectLoopFunctions.wl) - a Mathematica script that extracts all unique loop-functions (defined in Matchete's convention), their analytical forms, and degenerate limits and stores them into source (LF.cpp) and header (LF.h) C++ files.

4. [OperatorImport.cpp](./lib/OperatorImport.cpp) - which provides a suite of C++ helper functions that enable computations such as loop-function evaluation, Einstein-summation et cetera within the exported expressions. The corresponding .h file is stored in the [include](./include/) directory.

5. [complex_math.cpp](./lib/complex_math.cpp) - defines overloaded operators to enable arithmetic involving double and complex valued quantities within the method bodies. The corresponding .h file is stored in the [include](./include/) directory. 

6. [utils](./utils/) - a Python module designed to provide file input-output facilities and to interface with popular formats such as csv, yaml, wcxf etc.

To aid in compilation, we provide a template [meson.build](./meson.build) file. It specifies two build targets as 
```    
executable(`main.out', `src/main.cpp', ...)
py.extension\_module(`match\_to\_py', `src/pyBindings.cpp', ...)
```
The specification of the C++ executable should be modified suitably if the source file is **not** named `main.cpp` or if there is more than one source file. In case there is no C++ executable to be built, then the corresponding command should be removed from [meson.build](./meson.build). 

The commands for building and installing the targets are collected in the shell script [compileCpp.sh](./compileCpp.sh). Therefore, simply executing the command, 
```
./compileCpp.sh
```
in a terminal from the working directory (e.g. [Examples/MSSM-to-SMEFTd6](../Examples/MSSM-to-SMEFTd6) ) compiles all source files, links them and the dependencies (`libomp`, `pybind11`) and places the build artifacts in the install directories specified in [meson.build](./meson.build).
