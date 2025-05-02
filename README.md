# OperatorToC++

OperatorToC++ is an extensible, hybrid Mathematica and C++ based tool that facilitates next steps beyond the matching of parameters between a UV model and an Effective Field Theory. It efficiently remedies the complexities within the analytical matched expressions such as those due to loop-functions and lengthy sums and products involving tensor objects. It then translates and bundles the results into C++ classes and functions which provide a convenient platform for further numerical analyses.

Its two main components are  
1.  [OperatorExport.m](./OperatorExport.m) - a Mathematica package that systematically translates Matchete output to C++ classes and functions and generates model-specific .h and .cpp files.  
2.  [OperatorImport.cpp](./lib/OperatorImport.cpp) - a suite of C++ helper functions that enable computations such as loop-function evaluation, Einstein-summation et cetera.

## Sample input and output files

We demonstrate the basic functionality and ease of use of OperatorToC++ using the highly non-trivial case of matching the general MSSM with the dimension-6 SMEFT coefficients. The physics details regarding the matching procedure and the subtleties involved were outlined in the article [SUSY meets SMEFT](https://www.arxiv.org). 
 
The notebook [OpExp_MSSM.nb](./OpExp_MSSM.nb) which reads the matched expressions as input from [MSSM-matching-conditions.m](./MSSM-matching-conditions.m) and generates header and [MSSM.h](./include/MSSM.h) source [MSSM.cpp](./lib/MSSM.cpp) files. These files define an MSSM class, incorporate the gauge and Yukawa couplings and masses as member variables and the SMEFT Wilson coefficients as methods of the class.

We also provide an easy way to reproduce the plots from [SUSY meets SMEFT](https://www.arxiv.org). One can use the [makefile](./makefile) to compile the short C++ program [write_to_files.cpp](./src/write_to_files.cpp) to generate .txt files containing numerical values of the Wilson coefficients for specific parameter choices. These .txt files are then read as inputs by [plots.ipynb](./plots/plots.ipynb) which then generates the 2d and bar-plots from Figs.1 and 2 of [SUSY meets SMEFT](https://www.arxiv.org).

## Instructions for use beyond the MSSM matching

OperatorToC++ is useful beyond just cross-checking and reproducing the numerical results from [SUSY meets SMEFT](https://www.arxiv.org). In principle, one can clone this repository and suitably integrate [OperatorExport.m](./OperatorExport.m), [OperatorImport.cpp](./lib/OperatorImport.cpp), and [OperatorImport.h](./include/OperatorImport.h) within their wokflow. However, the *recommended option* is to download and use a suitable release version.

OperatorToC++, in its present iteration, is fully self-contained and does not rely on third-party libraries, either on the Mathematica end or on the C++ end. However, we use several modern C++ features and *the code requires a compiler that abides by the C++23 standard*.

In order to use OperatorToC++ for transpiling the matching conditions relating an arbitrary UV model with SMEFT coefficients, 
 - One requires the Matchete output stored (as a key-value pair with Warsaw basis Wilson coefficients as keys) in a .m file which can be read from a notebook such as OpExp_MSSM.nb (after incorporating model-specific changes to the notebook). 
 - One can then call the relevant functions to generate <model>.h and <model>.cpp files, which by default are placed inside the **include** and **lib** directories.
 - Any additional .cpp source files can be stored in the **src** directory and one must extend the makefile to account for such files.
 - If the directory structure is updated/modified, then the makefile should be suitably modified as well.
