This subdirectory contains the Mathematica component of OperatorToC++, which is comprised of the following files.

1.  [OperatorExport.m](./OperatorExport.m) - a Mathematica package that provides functions to systematically translate Matchete output to C++ classes and methods and generates model-specific .h and .cpp files, as well as a .pyi file containing the corresponding Python declarations of the model class and its methods.

2.  [CollectLoopFunctions.wl](./CollectLoopFunctions.wl) - a Mathematica script that extracts all unique loop-functions (defined in Matchete's convention), their analytical forms, and degenerate limits and stores them into source (LF.cpp) and header (LF.h) C++ files.

3. [OpExp_MSSM.nb](./OpExp_MSSM.nb) - is the model-specific notebook that enables the following steps:
   - Read the matched expressions as input from [MSSM-matching-conditions.m](./matching-results/MSSM-matching-conditions.m).
   - Execute the [CollectLoopFunctions.wl](./CollectLoopFunctions.wl) script and generate the files [LF.h](../include/LF.h) and [LF.cpp](../lib/LF.cpp).
   - Apply necessary modifications to the matching output (e.g. replacing Greek letters by a more C++ friendly form).
   - Generate header [MSSM.h](../include/MSSM.h) and source [MSSM.cpp](../lib/MSSM.cpp) files that define the MSSM class (along with a constructor and an updater), incorporate the model parameters as member variables and declare SMEFT Wilson coefficients as methods of the class. The definitions of the SMEFT coefficients (as C++ functions) are stored in separate cXX.cpp files within the [lib](../lib/) folder, further segregated based on the number of fermions.
   - Generate the [match_to_py.pyi](../py/match_to_py.pyi) file that stores Python-declarations for the MSSM class, the constructor, the updater and the Wilson coefficient methods.
