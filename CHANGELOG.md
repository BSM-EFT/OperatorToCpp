## Since v0.1

### Added

1. Automated parsing of EFT / Wilson Coefficient metadata. This extends the utility of the code beyond the SMEFT-Warsaw basis.
2. A new Python front-end and the interop with Python libraries through the [utils](./py/utils) module.
3. New dependencies: `libomp` and `pybind11`.
4. Automated extraction of loop-functions and their degeneracy limits using the [CollectLoopFunctions.wl](./mathematica/CollectLoopFunctions.wl) script.
5. The facility to deal with complex quantities, previously all variables and parameters were assumed to be real.
6. Build system based on `meson` and `ninja`.
7. Matching conditions for simpler UV models ([Singlet Scalar Extension of the SM](./Examples/SSE-to-SMEFTd6/) and [Type I Seesaw model](./Examples/SeesawTypeId6/)) added as examples, in addition to the [MSSM-to-SMEFT](./Examples/MSSM-to-SMEFTd6) matching.

### Updated

1. The directory structure and READMEs have been reorganized, to describe the Mathematica, C++ and Python aspects properly, while also separating the core from the examples.

### Removed

1. Hardcoded SMEFT-Warsaw basis information within [OperatorExport.m](./mathematica/OperatorExport.m)
1. FileIO.h and FileIO.cpp that facilitated reading and writing to files using native C++ code.
2. Build system based on a single `Makefile`.
