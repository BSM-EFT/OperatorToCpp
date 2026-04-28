## Since v0.1

### Added

1. A new Python front-end and the interop with Python libraries through the \texttt{utils} module.
2. New dependencies: `libomp` and `pybind11`.
3. Automated extraction of loop-functions and their degeneracy limits using the `CollectLoopFunctions.wl` script.
4. The facility to deal with complex quantities, previously all variables and parameters were assumed to be real.
5. Build system based on `meson` and `ninja`.

### Updated

1. The directory structure and READMEs hasve been reorganized across folders, to describe the Mathematic, C++ and Python aspects properly.


### Removed

1. FileIO.h and FileIO.cpp that facilitated reading and writing to files using native C++ code.
2. Build system based on a single `Makefile`.
