# OperatorToC++

[![DOI](https://zenodo.org/badge/948956654.svg)](https://doi.org/10.5281/zenodo.15599937)

OperatorToC++ is an extensible, hybrid Mathematica and C++ based tool that facilitates next steps beyond the matching of parameters between a UV model and an Effective Field Theory. It efficiently remedies the complexities within the analytical matched expressions such as those due to loop-functions and lengthy sums and products involving tensor objects. It then translates and bundles the results into C++ classes and methods which provide a convenient platform for further numerical analyses. 

Additionally, it also provides a convenient Python front-end from where the transpiled classes and functions can be called and incorporated within Python-based HEP workflows. 

It admits [Matchete](https://gitlab.com/matchete/matchete/) (>= v0.3.0) output, stored in a <matching-results>.m file in a key-value pair format, where the keys correspond to EFT coefficients and the values are the respective matching expressions, i.e. functions of the parameters of the specific UV model.

The code can be downloaded by cloning this repository to a local folder.

```
git clone https://github.com/BSM-EFT/OperatorToCpp.git
```
It is possible to interact with the code at 3 distinct levels, i.e. [Mathematica](./mathematica/README.md), [C++](./src/README.md) and [Python](./py/README.md). More details on the C++ class generated after translating the Matchete results and the provided helper functions can be found [here](./lib/README.md).


## Dependencies

To utilise the functionality of OperatorToC++, the following dependencies must be pre-installed. 

1. [Matchete](https://gitlab.com/matchete/matchete/) (>= v0.3.0).

2. A C++ compiler that abides by the C++17 standard or newer, e.g. [GCC](https://gcc.gnu.org/) >= v7.1 and [Clang](https://clang.llvm.org/) >= v4. The use of Clang as the compiler is recommended, on account of the notably faster compile times in comparison to GCC.

3. OpenMP: This is needed to ensure that repetitive operations such Einstein Summations across a large number of repeated indices benefit from parallelization. If one uses Clang as the compiler (default on MacOS) one would need to install `libomp`, while GCC comes bundled with `libgomp`. 

4. [pybind11](https://pybind11.readthedocs.io/en/stable/installing.html) (>= v3.0.0):  It provides the necessary headers for building a Python module using the C++ class and functions. 

5. [meson](https://mesonbuild.com/) + [ninja](https://ninja-build.org/): Our build system comprises of a [meson.build](./meson.build) file which takes care of dependency management as well as cross-platform compilation. 

6. Python libraries: [wcxf](https://wcxf.github.io/), [yaml](https://pypi.org/project/PyYAML/), [numpy](https://numpy.org/), and [pandas](https://pandas.pydata.org/) are required for the proper functioning of the code provided through the [utils](./py/utils/) module.

The development and testing of OperatorToC++ was largely done on a MacOS-based setup. For users working with a similar setup, the C++ dependencies can be installed using the [Homebrew](https://brew.sh/) package manager, using the command,

```
brew install libomp pybind11 meson ninja
```
Similarly, one can use package managers to install the dependencies on Linux distributions. For instance, on Ubuntu, one can use

```
sudo apt install libomp-dev meson ninja-build pybind11-dev
```
However, this can sometimes lead to older package versions, e.g. pybind11 v2.x instead of v3.x. Therefore, for truly cross-platform use, we recommend installing `pybind11`, `meson`and `ninja` within a Python virtual environment, along with the Python dependencies, simply through,

```
pip install pybind11 meson ninja wcxf pyyaml numpy pandas
```

## Citation

If you use this code, please cite:

> Suraj Prakash, *OperatorToC++*. Zenodo, 2025. DOI: [10.5281/zenodo.15599937](https://doi.org/10.5281/zenodo.15599937)

