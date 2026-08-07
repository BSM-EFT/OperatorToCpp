#!/bin/bash
set -e

# Target directory from command line argument
TARGET_DIR="$1"

# Check if target argument was provided
if [ -z "$TARGET_DIR" ]; then
    echo "Usage: $0 <target_directory>"
    exit 1
fi

# Create target subdirectories if they don't exist
mkdir -p "$TARGET_DIR/include"
mkdir -p "$TARGET_DIR/lib"
mkdir -p "$TARGET_DIR/py/utils"

# Copy C++ lib, and include files into target directory
cp -r "OperatorToCpp/include/." "$TARGET_DIR/include/"
cp -r "OperatorToCpp/lib/." "$TARGET_DIR/lib/"

# Copy all contents of the Python utils library
cp -r "OperatorToCpp/utils/." "$TARGET_DIR/py/utils/"

# 1. BUILD TOOLCHECK (MESON & NINJA)
MISSING_TOOLS=0
if ! command -v meson >/dev/null 2>&1; then MISSING_TOOLS=1; fi
if ! command -v ninja >/dev/null 2>&1; then MISSING_TOOLS=1; fi

if [ $MISSING_TOOLS -eq 1 ]; then
    echo "====================================================================="
    echo "  Missing dependencies: 'meson' and 'ninja' are needed to build the project."
    echo "  These can be installed by:"
    echo "   (i)  Using a package manager on Linux-based systems ('brew' for macOS)"
    echo "          https://mesonbuild.com/Getting-meson.html"
    echo "   (ii) Within your Python environment using 'pip':"
    echo "          pip install meson ninja"
    echo "====================================================================="
    exit 1
fi

# Copy build system files to the target directory
cp "OperatorToCpp/meson.build" "$TARGET_DIR/meson.build"
cp "OperatorToCpp/compileCpp.sh" "$TARGET_DIR/compileCpp.sh"

# 2. COMPILER & OPENMP CHECK

# Determine compiler choice early
if command -v clang++ >/dev/null 2>&1; then
    DEFAULT_CXX="clang++"
    export CXX=clang++
    export CC=clang
elif command -v g++ >/dev/null 2>&1; then
    DEFAULT_CXX="g++"
else
    echo "Error: Neither clang++ nor g++ was found on your system."
    exit 1
fi

if [ "$DEFAULT_CXX" == "g++" ]; then
    echo "GCC ('g++') detected as the default compiler. (Note: We recommend using Clang for faster compilation times)."
else
    echo "Clang ('clang++') detected as the default compiler."
    # Verify libomp presence for Clang
    LIBOMP_FOUND=0
    if [ "$(uname)" == "Darwin" ]; then
        if [ -d "/opt/homebrew/opt/libomp" ] || [ -d "/usr/local/opt/libomp" ]; then
            LIBOMP_FOUND=1
        fi
    else
        # Linux check via fallback compilation test or dpkg/pacman if needed, 
        # but standard Clang environments typically drop it in default search paths.
        if echo "#include <omp.h>" | clang++ -E -x c++ - >/dev/null 2>&1; then
            LIBOMP_FOUND=1
        fi
    fi

    if [ $LIBOMP_FOUND -eq 0 ]; then
        echo "====================================================================="
        echo "   Missing Dependency: 'libomp' is needed for parallel execution of code."
        echo "   It can be installed using Linux package managers (or 'brew' on macOS):"
        echo ""     
        echo "     brew install libomp"
        echo "====================================================================="
        exit 1
    fi
fi

# 3. PYBIND11 WARNING (NON-FATAL)

# Check if python can see pybind11 globally/locally
HAS_PYBIND11=0
if python3 -c "import pybind11" >/dev/null 2>&1 || python -c "import pybind11" >/dev/null 2>&1; then
    HAS_PYBIND11=1
fi

if [ $HAS_PYBIND11 -eq 0 ]; then
    echo "====================================================================="
    echo "   Missing dependency: 'pybind11' is required to create the Python front-end."
    echo "   The build system can automatically download it via WrapDB during setup,"
    echo "   otherwise it can also be installed within your Python environment via:"
    echo ""     
    echo "     pip install pybind11"
    echo "====================================================================="
    echo "Proceeding with automatic build system handling..."
    mkdir -p "$TARGET_DIR/subprojects"
    cp "OperatorToCpp/subprojects/." "$TARGET_DIR/subprojects/"
fi
