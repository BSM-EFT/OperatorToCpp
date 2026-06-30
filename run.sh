#!/bin/bash
set -e

BUILD_DIR="build"

# 0. CLEANUP
if [ "$1" == "clean" ]; then
    echo "Cleaning project..."
    rm -rf build
    rm -f ./main.out
    rm -f ./py/*.so
    echo "Done."
    exit 0
fi

# ==========================================
# 0. SANITY CHECK: REQUIRED DIRECTORIES
# ==========================================
REQUIRED_DIRS=("mathematica" "include" "lib" "src" "py")
for dir in "${REQUIRED_DIRS[@]}"; do
    if [ ! -d "./$dir" ]; then
        echo "Error: Missing project files - no '$dir' directory found."
        exit 1
    fi
done

# ==========================================
# 1. BUILD TOOLCHECK (MESON & NINJA)
# ==========================================
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

# ==========================================
# 2. COMPILER & OPENMP CHECK
# ==========================================
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

# ==========================================
# 3. PYBIND11 WARNING (NON-FATAL)
# ==========================================
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
fi

# ==========================================
# 4. INITIALIZE & COMPILE
# ==========================================
if [ ! -d "$BUILD_DIR" ]; then
    echo "--- Initializing Meson Build ---"
    meson setup "$BUILD_DIR" --prefix=$(pwd)
fi

echo "--- Raising system limits ---"
ulimit -s 65536 || echo "Could not increase stack limit"

echo "--- Compiling ---"
if [ -f /.dockerenv ]; then
    JOBS=1
else
    JOBS=4
fi

meson compile -C "$BUILD_DIR" -j $JOBS

# 5. INSTALL
meson install -C build