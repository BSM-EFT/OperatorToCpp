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


# 2. DEPENDENCY CHECK (PYBIND11)

# Check if python can see pybind11 globally/locally
HAS_PYBIND11=0
if python3 -c "import pybind11" >/dev/null 2>&1 || python -c "import pybind11" >/dev/null 2>&1; then
    HAS_PYBIND11=1
fi

if [ $HAS_PYBIND11 -eq 0 ]; then
    echo "====================================================================="
    echo "   Missing dependency: 'pybind11' is required to create the Python front-end."
    echo "   It can also be installed system-wide using a package manager or within" 
    echo "   your Python environment via:"
    echo ""     
    echo "     pip install pybind11"
    echo "====================================================================="
    exit 1
fi
