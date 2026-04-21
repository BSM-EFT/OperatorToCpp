#!/bin/bash
set -e

# -----------------------------
# Configuration
# -----------------------------
BUILD_DIR="build"
PY_INSTALL_DIR="py"
PYTHON_EXEC=$(which python)

EXECUTABLE_NAME="main.out"

# -----------------------------
# Clean
# -----------------------------
if [ "$1" = "clean" ]; then
    echo "Cleaning build artifacts..."
    rm -rf "$BUILD_DIR"
    rm -f "$PY_INSTALL_DIR"/*.so
    rm -f "$EXECUTABLE_NAME"
    exit 0
fi

# -----------------------------
# Configure
# -----------------------------
PYBIND_DIR=$($PYTHON_EXEC -m pybind11 --cmakedir)
if [ ! -d "$BUILD_DIR" ]; then
    echo "Configuring project..."
    cmake -B "$BUILD_DIR" -S . -DPython_EXECUTABLE="$PYTHON_EXEC" -DCMAKE_PREFIX_PATH="$PYBIND_DIR"
fi

# -----------------------------
# Build incrementally
# -----------------------------
echo "Building project..."
cmake --build "$BUILD_DIR" --config Debug -- -j$(nproc)

# -----------------------------
# Install
# -----------------------------
echo "Installing targets..."

# Ensure python install directory exists
mkdir -p "$PY_INSTALL_DIR"

# Install main executable to project root
cp "$BUILD_DIR/$EXECUTABLE_NAME" "./$EXECUTABLE_NAME"

# Install Python module to py/
# Note: CMake install already knows the install path from CMakeLists.txt
cmake --install "$BUILD_DIR" --prefix "$PY_INSTALL_DIR"

echo "Build and install complete."
