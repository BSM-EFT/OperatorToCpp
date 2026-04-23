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

# 1. SETUP
if [ ! -d "$BUILD_DIR" ]; then
    echo "--- Initializing Meson Build ---"

    # Check if clang++ exists, otherwise default to g++
    if command -v clang++ >/dev/null 2>&1; then
        export CXX=clang++
        export CC=clang
    fi

    meson setup "$BUILD_DIR" --prefix=$(pwd)
fi

echo "--- Raising system limits ---"
ulimit -s 65536 || echo "Could not increase stack limit"

echo "--- Compiling ---"
# If in Docker/Low RAM, use fewer cores
if [ -f /.dockerenv ]; then
    JOBS=1
else
    JOBS=4
fi

meson compile -C "$BUILD_DIR" -j $JOBS

# 2. COMPILE
# echo "--- Compiling ---"
# meson compile -C "$BUILD_DIR" -j 4

# 3. INSTALL
meson install -C build