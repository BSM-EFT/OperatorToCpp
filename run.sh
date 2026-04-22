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
    meson setup "$BUILD_DIR" --prefix=$(pwd)
fi

# 2. COMPILE
echo "--- Compiling ---"
meson compile -C "$BUILD_DIR" -j 4

# 3. INSTALL
meson install -C build