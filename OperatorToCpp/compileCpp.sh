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

# 1. INITIALIZE & COMPILE
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