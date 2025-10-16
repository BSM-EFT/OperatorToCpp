# Detect platform
UNAME_S = $(shell uname -s)

# Default to c++ (will be overwritten if GCC found)
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wno-sign-compare -Iinclude
LIBS =

# Enable OpenMP
ifeq ($(UNAME_S),Darwin)
    # Apple Clang + Homebrew libomp
    CXXFLAGS += -Xpreprocessor -fopenmp -I/opt/homebrew/Cellar/libomp/21.1.3/include
    LIBS += -L/opt/homebrew/Cellar/libomp/21.1.3/lib -lomp
else
    # Linux GCC
    CXXFLAGS += -fopenmp
    LIBS += -lgomp
endif

# Precompiled header
PCH = include/pch.h
PCH_GCH = include/pch.pch

# Directories
SRC_DIR = src
LIB_DIR = lib
OBJ_DIR = obj
OUT_DIR = .

LIB_CPP_FILES := $(wildcard $(LIB_DIR)/*.cpp)
LIB_OBJ_FILES := $(patsubst $(LIB_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(LIB_CPP_FILES))

# Object files for each executable (add LIB_OBJ_FILES to all)
OBJ_FILES = $(LIB_OBJ_FILES) $(OBJ_DIR)/main.o

# Executables
TARGET = $(OUT_DIR)/main.out

# Build the precompiled header
$(PCH_GCH): $(PCH)
	$(CXX) $(CXXFLAGS) -x c++-header $< -o $@

# Create obj folder if it does not exist
$(OBJ_DIR):
	@echo "Building for $(UNAME_S)"
	mkdir -p $(OBJ_DIR)

# Rule to link object files for main.cpp
$(TARGET): $(OBJ_FILES) | $(OBJ_DIR)
	$(CXX) $(OBJ_FILES) -o $(TARGET) $(LIBS)

# Compile .cpp files to .o files in src/ directory
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp $(PCH_GCH) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -include $(basename $(PCH)) -c $< -o $@

# Compile .cpp files to .o files in lib/ directory
$(OBJ_DIR)/%.o: $(LIB_DIR)/%.cpp $(PCH_GCH) | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -include $(basename $(PCH)) -c $< -o $@

# Clean object files and executables
clean:
	rm -f $(OBJ_DIR)/*.o $(TARGET) $(PCH_GCH)
	rm -rf obj

# Phony targets
.PHONY: all clean
