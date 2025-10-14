# Detect platform
UNAME_S = $(shell uname -s)

# Default to c++ (will be overwritten if GCC found)
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wno-sign-compare -Iinclude
LIB =

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

# Directories
SRC_DIR = src
LIB_DIR = lib
OBJ_DIR = obj
OUT_DIR = .

# Object files for each executable
OBJ_FILES_1 = $(OBJ_DIR)/LF.o $(OBJ_DIR)/OperatorImport.o $(OBJ_DIR)/MSSM.o $(OBJ_DIR)/FileIO.o $(OBJ_DIR)/main.o
OBJ_FILES_2 = $(OBJ_DIR)/LF.o $(OBJ_DIR)/OperatorImport.o $(OBJ_DIR)/MSSM.o $(OBJ_DIR)/FileIO.o $(OBJ_DIR)/light_higgsinos.o
OBJ_FILES_3 = $(OBJ_DIR)/LF.o $(OBJ_DIR)/OperatorImport.o $(OBJ_DIR)/MSSM.o $(OBJ_DIR)/FileIO.o $(OBJ_DIR)/benchmark.o

# Executables
TARGET1 = $(OUT_DIR)/main.out
TARGET2 = $(OUT_DIR)/light_higgsinos.out
TARGET3 = $(OUT_DIR)/benchmark.out

# Create obj folder if it does not exist
$(OBJ_DIR):
	@echo "Building for $(UNAME_S)"
	mkdir -p $(OBJ_DIR)

# Rule to link object files for main.cpp
$(TARGET1): $(OBJ_FILES_1) | $(OBJ_DIR)
	$(CXX) $(OBJ_FILES_1) -o $(TARGET1) $(LIBS)

# Rule to link object files for light_higgsinos.cpp
$(TARGET2): $(OBJ_FILES_2) | $(OBJ_DIR)
	$(CXX) $(OBJ_FILES_2) -o $(TARGET2) $(LIBS)

# Rule to link object files for benchmark.cpp
$(TARGET3): $(OBJ_FILES_3) | $(OBJ_DIR)
	$(CXX) $(OBJ_FILES_3) -o $(TARGET3) $(LIBS)

# Compile .cpp files to .o files in src/ directory
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile .cpp files to .o files in lib/ directory
$(OBJ_DIR)/%.o: $(LIB_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean object files and executables
clean:
	rm -f $(OBJ_DIR)/*.o $(TARGET1) $(TARGET2) $(TARGET3)
	rm -rf obj

# Phony targets
.PHONY: all clean
