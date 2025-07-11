# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Iinclude

# Directories
SRC_DIR = src
LIB_DIR = lib
OBJ_DIR = obj
OUT_DIR = .

# Object files for each executable
OBJ_FILES_1 = $(OBJ_DIR)/OperatorImport.o $(OBJ_DIR)/MSSM.o $(OBJ_DIR)/FileIO.o $(OBJ_DIR)/main.o
OBJ_FILES_2 = $(OBJ_DIR)/OperatorImport.o $(OBJ_DIR)/MSSM.o $(OBJ_DIR)/FileIO.o $(OBJ_DIR)/light_higgsinos.o

# Executables
TARGET1 = $(OUT_DIR)/main.out
TARGET2 = $(OUT_DIR)/light_higgsinos.out

# Create obj folder if it does not exist
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Rule to link object files for main.cpp
$(TARGET1): $(OBJ_FILES_1) | $(OBJ_DIR)
	$(CXX) $(OBJ_FILES_1) -o $(TARGET1)

# Rule to link object files for light_higgsinos.cpp
$(TARGET2): $(OBJ_FILES_2) | $(OBJ_DIR)
	$(CXX) $(OBJ_FILES_2) -o $(TARGET2)

# Compile .cpp files to .o files in src/ directory
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile .cpp files to .o files in lib/ directory
$(OBJ_DIR)/%.o: $(LIB_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean object files and executables
clean:
	rm -f $(OBJ_DIR)/*.o $(TARGET1)
	rm -rf obj

# Phony targets
.PHONY: all clean
