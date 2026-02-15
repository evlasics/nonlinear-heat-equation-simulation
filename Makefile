# ============================================
# CONFIG
# ============================================

CXX = g++
CXXFLAGS = -std=c++17 -O3 -march=native -ffast-math -DNDEBUG -Iinclude

SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin

TARGET = $(BIN_DIR)/simulation

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# ============================================
# DEFAULT
# ============================================

all: $(TARGET)

# ============================================
# LINK
# ============================================

$(TARGET): $(OBJ)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(OBJ) -o $(TARGET)

# ============================================
# COMPILE
# ============================================

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ============================================
# RUN
# ============================================

run: all
	@mkdir -p output
	./$(TARGET)

# ============================================
# CLEAN OUTPUT
# ============================================

clean_output:
	@mkdir -p output
	@rm -rf output/*

# ============================================
# CLEAN EVERYTHING
# ============================================

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR) output
