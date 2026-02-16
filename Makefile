# ============================================
# CONFIG
# ============================================

CXX = g++
BASE_CXXFLAGS = -std=c++17 -O3 -march=native -ffast-math -DNDEBUG -Iinclude
OPENMP_CXXFLAGS = -DUSE_OPENMP -fopenmp

MT ?= 0
THREADS ?= 0

ifeq ($(MT),1)
CXXFLAGS = $(BASE_CXXFLAGS) $(OPENMP_CXXFLAGS)
else
CXXFLAGS = $(BASE_CXXFLAGS)
endif

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
ifeq ($(MT),1)
	@if [ "$(THREADS)" = "0" ]; then \
		echo "Running multithreaded (OpenMP default thread count)."; \
		./$(TARGET); \
	else \
		echo "Running multithreaded with OMP_NUM_THREADS=$(THREADS)."; \
		OMP_NUM_THREADS=$(THREADS) ./$(TARGET); \
	fi
else
	@echo "Running single-threaded build."
	./$(TARGET)
endif

# ============================================
# SAVE LATEST OUTPUT SNAPSHOT
# Usage: make save <name>
# ============================================

save:
	@name="$(word 2,$(MAKECMDGOALS))"; \
	if [ -z "$$name" ]; then \
		echo "Usage: make save <name>"; \
		exit 1; \
	fi; \
	if [ ! -d output/latest ]; then \
		echo "No latest output found. Run 'make run' first."; \
		exit 1; \
	fi; \
	mkdir -p output/saved; \
	rm -rf "output/saved/$$name"; \
	cp -r output/latest "output/saved/$$name"; \
	echo "Saved latest output to output/saved/$$name"

%:
	@:

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
