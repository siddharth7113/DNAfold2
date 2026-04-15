# DNAfold2 Makefile
# Compiles C source files with optimizations and OpenMP support

CC = gcc
CXX = g++
CFLAGS = -O3 -Wall
CFLAGS_OMP = -O3 -Wall -fopenmp
LDFLAGS = -lm

# Directories
SRC_DIR = src
BIN_DIR = bin
CORE_DIR = $(SRC_DIR)/core
ANALYSIS_DIR = $(SRC_DIR)/analysis
UTILS_DIR = $(SRC_DIR)/utils
REBUILD_DIR = $(SRC_DIR)/rebuild
INITIAL_DIR = $(SRC_DIR)/initial

# Core binaries (buildable from source)
# Note: TiRNA_remc is a pre-compiled binary (no source available)
#       It must be preserved in bin/ - do not delete it
CORE_BINS = $(BIN_DIR)/TiRNA_sa $(BIN_DIR)/TiRNA_optimize

# Analysis binaries
ANALYSIS_BINS = $(BIN_DIR)/secondary $(BIN_DIR)/wham $(BIN_DIR)/A_state $(BIN_DIR)/A_state1

# Utility binaries
UTIL_BINS = $(BIN_DIR)/center $(BIN_DIR)/tc $(BIN_DIR)/cat $(BIN_DIR)/t1

.PHONY: all clean install core analysis utils help

all: $(BIN_DIR) core analysis utils
	@echo "Build complete. Binaries are in $(BIN_DIR)/"

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Core algorithms
core: $(CORE_BINS)

$(BIN_DIR)/TiRNA_sa: $(CORE_DIR)/TiRNA_sa.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

$(BIN_DIR)/TiRNA_optimize: $(CORE_DIR)/TiRNA_optimize.c | $(BIN_DIR)
	$(CC) $(CFLAGS_OMP) $< -o $@ $(LDFLAGS)

# Analysis tools
analysis: $(ANALYSIS_BINS)

$(BIN_DIR)/secondary: $(ANALYSIS_DIR)/secondary.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

$(BIN_DIR)/wham: $(ANALYSIS_DIR)/wham.c | $(BIN_DIR)
	$(CXX) $(CFLAGS) $< -o $@

$(BIN_DIR)/A_state: $(ANALYSIS_DIR)/A_state.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

$(BIN_DIR)/A_state1: $(ANALYSIS_DIR)/A_state1.c | $(BIN_DIR)
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

# Utility programs
utils: $(UTIL_BINS)

$(BIN_DIR)/center: $(UTILS_DIR)/center.c | $(BIN_DIR)
	$(CXX) $(CFLAGS) $< -o $@

$(BIN_DIR)/tc: $(UTILS_DIR)/tc.c | $(BIN_DIR)
	$(CXX) $(CFLAGS) $< -o $@

$(BIN_DIR)/cat: $(UTILS_DIR)/cat.c | $(BIN_DIR)
	$(CXX) $(CFLAGS) $< -o $@

$(BIN_DIR)/t1: $(UTILS_DIR)/t1.c | $(BIN_DIR)
	$(CXX) $(CFLAGS) $< -o $@

# Rebuild tools (in subdirectory)
rebuild:
	cd $(REBUILD_DIR) && $(CC) $(CFLAGS) rebuild.c -o rebuild $(LDFLAGS)
	cd $(REBUILD_DIR) && $(CC) $(CFLAGS) secondary.c -o secondary $(LDFLAGS)

# Initial sequence processing
initial:
	cd $(INITIAL_DIR) && $(CXX) $(CFLAGS) seq_initial.c -o seq_initial

# Install Python package
install:
	pip install -e .

# Install with development dependencies
install-dev:
	pip install -e ".[dev]"

# Run tests
test:
	pytest tests/ -v

# Clean build artifacts
clean:
	rm -rf $(BIN_DIR)/*.o
	rm -rf $(BIN_DIR)/a.out
	find . -name "*.o" -delete
	find . -name "a.out" -delete
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
	find . -name "*.pyc" -delete

# Full clean including binaries
distclean: clean
	rm -rf $(BIN_DIR)
	rm -rf build/ dist/ *.egg-info/

help:
	@echo "DNAfold2 Build System"
	@echo ""
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  all         Build all binaries (default)"
	@echo "  core        Build core algorithms (TiRNA_sa, TiRNA_optimize)"
	@echo "  analysis    Build analysis tools (secondary, wham, A_state)"
	@echo "  utils       Build utility programs (center, tc, cat, t1)"
	@echo "  rebuild     Build rebuild tools"
	@echo "  initial     Build sequence initialization"
	@echo "  install     Install Python package"
	@echo "  install-dev Install with development dependencies"
	@echo "  test        Run tests"
	@echo "  clean       Remove object files"
	@echo "  distclean   Remove all build artifacts"
	@echo "  help        Show this help message"
