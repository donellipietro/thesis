# Define the Rscript command
RSCRIPT := Rscript

# Directories
SCRIPTS_DIR := scripts
FUNCTIONS_DIR := scripts/functions

# Get a list of all R scripts
TEST_SCRIPTS := $(wildcard $(SCRIPTS_DIR)/test*.R)
ANALYSIS_SCRIPTS := $(wildcard $(SCRIPTS_DIR)/analysis*.R)
FUNCTIONS_SCRIPTS := $(wildcard $(FUNCTIONS_DIR)/*.R)

# List of the calls
TESTS := $(patsubst %.R,%.RExec,$(TEST_SCRIPTS))
ANALYSIS := $(patsubst %.R,%.RExec,$(ANALYSIS_SCRIPTS))
FUNCTIONS := $(patsubst %.R,%.RData,$(FUNCTIONS_SCRIPTS))

# Targets
.PHONY: all run_tests run_analysis build clean distclean

# Default target (all): Run all the main tasks
all: run_tests run_analysis

# Target: run_tests
# Execute test R scripts
run_tests: build $(TESTS)

# Target: run_analysis
# Execute analysis R scripts
run_analysis: build $(ANALYSIS)

# Target: build
# Build the functions and directories needed
build: $(FUNCTIONS)
	make -C ./../../. build
	mkdir -p images results data
	
# Target: clean
# Clean intermediate files
clean:
	$(RM) Rplots.pdf
	$(RM) $(FUNCTIONS_DIR)/*.RData
	
# Target: distclean
# Clean all generated files
distclean: clean
	$(RM) data/ -r
	$(RM) images/ -r
	$(RM) results/ -r
	
# Rule for R functions
%.RData: $(patsubst %.RData,%.R, $@)
	$(RSCRIPT) $(patsubst %.RData,%.R, $@)

# Rule for R scripts	
%.RExec: $(patsubst %.RExec,%.R, $@)
	$(RSCRIPT) $(patsubst %.RExec,%.R, $@)
	

