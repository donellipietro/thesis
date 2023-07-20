# Define the Rscript command
RSCRIPT := Rscript

# Directories
UTILS_DIR := utils/functions

# Get a list of all R scripts
UTILS_SCRIPTS := $(wildcard $(UTILS_DIR)/*.R)

# List of the calls
UTILS := $(patsubst %.R,%.RData,$(UTILS_SCRIPTS))

# Targets
.PHONY: all build clean

# Default target (all): Run all the main tasks
all: build

# Target: build
# Build the functions needed
build: $(UTILS)

# Target: clean
# Clean intermediate files
clean:
	$(RM) $(UTILS_DIR)/*.RData
	
# Rule for R functions
%.RData: $(patsubst %.RData,%.R, $@)
	$(RSCRIPT) $(patsubst %.RData,%.R, $@)

# Rule for R scripts	
%.RExec: $(patsubst %.RExec,%.R, $@)
	$(RSCRIPT) $(patsubst %.RExec,%.R, $@)