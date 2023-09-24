# Define the Rscript command
RSCRIPT := Rscript

# Directories
UTILS_FUNCTIONS_DIR := utils/functions
UTILS_DIR := utils/

# Get a list of all R scripts
UTILS_FUNCTIONS_SCRIPTS := $(wildcard $(UTILS_FUNCTIONS_DIR)/*.R)

# List of the calls
UTILS_FUNCTIONS := $(patsubst %.R,%.RData,$(UTILS_FUNCTIONS_SCRIPTS))

# Targets
.PHONY: all build clean

# Default target (all): Run all the main tasks
all: build

# Target: build
# Build the functions needed
build: $(UTILS_FUNCTIONS)

# Target: clean
# Clean intermediate files
clean:
	$(RM) $(UTILS_FUNCTIONS_DIR)/*.RData

# Target: install_fdaPDE
# Install fdaPDE core
DIR_TO_CHECK := fdaPDE
install_fdaPDE:
	@if [ -d $(DIR_TO_CHECK) ]; then \
    	echo "fdaPDE repository already exists."; \
    else \
		echo "Cloning fdaPDE repository."; \
        git clone https://github.com/donellipietro/fdaPDE.git; \
    fi
	cd fdaPDE && git checkout develop_donelli

# Target: install_libraries
# Install fdaPDE2 and pR1FPLS
install_fdaPDE2: install_fdaPDE utils/fdaPDE2_install.RExec
install_pR1FPLS: utils/pR1FPLS_install.RExec
install_libraries: install_fdaPDE2 install_pR1FPLS
	
# Rule for R functions
%.RData: $(patsubst %.RData,%.R, $@)
	$(RSCRIPT) $(patsubst %.RData,%.R, $@)

# Rule for R scripts	
%.RExec: $(patsubst %.RExec,%.R, $@)
	$(RSCRIPT) $(patsubst %.RExec,%.R, $@)