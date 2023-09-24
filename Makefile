# Define the Rscript command
RSCRIPT := Rscript

# Directories
UTILS_FUNCTIONS_DIR := utils/functions
UTILS_DIR := utils/

# Get a list of all R scripts
UTILS_FUNCTION_SCRIPTS := $(wildcard $(UTILS_FUNCTION_DIR)/*.R)
UTILS_SCRIPTS := $(wildcard $(UTILS_DIR)/*.R)

# List of the calls
UTILS_FUNCTION := $(patsubst %.R,%.RData,$(UTILS_FUNCTION_SCRIPTS))
UTILS := $(patsubst %.R,%.RData,$(UTILS_SCRIPTS))

# Targets
.PHONY: all build clean

# Default target (all): Run all the main tasks
all: build

# Target: build
# Build the functions needed
build: $(UTILS_FUNCTION)

# Target: clean
# Clean intermediate files
clean:
	$(RM) $(UTILS_FUNCTION_DIR)/*.RData

# Target: install_fdaPDE
# Install fdaPDE core
install_fdaPDE:
	@if [ -d $("/fdaPDE") ]; then \
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