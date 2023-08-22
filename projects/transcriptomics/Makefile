# Define the Rscript command
RSCRIPT := Rscript

# Directories
SCRIPTS_DIR := scripts
FUNCTIONS_DIR := scripts/functions

# Get a list of all R scripts
FUNCTIONS_SCRIPTS := $(wildcard $(FUNCTIONS_DIR)/*.R)

# List of the calls
FUNCTIONS := $(patsubst %.R,%.RData,$(FUNCTIONS_SCRIPTS))

# Targets
.PHONY: all generate_mesh hyperparameters_tuning build clean distclean

# Target: generate_mesh
# Generate the mesh for all the datasets
generate_mesh: build scripts/generate_mesh.RExec

# Target: hyperparameters_tuning
# Perform hyper-parameters tuning for DLPFC data
hyperparameters_tuning: build scripts/hyperparameters_tuning.RExec

# Target: preprocess_data
# Preprocess the data and export csv files
# preprocess_data: build scripts/preprocess_data.RExec

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
	

