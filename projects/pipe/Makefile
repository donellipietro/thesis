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

# Target: generate_mesh
# Generate the two meshes (alpha = 10, alpha = 20)
generate_mesh: build scripts/generate_mesh.RExec

# Target: preprocess_data
# Preprocess the data and export csv files
preprocess_data: build scripts/preprocess_data.RExec

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