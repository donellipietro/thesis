# Load the devtools package if not already installed
if (!require(devtools)) {
  install.packages("devtools")
}

# Set your GitHub PAT
auth_token <- "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"

# Install the package from GitHub
devtools::install_github("hhroig/pR1FPLS", auth_token = auth_token)