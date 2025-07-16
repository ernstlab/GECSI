#!/bin/bash

# Ensure Conda is properly initialized
source ~/mambaforge/etc/profile.d/conda.sh


# Check for conda
if ! command -v conda >/dev/null 2>&1; then
  echo "conda is not installed. Please install Miniconda or Anaconda first."
  exit 1
fi

# Set the working directory to the directory of the script
SCRIPT_DIR="$(dirname "$0")"

# Change to that directory
cd "$SCRIPT_DIR" || exit

conda install mamba -n base -c conda-forge


ENV_NAME="gecsi-r"
ENV_FILE="./scripts/env/environment.yml"

# Check if environment exists
if conda env list | grep -qE "^${ENV_NAME}[[:space:]]"; then
    echo "Environment '$ENV_NAME' exists. Updating..."
    mamba env update -n "$ENV_NAME" -f "$ENV_FILE" --prune
else
    echo "Creating new environment '$ENV_NAME'..."
    mamba env create -n "$ENV_NAME" -f "$ENV_FILE"
fi

# echo "Done. You can now run: conda activate $ENV_NAME"
mamba init

source ~/.bashrc

mamba activate gecsi-r

# Install R packages listed in R_environment.txt
Rscript ./scripts/env/install.R