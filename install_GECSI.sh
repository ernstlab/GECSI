#!/bin/bash
set -e

# Check for conda
if ! command -v conda >/dev/null 2>&1; then
  echo "conda is not installed. Please install Miniconda or Anaconda first."
  exit 1
fi

source $(conda info --base)/etc/profile.d/conda.sh

# Check if mamba is installed
if ! conda list mamba | grep -q mamba; then
  echo "Mamba is not installed. Installing mamba..."
  conda install -y -n base -c conda-forge mamba
fi

# Set the working directory to the directory of the script
SCRIPT_DIR="$(dirname "$0")"

# Change to that directory
cd "$SCRIPT_DIR" || exit

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


# >>> conda initialize >>>
if command -v conda >/dev/null 2>&1; then
    __conda_setup="$(conda 'shell.bash' 'hook' 2> /dev/null)"
    if [ $? -eq 0 ]; then
        eval "$__conda_setup"
    else
        if [ -f "$(dirname \"$(command -v conda)\")/../etc/profile.d/conda.sh" ]; then
            . "$(dirname \"$(command -v conda)\")/../etc/profile.d/conda.sh"
        fi
    fi
    unset __conda_setup
else
    log_error "Conda not found in PATH. Please install or load Conda."
    exit 1
fi
# <<< conda initialize <<<

conda activate gecsi-r

# Install R packages listed in R_environment.txt
Rscript ./scripts/env/install.R

chmod +x ./GECSI.sh

chmod -R +x ./scripts/*.R

chmod -R +x ./preprocessing/*.sh

chmod +x ./run_GECSI_example.sh
