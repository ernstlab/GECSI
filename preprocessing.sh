#!/bin/bash

# Set the working directory to the directory of the script
SCRIPT_DIR="$(dirname "$0")"

# Change to that directory
cd "$SCRIPT_DIR" || exit

./preprocessing/bin_bed.sh "../data/training_ct_chromhmm" 200 "../data/training_ct_chromhmm_bins" hg38