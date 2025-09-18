#!/bin/bash

# Usage:
# ./bin_bed.sh input_dir bin_size output_dir [genome]

set -euo pipefail
trap 'echo "Error at line $LINENO: command \`$BASH_COMMAND\` failed."' ERR


# Set the working directory to the directory of the script
SCRIPT_DIR="$(dirname "$0")"

# Change to that directory
cd "$SCRIPT_DIR" || exit


# Parse arguments
input_dir="$1"
bin_size="$2"
output_dir="$3"
genome="${4:-hg38}"  # default to hg38
chrom_col="${5:-}"           # chromatin state column number (1-based)



# Validate input
if [[ -z "$input_dir" || -z "$bin_size" || -z "$output_dir" || -z "$chrom_col" ]]; then
  echo "Usage: $0 <input_dir> <bin_size> <output_dir> [genome] <chrom_col>"
  echo "Example: $0 input_dir 200 output_dir hg38 7"
  exit 1
fi

# Check if chrom_col is a positive integer
if ! [[ "$chrom_col" =~ ^[1-9][0-9]*$ ]]; then
  echo "Error: <chrom_col> must be a positive integer."
  exit 1
fi

# Create output directory
mkdir -p "$output_dir"

# Set genome file (chromosome sizes)
if [[ "$genome" == "hg38" ]]; then
  genome_file="../data/hg38.chrom.sizes"
elif [[ "$genome" == "hg19" ]]; then
  genome_file="../data/hg19.chrom.sizes"
else
  echo "Unsupported genome: $genome"
  exit 1
fi

# Filter genome file to chr1–chr22, X (edit if you want Y/M included)
filtered_genome="${output_dir}/filtered.chrom.sizes"
awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X)$/' "$genome_file" > "$filtered_genome"

# Bin generation using bedtools makewindows
bin_file="${output_dir}/bins_${genome}_${bin_size}.bed"
bedtools makewindows -g "$filtered_genome" -w "$bin_size" > "$bin_file"

# Loop through BED files and bin each
for bed in "$input_dir"/*.bed.gz; do
  filename=$(basename "$bed")
  echo "Processing $filename ..."


  # Intersect input BED with bins
  bedtools intersect -a "$bin_file" -b <(zcat "$bed") -wa -wb | \
  awk -v col="$chrom_col" '{print $1, $2, $3, $(col+3)}' OFS='\t' > "${output_dir}/${filename%.bed.gz}_binned.bed"

  bgzip -f "${output_dir}/${filename%.bed.gz}_binned.bed"
done

echo "Binning completed."