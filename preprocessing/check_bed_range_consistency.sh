#!/bin/bash
set -euo pipefail

input_dir="$1"

if [[ ! -d "$input_dir" ]]; then
  echo "Directory not found: $input_dir"
  exit 1
fi

temp_dir=$(mktemp -d)
echo "Checking chromosomal ranges in: $input_dir"

# Loop through each file and extract min/max per chromosome
for file in "$input_dir"/*.bed.gz; do
  ct=$(basename "$file" .bed.gz)
  zcat "$file" | awk '{
    chr = $1;
    start = $2;
    end = $3;
    if (chr in min && start < min[chr]) min[chr] = start;
    else if (!(chr in min)) min[chr] = start;
    if (chr in max && end > max[chr]) max[chr] = end;
    else if (!(chr in max)) max[chr] = end;
  } END {
    for (c in min)
      print c, min[c], max[c];
  }' | sort > "$temp_dir/${ct}.range"
done

# Compare ranges
ref_file=$(ls "$temp_dir" | head -n 1)
ref_path="$temp_dir/$ref_file"
consistent=true

for file in "$temp_dir"/*.range; do
  if ! diff -q "$ref_path" "$file" >/dev/null; then
    echo "Inconsistent ranges found in: $(basename "$file" .range)"
    consistent=false
  fi
done

if $consistent; then
  echo "All files have consistent chromosome ranges."
else
  echo "Some files have inconsistent ranges. Please review."
fi

rm -r "$temp_dir"