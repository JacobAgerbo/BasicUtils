#!/bin/sh
# Iterate over each sample in the list
mapfile -t sample_list < $1
for sample in "${sample_list[@]}"; do
  # Check if both files exist for the current sample
  if [[ -e noEUK_"${sample}_1.fq.gz" && -e noEUK_"${sample}_2.fq.gz" ]]; then
    echo "$sample" >> start_assembly_samples.txt
  fi
done