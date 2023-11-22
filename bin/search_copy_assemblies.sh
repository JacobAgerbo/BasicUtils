#!/bin/bash
# Find immediate subfolders in the current directory
subfolders=$(find . -maxdepth 1 -type d)
# Loop through each subfolder
for subfolder in $subfolders; do
    # Check if contigs.fasta exists in the subfolder
    if [ -f "$subfolder/contigs.fasta" ]; then

        # Get the folder name
        folder_name="${subfolder##*/}"
	echo "$folder_name" >> temp
        # Create a new filename with folder name and contigs.fasta
        new_filename="${folder_name}_contigs.fasta"

        # Check if the renamed file already exists in the parent folder
        if [ -f "./$new_filename" ]; then
            echo "The file $new_filename already exists in this folder. Skipping..."
        else
            # Copy and rename contigs.fasta to the new filename
            cp "$subfolder/contigs.fasta" ./"$new_filename"
            echo "contigs.fasta copied and renamed to $new_filename from $subfolder to current folder."
        fi

    else
        echo "contigs.fasta does not exist in $subfolder."

    fi

done
sort temp > Assembly_Overview.txt
rm temp