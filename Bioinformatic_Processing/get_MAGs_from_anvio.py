#!/usr/bin/env python3

import sys
import csv
import shutil
import argparse

parser = argparse.ArgumentParser(description='Copy contigs.fa from bins from an anvio summary \n with a known SCG domain and a completion above 50 %.')
parser.add_argument('--anvio_summary', metavar='MAGS_Summary_path', type=str, help='Path to MAGS_Summary directory')
parser.add_argument('--bin_list', metavar='input_file_path', type=str, help='Path to the input file')
args = parser.parse_args()

MAGS_SUMMARY_PATH = args.anvio_summary
INPUT_FILE_PATH = args.bin_list

# Extract new_bin_name values to potential_MAGs.txt
with open(INPUT_FILE_PATH, 'r') as file:
    reader = csv.DictReader(file, delimiter='\t')
    potential_MAGs = [row['new_bin_name'] for row in reader if row['SCG_domain'] != 'blank' and float(row['completion']) > 50]


with open('potential_MAGs.txt', 'w') as file:
    file.write('\n'.join(potential_MAGs))

# Loop through each MAG in potential_MAGs.txt and copy corresponding contigs file
with open('potential_MAGs.txt', 'r') as file:
    for i, MAG in enumerate(file):
        MAG = MAG.strip()
        shutil.copy(f"{MAGS_SUMMARY_PATH}/bin_by_bin/{MAG}/{MAG}-contigs.fa", f"./{MAG}-contigs.fa")
        print(f"\033[32mProcessed {i+1} out of {len(potential_MAGs)}: {MAG}-contigs.fa copied.\033[0m")


print("Done!")