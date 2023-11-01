import argparse
from Bio import SeqIO

def subtract_reads(header_file, fasta_file):
    # Read the header names from the file
    with open(header_file, 'r') as file:
        header_names = [line.strip() for line in file]

    # Read the FASTA file and filter based on header names
    sequences = SeqIO.parse(fasta_file, 'fasta')
    filtered_sequences = [seq for seq in sequences if seq.id in header_names]

    # Write the filtered sequences to a new FASTA file
    output_file = 'viral_contigs.fasta'
    SeqIO.write(filtered_sequences, output_file, 'fasta')

    print(f"Filtered reads have been saved to {output_file}.")

if __name__ == '__main__':
    # Create argument parser
    parser = argparse.ArgumentParser(description='Subtract reads from a FASTA file based on header names')
    parser.add_argument('header_file', type=str, help='Path to the file containing header names')
    parser.add_argument('fasta_file', type=str, help='Path to the FASTA file')

    # Parse command-line arguments
    args = parser.parse_args()

    # Call the function with provided arguments
    subtract_reads(args.header_file, args.fasta_file)

