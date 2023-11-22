from Bio import SeqIO

def split_contigs(contigs, chunk_size=10000):
    split_contigs = []

    for contig in contigs:
        contig_length = len(contig)
        num_chunks = contig_length // chunk_size

        for i in range(num_chunks):
            start = i * chunk_size
            end = start + chunk_size
            chunk = contig[start:end]
            split_contigs.append(chunk)

        # Handle the remaining part of the contig if it's shorter than chunk_size
        if contig_length % chunk_size != 0:
            remaining_chunk = contig[num_chunks * chunk_size:]
            split_contigs.append(remaining_chunk)

    return split_contigs

# Input FASTA file
input_file = input("Enter the input FASTA file path: ")

# Output FASTA file
output_file = input("Enter the output FASTA file path: ")

# Chunk size
chunk_size = int(input("Enter the chunk size (default is 10000): ") or "10000")

# Parse the input FASTA file using Biopython's SeqIO module
contigs = [str(record.seq) for record in SeqIO.parse(input_file, "fasta")]

# Split the contigs into chunks
split_contigs = split_contigs(contigs, chunk_size=chunk_size)

# Write the split contigs to the output FASTA file
with open(output_file, "w") as f:
    for i, chunk in enumerate(split_contigs):
        f.write(">Split{}\n{}\n".format(i+1, chunk))