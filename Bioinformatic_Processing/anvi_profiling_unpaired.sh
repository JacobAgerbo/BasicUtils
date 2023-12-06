#!/bin/sh
module load anvio bwa samtools
# Set the paths to the input files and output directory
READS="/projects/mjolnir1/people/bfg522/666_Salmon_ARGs/03_MAGs/00_NoEUK"
CONTIGS="/projects/mjolnir1/people/bfg522/666_Salmon_ARGs/03_MAGs/ANVI_TEST/AQ1"
OUT="/projects/mjolnir1/people/bfg522/666_Salmon_ARGs/03_MAGs/ANVI_TEST/AQ1"
# Set the input file names based on the user input
read1="$READS/noEUK_$1_1.fq.gz"
read2="$READS/noEUK_$1_2.fq.gz"

# Print the starting message
echo "Starting PROFILING on: $1"

#So apparently reads are unpaired. Reads will be mapped individually and merged subsequently.
# Align reads to contigs using bwa
bwa mem -t "$SLURM_CPUS_PER_TASK" "$CONTIGS/contigs-fixed.fa" "$read1" > "$OUT/$1_1_aln_pe.sam"
bwa mem -t "$SLURM_CPUS_PER_TASK" "$CONTIGS/contigs-fixed.fa" "$read2" > "$OUT/$1_2_aln_pe.sam"

# Convert SAM to BAM
samtools view -@ "$SLURM_CPUS_PER_TASK" -bS "$OUT/$1_1_aln_pe.sam" > "$OUT/$1_1_aln_pe.bam"
samtools view -@ "$SLURM_CPUS_PER_TASK" -bS "$OUT/$1_2_aln_pe.sam" > "$OUT/$1_2_aln_pe.bam"

# Sort and merge BAMs
# Set the output file name
merged_bam="$1_merged.bam"

# Set the input BAM file names
bam_files=("$1_1_aln_pe.bam" "$1_2_aln_pe.bam")

# Sort and index the input BAM files
for bam in "${bam_files[@]}"
do
    sorted_bam="${bam%.*}.sorted.bam"
    samtools sort -@ "$SLURM_CPUS_PER_TASK" -o $OUT/$sorted_bam $OUT/$bam
    samtools index -@ "$SLURM_CPUS_PER_TASK" $OUT/$sorted_bam
    sorted_bam_files+=($sorted_bam)
done

# Merge the sorted BAM files
samtools merge -@ "$SLURM_CPUS_PER_TASK" "$OUT/$merged_bam" $OUT/"$1"_1_aln_pe.sorted.bam $OUT/"$1"_2_aln_pe.sorted.bam

# Initialize BAM file for anvi'o
anvi-init-bam "$OUT/$merged_bam" -o "$OUT/$1.bam"

# Profile the BAM file using anvi'o
anvi-profile -i "$OUT/$1.bam" -c "$CONTIGS/AQ1_CONTIGS.db" -o "$OUT/$1/" -T "$SLURM_CPUS_PER_TASK" --profile-SCVs

# Print the completion message
echo "Done with: $1"