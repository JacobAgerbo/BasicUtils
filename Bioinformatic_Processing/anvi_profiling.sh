profiling() {
    # Load necessary modules
    module load anvio bwa samtools

    # Parse named arguments
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --sample)
                sample_name="$2"
                shift 2
                ;;
            --reads)
                reads_dir="$2"
                shift 2
                ;;
            --contigs)
                contigs_dir="$2"
                shift 2
                ;;
            --output)
                output_dir="$2"
                shift 2
                ;;
            *)
                echo "Invalid argument: $1"
                exit 1
                ;;
        esac
    done

    # Set the paths to the input files and output directory
    READS="$reads_dir"
    CONTIGS="$contigs_dir"
    OUT="$output_dir"

    # Set the input file names based on the user input
    read1="$READS/noEUK_${sample_name}_1.fq.gz"
    read2="$READS/noEUK_${sample_name}_2.fq.gz"

    # Print the starting message
    echo "Starting PROFILING on: $sample_name"

    #So apparently reads are unpaired. Reads will be mapped individually and merged subsequently.
    # Align reads to contigs using bwa
    bwa mem -t "$SLURM_CPUS_PER_TASK" "$CONTIGS/contigs-fixed.fa" "$read1" "$read2" > "$OUT/${sample_name}_aln_pe.sam"

    # Convert SAM to BAM
    samtools view -@ "$SLURM_CPUS_PER_TASK" -bS "$OUT/${sample_name}_aln_pe.sam" > "$OUT/${sample_name}_aln_pe.bam"

    # Initialize BAM file for anvi'o
    anvi-init-bam "$OUT/${sample_name}_aln_pe.bam" -o "$OUT/${sample_name}.bam"

    # Profile the BAM file using anvi'o
    anvi-profile -i "$OUT/${sample_name}.bam" -c "$CONTIGS/CONTIGS.db" -o "$OUT/${sample_name}/" -T "$SLURM_CPUS_PER_TASK" --profile-SCVs

    # Print the completion message
    echo "Done with: $sample_name"
}
