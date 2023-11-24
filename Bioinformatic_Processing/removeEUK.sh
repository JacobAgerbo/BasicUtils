#!/bin/bash

run_analysis() {
    # Load necessary modules
    module load trimmomatic kraken2/2.1.2 pigz/2.6.0

    # Parse named arguments
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --input)
                DATA="$2"
                shift 2
                ;;
            --output)
                OUT="$2"
                shift 2
                ;;
            --sample)
                sample_name="$2"
                shift 2
                ;;
            --kraken-db)
                kraken_db="$2"
                shift 2
                ;;
            *)
                echo "Invalid argument: $1"
                exit 1
                ;;
        esac
    done

    # Create temporary directory
    if [[ -O /home/$USER/tmp && -d /home/$USER/tmp ]]; then
        TMPDIR=/home/$USER/tmp
    else
        rm -rf /home/$USER/tmp 2> /dev/null
        mkdir -p /home/$USER/tmp
        TMPDIR=$(mktemp -d /home/$USER/tmp/XXXX)
    fi

    TMP=$TMPDIR
    TEMP=$TMPDIR

    export TMPDIR TMP TEMP

    # Set input and output file names
    input_file_1=$DATA/${sample_name}_1.fq.gz
    input_file_2=$DATA/${sample_name}_2.fq.gz
    output_file_1=$OUT/${sample_name}_filtered_1.fastq.gz
    output_file_2=$OUT/${sample_name}_filtered_2.fastq.gz

    # Run Trimmomatic to trim low quality reads
    trimmomatic PE -threads $SLURM_CPUS_PER_TASK -phred33 $input_file_1 $input_file_2 $output_file_1 /dev/null $output_file_2 /dev/null LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50

    # Remove eukaryotic data from the metagenome using Kraken2
    DATA="$OUT"
    
    kraken2 --db $kraken_db --threads $SLURM_CPUS_PER_TASK --gzip-compressed --paired $output_file_1 $output_file_2 --output "$sample_name"_12_krk2_500.kraken2 > "$sample_name"_12_krk2_500.kraken2.report

    pigz -p $SLURM_CPUS_PER_TASK $OUT/"$sample_name"_12_krk2_500.kraken2 $OUT/"$sample_name"_12_krk2_500.kraken2.report

    /projects/mjolnir1/apps/bin/extract_kraken_readsBB.py -k $OUT/"$sample_name"_12_krk2_500.kraken2.gz -s1 $DATA/"$sample_name"_filtered_1.fastq.gz -s2 $DATA/"$sample_name"_filtered_2.fastq.gz -o $OUT/noEUK_"$sample_name"_1.fq -o2 $OUT/noEUK_"$sample_name"_2.fq --fastq-output  --exclude --taxid 2759 --include-children  -r $OUT/"$sample_name"_12_krk2_500.kraken2.report.gz

    pigz -p $SLURM_CPUS_PER_TASK $OUT/noEUK_"$sample_name"_1.fq $OUT/noEUK_"$sample_name"_2.fq

    # Check if noEUK files exist before removing filtered files
    if [[ -e $OUT/noEUK_"$sample_name"_1.fq.gz && -e $OUT/noEUK_"$sample_name"_2.fq.gz ]]; then
        # Remove filtered files
        rm $OUT/${sample_name}_filtered_1.fastq.gz
        rm $OUT/${sample_name}_filtered_2.fastq.gz

        # Create done file
        echo ""$sample_name" is done" > $OUT/done."$sample_name"
    else
        echo "Error: noEUK files do not exist."
        exit 1
    fi
}
