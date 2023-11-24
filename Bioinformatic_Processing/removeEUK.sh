#!/bin/sh
conda activate Basic_Utils # load biopython
module load trimmomatic kraken2/2.1.2 pigz/2.6.0
DATA='/projects/mjolnir1/people/bfg522/01_FlyProject/00_RAWDATA'
OUT='/projects/mjolnir1/people/bfg522/01_FlyProject/01_NOEUK'

if [[ -O /home/$USER/tmp && -d /home/$USER/tmp ]]; then
        TMPDIR=/home/$USER/tmp
else
        # You may wish to remove this line, it is there in case
        # a user has put a file 'tmp' in there directory or a
        rm -rf /home/$USER/tmp 2> /dev/null
        mkdir -p /home/$USER/tmp
        TMPDIR=$(mktemp -d /home/$USER/tmp/XXXX)
fi

TMP=$TMPDIR
TEMP=$TMPDIR

export TMPDIR TMP TEMP

# Trim data to remove low quality reads
# Define the quality filtering function
  # Set the sample name and input file names
  sample_name=$1
  input_file_1=$DATA/${sample_name}_1.fq.gz
  input_file_2=$DATA/${sample_name}_2.fq.gz

  # Set the output file names
  output_file_1=$OUT/${sample_name}_filtered_1.fastq.gz
  output_file_2=$OUT/${sample_name}_filtered_2.fastq.gz

  # Run Trimmomatic with the desired settings
  trimmomatic PE -threads $SLURM_CPUS_PER_TASK -phred33 $input_file_1 $input_file_2 $output_file_1 /dev/null $output_file_2 /dev/null LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50

# Remove eukaryotic data from the metagenome
DATA=$OUT
DB="/projects/mjolnir1/data/databases/kraken2/kraken2_standard/20220926/"
/projects/mjolnir1/apps/conda/kraken2-2.1.2/bin/kraken2 --db $DB \
	--paired --use-names \
	--threads $SLURM_CPUS_PER_TASK \
	--output $OUT/"$1"_12_krk2_500.kraken2 --report $OUT/"$1"_12_krk2_500.kraken2.report \
	$DATA/"$1"_filtered_1.fastq.gz  $DATA/"$1"_filtered_2.fastq.gz

pigz -p $SLURM_CPUS_PER_TASK $OUT/"$1"_12_krk2_500.kraken2 $OUT/"$1"_12_krk2_500.kraken2.report

if [[ -O /home/$USER/tmp && -d /home/$USER/tmp ]]; then
        TMPDIR=/home/$USER/tmp
else
        # You may wish to remove this line, it is there in case
        # a user has put a file 'tmp' in there directory or a
        rm -rf /home/$USER/tmp 2> /dev/null
        mkdir -p /home/$USER/tmp
        TMPDIR=$(mktemp -d /home/$USER/tmp/XXXX)
fi

TMP=$TMPDIR
TEMP=$TMPDIR

export TMPDIR TMP TEMP

# Remove eukaryotic data from the metagenome

/projects/mjolnir1/apps/bin/extract_kraken_readsBB.py -k $OUT/"$1"_12_krk2_500.kraken2.gz -s1 $DATA/"$1"_filtered_1.fastq.gz -s2 $DATA/"$1"_filtered_2.fastq.gz -o $OUT/noEUK_"$1"_1.fq -o2 $OUT/noEUK_"$1"_2.fq --fastq-output  --exclude --taxid 2759 --include-children  -r $OUT/"$1"_12_krk2_500.kraken2.report.gz

pigz -p $SLURM_CPUS_PER_TASK $OUT/noEUK_"$1"_1.fq $OUT/noEUK_"$1"_2.fq
echo ""$1" is done" > $OUT/done."$1"
#rm $DATA/"$1"_filtered_*.fastq.gz