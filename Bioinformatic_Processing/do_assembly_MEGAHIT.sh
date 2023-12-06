#!/bin/sh
module load miniconda megahit
conda activate vaplinev2
# Default parameter values
DATA=''
OUT=''

# Parse command line arguments
while getopts "d:o:s:" opt; do
  case ${opt} in
    d )
      DATA=${OPTARG}
      ;;
    o )
      OUT=${OPTARG}
      ;;
    s )
      SAMPLE=${OPTARG}
      ;;
    \? )
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    : )
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Check if required parameters are provided
if [ -z "$DATA" ] || [ -z "$OUT" ] || [ -z "$SAMPLE" ]; then
  echo "Usage: bash do_MEGAHIT.sh -d /path/to/input/files -o /path/to/output/files -s sample" 
  exit 1
fi


# Replace <file1> and <file2> with the names of your paired-end read files
file1=noEUK_"$SAMPLE"_1.fq
file2=noEUK_"$SAMPLE"_2.fq

# Count the number of reads in each file
count1=$(zcat $DATA/$file1.gz | echo $((`wc -l`/4)))
count2=$(zcat $DATA/$file2.gz | echo $((`wc -l`/4)))
echo "$file1 contained $count1 reads"
echo "$file2 contained $count2 reads"

# If the files have an unequal number of reads, use BBMap's repair.sh tool to fix them
if [ $count1 -ne $count2 ]; then
    echo "The reads have unequal amount of reads, which will corrected now"
    #repair.sh in1=$DATA/$file1 in2=$DATA/$file2 out1=$DATA/fixed_$file1.fq.gz out2=$DATA/fixed_$file2.fq.gz
        if [[ ! -f $DATA/$file1 && ! -f $DATA/$file2 ]]; then
            THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 1)decompressed files dont exist. Starting decompressing. $THUMBS_UP $(tput sgr 0)"
            pigz -d $DATA/$file1.gz -p $SLURM_CPUS_PER_TASK
            pigz -d $DATA/$file2.gz -p $SLURM_CPUS_PER_TASK
        else
            UNICORN='\U1F984'; echo -e "$(tput setaf 1)Files are ready for fastq_pair ${UNICORN}$(tput sgr 0)"
        fi
    #
        if [[ ! -f $DATA/fixed_"$SAMPLE"_1.fq.gz && ! -f $DATA/fixed_"$SAMPLE"_2.fq.gz ]]; then
            THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 1)paired files dont exist. Starting re-pairing. $THUMBS_UP $(tput sgr 0)"
            fastq_pair $DATA/$file1 $DATA/$file2
            mv $DATA/$file1.paired.fq $DATA/fixed_"$SAMPLE"_1.fq
            mv $DATA/$file2.paired.fq $DATA/fixed_"$SAMPLE"_2.fq
            pigz $DATA/fixed_"$SAMPLE"_1.fq -p $SLURM_CPUS_PER_TASK
            pigz $DATA/fixed_"$SAMPLE"_2.fq -p $SLURM_CPUS_PER_TASK
            pigz $DATA/$file1 -p $SLURM_CPUS_PER_TASK
            pigz $DATA/$file2 -p $SLURM_CPUS_PER_TASK
         else
            UNICORN='\U1F984'; echo -e "$(tput setaf 1)Files are ready for assembly ${UNICORN}$(tput sgr 0)"
        fi
    
    echo "Now I will continue an assembly with fixed reads"
    
    megahit -1 $DATA/fixed_$SAMPLE_1.fq.gz \
            -2 $DATA/fixed_$SAMPLE_2.fq.gz \
        --min-contig-len 1000 \
        -t $SLURM_CPUS_PER_TASK \
        --presets meta-sensitive \
        -o $OUT/$SAMPLE

else
    echo "The paired-end read files already have an equal number of reads."
    echo "Now I will continue an assembly with original reads"
    module load megahit
    megahit -1 $DATA/$file1.gz \
        -2 $DATA/$file2.gz \
        --min-contig-len 1000 \
        -t $SLURM_CPUS_PER_TASK \
        --presets meta-sensitive \
        -o $OUT/$SAMPLE


fi