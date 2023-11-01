#!/bin/bash
module load trimmomatic flye medaka busco spades
# Define paths to input and output files
# Default parameter values
LONG_DATA=""
SHORT_DATA=""
LR_input_dir=$LONG_DATA/$1
SR_input_dir=$SHORT_DATA/$1
output_dir=""

# Parse command line arguments
while getopts "l:s:o:" opt; do
  case ${opt} in
    l )
        LONG_DATA=${OPTARG}
      ;;
    s )
        SHORT_DATA=${OPTARG}
      ;;
    o )
      output_dir="${OPTARG}/$1"
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
if [ -z "$LONG_DATA" ] || [ -z "$SHORT_DATA" ] || [ -z "$output_dir" ]; then
        echo -e "$(tput setaf 2)Usage: bash do_hybrid_assembly.sh -l long_reads_dir -s short_reads_dir -o output_dir $(tput sgr 0)"
        echo -e "$(tput setaf 1)Please provide the required parameters: -l for long reads directory, -s for short reads directory, and -o for output directory.$(tput sgr 0)"
exit 1
fi

# Define command-line tools and parameters
trimmomatic="trimmomatic"
flye="flye"
busco="busco"
spades="metaspades.py"
quickmerge="merge_wrapper.py"

BIN="$USER/.bin"

if [[ ! -d $BIN ]]; then
    THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 1)misc files for this script doent exist.. $(tput sgr 0)"
    THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 2)Files will be available at: $BIN $THUMBS_UP $(tput sgr 0)"
    wget -O $BIN 

fi

# Create output directory if it does not exist
mkdir -p $output_dir

# Perform error correction on Illumina reads using Trimmomatic
echo "$(tput setaf 2)###################################################$(tput sgr 0)"
echo "$(tput setaf 2)## Performing error correction on Illumina reads.##$(tput sgr 0)"
echo "$(tput setaf 2)###################################################$(tput sgr 0)"

if [[ ! -f $output_dir/"$1"_paired_1.fq.gz && ! -f $output_dir/"$1"_paired_2.fq.gz ]]; then
    THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 1)Paired files dont exist. Starting error correction. $THUMBS_UP $(tput sgr 0)"
    trimmomatic PE -threads $SLURM_CPUS_PER_TASK \
        $SR_input_dir/"$1"_1.fq.gz \
        $SR_input_dir/"$1"_2.fq.gz \
        $output_dir/"$1"_paired_1.fq.gz $output_dir/"$1"_unpaired_1.fq.gz \
        $output_dir/"$1"_paired_2.fq.gz $output_dir/"$1"_unpaired_2.fq.gz \
        -phred33 \
        LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50
else
    UNICORN='\U1F984'; echo -e "$(tput setaf 1)Paired files already exist. Skipping error correction. ${UNICORN}$(tput sgr 0)"

fi

# Assemble long reads using Flye
echo "$(tput setaf 2)####################################################$(tput sgr 0)"
echo "$(tput setaf 2)##         Assembling long reads with Flye.       ##$(tput sgr 0)"
echo "$(tput setaf 2)####################################################$(tput sgr 0)"

if [ ! -f "$output_dir/flye_assembly/assembly.fasta" ]; then
    THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 1)Performing assembly using Flye. $THUMBS_UP $(tput sgr 0)"
    $flye --meta --nano-raw $LR_input_dir/*_long.fastq.gz --out-dir $output_dir/flye_assembly \
        -t $SLURM_CPUS_PER_TASK

    echo "Assembly completed."
else
    UNICORN='\U1F984'; echo -e "$(tput setaf 1)Flye assembly directory already exists. Skipping assembly. ${UNICORN}$(tput sgr 0)"
fi

# Mapping data to assembly for Racon
## Starting with indexing, if needed.
echo "$(tput setaf 2)####################################################$(tput sgr 0)"
echo "$(tput setaf 2)##             Polishing with Racon               ##$(tput sgr 0)"
echo "$(tput setaf 2)####################################################$(tput sgr 0)"
if [ ! -f "$output_dir/flye_assembly/ref.mmi" ]; then
    module load minimap2 samtools
    THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 1)Indexing assembly for minimap2. $THUMBS_UP $(tput sgr 0)"
    minimap2 -d "$output_dir/flye_assembly/ref.mmi" "$output_dir/flye_assembly/assembly.fasta"
    echo "Assembly is now indexed"
else
    UNICORN='\U1F984'; echo -e "$(tput setaf 1)Assembly is indexed. I will continue mapping ${UNICORN}$(tput sgr 0)"
fi
# Mapping of reads to assembly
if [ ! -f "$output_dir/flye_assembly/read_2_assembly.sam" ]; then
    THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 1)Align reads to assembly with minimap2. $THUMBS_UP $(tput sgr 0)"
    module load minimap2 samtools
    minimap2 -a -t $SLURM_CPUS_PER_TASK "$output_dir/flye_assembly/assembly.fasta" $LR_input_dir/*_long.fastq.gz  > "$output_dir/flye_assembly/read_2_assembly.sam"
    echo "Read are now alligned to assembly"
else
    UNICORN='\U1F984'; echo -e "$(tput setaf 1)Reads are alligned to assembly. I will continue with Racon. ${UNICORN}$(tput sgr 0)"
fi

# Polish the assembly using Racon
if [ ! -f "$output_dir/racon.fasta" ]; then
    echo "Polish with Racon"
    module load minimap2 samtools racon
    racon $LR_input_dir/*_long.fastq.gz "$output_dir/flye_assembly/read_2_assembly.sam" "$output_dir/flye_assembly/assembly.fasta" \
        -t $SLURM_CPUS_PER_TASK \
        -m 8 -x -6 -g -8 -w 500 > "$output_dir/racon.fasta"
    THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 1)Assembly is now polished with Racon. $(tput sgr 0)"
    THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 1)Will continue with hybrid assembly with SPAdes. $THUMBS_UP $(tput sgr 0)"
else
    UNICORN='\U1F984'; echo -e "$(tput setaf 1)Assembly is already polished with Racon.$(tput sgr 0)"
    UNICORN='\U1F984'; echo -e "$(tput setaf 1)Will continue with hybrid assembly with SPAdes. ${UNICORN}$(tput sgr 0)"
fi

# Assemble short reads using SPAdes, together with assembled long reads
echo "$(tput setaf 2)####################################################$(tput sgr 0)"
echo "$(tput setaf 2)##         Assemble short reads using SPAdes,     ##$(tput sgr 0)"
echo "$(tput setaf 2)##         together with assembled long reads     ##$(tput sgr 0)"
echo "$(tput setaf 2)####################################################$(tput sgr 0)"
if [ ! -f "$output_dir/spades_assembly/scaffolds.fasta" ]; then
    echo "Assembling short reads with SPAdes.."
$spades -t $SLURM_CPUS_PER_TASK \
        -k 21,33,99,127 \
        --pe1-1 $output_dir/"$1"_paired_1.fq.gz \
        --pe1-2 $output_dir/"$1"_paired_2.fq.gz \
        --nanopore $output_dir/racon.fasta \
        -o $output_dir/spades_assembly
    THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 1)Short reads are now assembled with SPAdes. $THUMBS_UP $(tput sgr 0)"
else
    UNICORN='\U1F984'; echo -e "$(tput setaf 1)Short reads are already assembled with SPAdes. ${UNICORN}$(tput sgr 0)"
    UNICORN='\U1F984'; echo -e "$(tput setaf 1)Will continue to merging of assemblies. ${UNICORN}$(tput sgr 0)"
fi

# Merge the assemblies using Quickmerge
echo "Merging assemblies with Quickmerge..."
conda activate quickmerge
$quickmerge $output_dir/spades_assembly/scaffolds.fasta $output_dir/racon.fasta -pre $1

# Evaluate the quality of the assembly using BUSCO
#echo "Evaluating assembly quality with BUSCO..."
$busco -i $output_dir/racon.fasta -o busco_results -l bacteria_odb10 -m genome --cpu $SLURM_CPUS_PER_TASK

### Make log report
bash $BIN/make_log.sh $1 $output_dir > $output_dir/hybrid_assembly.log