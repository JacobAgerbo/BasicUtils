#!/bin/sh
#####
# Define paths to input and output files
output_dir=$2
BIN="$USER/.bin"
echo "$(tput setaf 2)All generated data can found here: $(tput sgr 0)"
echo "$output_dir"
echo ""
echo "$(tput setaf 2)###################################################$(tput sgr 0)"
echo "$(tput setaf 2)## Performing error correction on Illumina reads.##$(tput sgr 0)"
echo "$(tput setaf 2)###################################################$(tput sgr 0)"
echo ""
if [[ ! -f $output_dir/"$1"_paired_1.fq.gz && ! -f $output_dir/"$1"_paired_2.fq.gz ]]; then
    SAD='\U1F641'; echo -e "$(tput setaf 1)Short reads has not been trimmed. $SAD $(tput sgr 0)"
else
    THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 2)Short reads has been corrected $THUMBS_UP $(tput sgr 0)"
fi
echo ""
# Assemble long reads using Flye
echo "$(tput setaf 2)####################################################$(tput sgr 0)"
echo "$(tput setaf 2)##         Assembling long reads with Flye.       ##$(tput sgr 0)"
echo "$(tput setaf 2)####################################################$(tput sgr 0)"
echo ""
if [ ! -f "$output_dir/flye_assembly/assembly.fasta" ]; then
    SAD='\U1F641'; echo -e "$(tput setaf 1)Long read-based assembly has not been assembled. $SAD $(tput sgr 0)"
else
    THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 2)Long reads has been assembled $THUMBS_UP $(tput sgr 0)"
    echo ""
    bash $BIN/get_long_contig_info.sh "$output_dir/flye_assembly/40-polishing/contigs_stats.txt"
fi
echo ""
# Mapping data to assembly for Racon
## Starting with indexing, if needed.
echo "$(tput setaf 2)####################################################$(tput sgr 0)"
echo "$(tput setaf 2)##         Polishing with Racon and Medaka        ##$(tput sgr 0)"
echo "$(tput setaf 2)####################################################$(tput sgr 0)"
echo ""
# Polish the assembly using Racon
if [ ! -f "$output_dir/racon.fasta" ]; then
    SAD='\U1F641'; echo -e "$(tput setaf 1)Long read-based assembly has not been polished. $SAD $(tput sgr 0)"
else
    THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 2)Assembly is now polished with Racon. $THUMBS_UP $(tput sgr 0)"
fi

# Further polish the assembly using Medaka2
if [ ! -f "$output_dir/polished_assembly.fasta" ]; then
     SAD='\U1F641'; echo -e "$(tput setaf 1)Long read-based assembly has not $(tput sgr 0)"
     SAD='\U1F641'; echo -e "$(tput setaf 1)been polished with medaka2. $SAD $(tput sgr 0)"
else
    THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 2)Assembly is now polished with Medaka2. $THUMBS_UP $(tput sgr 0)"
fi
echo ""
echo "$(tput setaf 2)####################################################$(tput sgr 0)"
echo "$(tput setaf 2)##         Assemble short reads using SPAdes,     ##$(tput sgr 0)"
echo "$(tput setaf 2)##         together with assembled long reads     ##$(tput sgr 0)"
echo "$(tput setaf 2)####################################################$(tput sgr 0)"
echo ""
if [ ! -f "$output_dir/spades_assembly/scaffolds.fasta" ]; then
    SAD='\U1F641'; echo -e "$(tput setaf 1)Short read-based assembly is missing. $SAD $(tput sgr 0)"
    echo ""
else
    echo ""
    THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 2)Short reads are assembled with SPAdes. $THUMBS_UP $(tput sgr 0)"
    echo ""
    THUMBS_UP='\U1F44D'; echo -e "$(tput setaf 2)Job is done $THUMBS_UP $(tput sgr 0)"
fi