#!/bin/sh
# Number of contigs (rows of the first column)
num_contigs=$(wc -l < $1)

# Mean contig length (mean of column two)
mean_length=$(awk '{sum += $2} END {print sum/NR}' $1)

# Mean coverage of contigs (mean of column three)
mean_coverage=$(awk '{sum += $3} END {print sum/NR}' $1)

# Display the summary
echo -e "  #########################################"
echo -e "  ##       Long assembly information     ##"
echo -e "  #########################################"
echo "      Number of contigs: $(tput setaf 2)       $num_contigs$(tput sgr 0)"
echo "      Mean contig length: $(tput setaf 2)      $mean_length$(tput sgr 0)"
echo "      Mean coverage of contigs: $(tput setaf 2)$mean_coverage$(tput sgr 0)"