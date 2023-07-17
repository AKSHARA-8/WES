i!/usr/bin/bash

bamfile=$1
output_prefix=$2
bedfile=$3

source activate new_base

mosdepth -n -t 3 --fast-mode -m --by ${bedfile} ${output_prefix} ${bamfile}
gunzip ${output_prefix}.regions.bed.gz

conda deactivate
