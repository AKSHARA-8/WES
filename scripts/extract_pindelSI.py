#! /usr/bin/env python3

import pandas as pd
import sys

inputFile = sys.argv[1] #input file pindel_SI.vcf from pindel
outputFile =sys.argv[2]

with open(inputFile) as f:
    for line in f:
        if line.startswith('#CHROM'):
            columns = line.strip().split('\t')
            columns[0] = 'CHROM'  # Change the column name from '#CHROM' to 'CHROM'
            break

data = pd.read_csv(inputFile, sep='\t', comment='#', skiprows=45, names=columns)
extractedData = data[['CHROM', 'POS', 'REF', 'ALT']]

extractedData.to_csv(outputFile, sep='\t', header=True, index=False)
