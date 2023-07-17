#! /usr/bin/env python3

import sys
import pandas as pd

input1 = sys.argv[1] # file extracted from somatic.vcf or pindel_SI.vcf
input2 = sys.argv[2] # file extracted from somatic_VEP_delheaders or pindel_VEP_delheaders
outputFile = sys.argv[3]

#To read the input files
df1 = pd.read_csv(input1, sep='\t')
df2 = pd.read_csv(input2, sep='\t')

# Convert 'POS' column to a common data type (e.g., str)
df1['POS'] = df1['POS'].astype(str)
df2['POS'] = df2['POS'].str.split('-').str[0]

# Merge based on 'CHROM' and 'POS'
merged = pd.merge(df1, df2, on=['CHROM', 'POS'], how='inner')
merged = merged.drop_duplicates() #to remove duplicates
merged.to_csv(outputFile, sep = '\t', header=True, index=False)
