#! /usr/bin/env python3

import sys
import pandas as pd

input1 = sys.argv[1] # file extracted from somatic.vcf
input2 = sys.argv[2] # file extracted from somatic_VEP_delheaders
outputFile = sys.argv[3]

#To read the input files
df1 = pd.read_csv(input1, sep='\t')
df2 = pd.read_csv(input2, sep='\t')

# Convert 'POS' column to a common data type (e.g., str)
df1['POS'] = df1['POS'].astype(str)
df2['POS'] = df2['POS'].astype(str)

# Merge based on 'CHROM' and 'POS'
#output = pd.merge(df1, df2, on=['CHROM', 'POS'], how='inner')
#output.to_csv(outputFile, sep = '\t', header=True, index=False)

# Merge based on 'CHROM' and 'POS'
merged = pd.merge(df1, df2, on=['CHROM', 'POS'], how='inner')

# Group by 'CHROM' and 'POS', and combine 'SYMBOL' entries
grouped = merged.groupby(['CHROM', 'POS']).agg({'SYMBOL': lambda x: ', '.join(x)}).reset_index()

# Merge grouped data with other columns
output = pd.merge(grouped, merged.drop('SYMBOL', axis=1), on=['CHROM', 'POS'], how='inner')
#output['VAF'] = output['VAF']*100 # Converting VAF TO VAF%
output = output.drop_duplicates() #to remove duplicates

extractedData = output[['CHROM',
						'POS' ,
						'REF',
						'ALT',
						'SYMBOL',
						'ID',
						'ref_count',
						'alt_count',
						'VAF%',
						'AF',
						'MAX_AF_POPS',
						'ALLELE_NUM',
						'Variant_callers',
						'VARIANT_CLASS',
						'CLIN_SIG',
						'Feature',
						'Feature_type',
						'Consequence',
						'cDNA_position',
						'CDS_position',
						'Protein_position',
						'Amino_acids',
						'Codons',
						'ENSP']]
extractedData.to_csv(outputFile, sep = '\t', header=True, index=False)
