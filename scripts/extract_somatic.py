#! /usr/bin/env python3

import sys
import csv
import pandas as pd

vcf_file = sys.argv[1]
outputFile = sys.argv[2]

# Read VCF file into a DataFrame
vcf_df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None)
#using comment"#" it will read anything starting with'#' should be treated as comments && header=None indicates that the file doesn't contain column headers
#basically removing #info and latter specifilying the column name

#for extracting the specified first 5 columns
extracted_df = vcf_df.iloc[:, :5]

# to extract variantcaller tool


vlsfph_values = vcf_df[7].apply(lambda x: x.split('VLSFPH=')[1].split(';')[0].split(','))
Variant_callers = []
for values in vlsfph_values:
    callers = []
    for i, value in enumerate(values):
        if i < len(['VarScan2', 'LoFreq', 'Strelka', 'Freebayes', 'Platypus', 'Haplotypecaller']):
            if value == '1':
                callers.append(['VarScan2', 'LoFreq', 'Strelka', 'Freebayes', 'Platypus', 'Haplotypecaller'][i])
    Variant_callers.append(", ".join(callers))


extracted_df['Variant_callers'] = Variant_callers

#to extract ref_count and alt_count from CD4
format_column = vcf_df[9]
cd4_values = format_column.str.split(':').str[2].str.split(',')
ref_count = cd4_values.str[1]
alt_count = cd4_values.str[3]

extracted_df['ref_count'] = ref_count
extracted_df['alt_count'] = alt_count

vaf_values = format_column.str.split(':').str[-1]  #extracting vaf valu and converting it into percentage

vaf_values = pd.to_numeric(vaf_values) * 100

#Add the VAF% values to the extracted DataFrame
extracted_df['VAF%'] = vaf_values

#Set column names
extracted_df.columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT','Variant_callers', 'ref_count', 'alt_count', 'VAF%']

extracted_df.to_csv(outputFile, sep ='\t', header=True, index=False)
