#!/usr/bin/env python3
# This script will check the pharmacogenomics markers in the final.concat.csv files

import sys
import re
import pandas as pd
import csv

input_xlsx_file = sys.argv[3]	# The first file is .xlsx containing the pharmacogenic markers to be mapped
sample = sys.argv[1]
filepath = sys.argv[2]
output = sys.argv[4]

map = dict()			# Dictionary for markers
map_counts = dict()		# Dictionary for the marker occurence
map_comments = dict()	# Dictionary for comments
xlsx_file_data = pd.read_excel(input_xlsx_file, engine='openpyxl', usecols = [0,1,2,3,4,5])	# Extracting the chrom,pos,rsid,ref and alt columns
for i in range(len(xlsx_file_data)):
	chromosome = xlsx_file_data.iloc[i, 0]
	position = xlsx_file_data.iloc[i, 1]
	ref = xlsx_file_data.iloc[i, 2]
	alt = xlsx_file_data.iloc[i, 3]
	rsid = xlsx_file_data.iloc[i, 4]
	comment = str (xlsx_file_data.iloc[i, 5]).replace(",", "&")
	comment = comment.lower()

	if re.search('[a-zA-Z]', str(chromosome)):
		chromosome = re.sub ("chr","", chromosome, flags = re.IGNORECASE)
		chromosome = re.sub ("X","23", chromosome, flags = re.IGNORECASE)
		chromosome = re.sub ("Y","y", chromosome, flags = re.IGNORECASE)
		
	id = str(chromosome) + ':' +''.join(str(position)) + ':' +''.join(str(ref)) + ':' +''.join(str(alt))
	# print ( chromosome, position)
	map[id] = rsid
	map_counts[id] = 0
	map_comments[id] = comment

if map:
	pass
else:
	print ("Exiting, nothing to be mapped!")
	quit ()

# Analyzing the csv file 
csv_file = filepath + sample + '.final.concat.csv'
outfile = open (output,'w')
#print ("Chr", "Pos" , "Ref allele", "Alt allele", "No_of_variant_callers", "Variant_callers", "REF_COUNT", "ALT_COUNT", "VAF", "PopFreqMax","Gene.refGene","rsid","COMMENTS", file=outfile, sep="," )
print ("Chr", "Pos" , "No_of_variant_callers", "Variant_callers", "PopFreqMax", "Gene", "REF_COUNT", "ALT_COUNT", "Ref allele", "Alt allele", "VAF", "rsid","COMMENTS", file=outfile, sep="," )
with open (csv_file,'r') as vcf:
	csv_handle = csv.reader(vcf)
	header = next (csv_handle)	# Removing the header
	for str_lines in csv_handle:
		# print (str_lines)
		vcf_chr = str_lines[0]
		vcf_pos = str_lines[1]
		# vcf_rsids = str_lines[3]
		vcf_ref = str_lines[3]
		vcf_alt = str_lines[4]
		vcf_variant_callers = str_lines[5]
		no_of_callers = len (vcf_variant_callers.split('|')) - 1
		vcf_REF_COUNT = str_lines[9]
		vcf_ALT_COUNT = str_lines[10]
		vcf_VAF = str_lines[11]
		vcf_PopFreqMax = str_lines[20]
		vcf_GenerefGene = str_lines[13]

		if re.search('[a-zA-Z]', str(vcf_chr)):
			vcf_chr = re.sub ("chr","", vcf_chr, flags = re.IGNORECASE)
			vcf_chr = re.sub ("X","23", vcf_chr, flags = re.IGNORECASE)
			vcf_chr = re.sub ("Y","y", vcf_chr, flags = re.IGNORECASE)

		vcf_id = str(vcf_chr) + ':' +''.join(vcf_pos) + ':' +''.join(vcf_ref) + ':' +''.join(vcf_alt)
		if vcf_id in map:
			map_counts[vcf_id] = map_counts[vcf_id] + 1
			#print (vcf_chr, vcf_pos, vcf_ref, vcf_alt, no_of_callers, vcf_variant_callers, vcf_REF_COUNT, vcf_ALT_COUNT, vcf_VAF, vcf_PopFreqMax, vcf_GenerefGene, map[vcf_id], map_comments[vcf_id], file=outfile, sep=",")
			print (vcf_chr, vcf_pos, no_of_callers, vcf_variant_callers, vcf_PopFreqMax, vcf_GenerefGene, vcf_REF_COUNT, vcf_ALT_COUNT, vcf_ref, vcf_alt, vcf_VAF, map[vcf_id], map_comments[vcf_id], file=outfile, sep=",")

# for ids in map:
# 	print (ids, map[ids], map_counts[ids], file=outfile)

outfile.close()
