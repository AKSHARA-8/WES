#! /usr/bin/env python3

import pandas as pd
import os, sys
import re

args = sys.argv
sample = args[1]
filepath = args[2]
outfile = args[3]
#cava_path = args
mosdepth_path = args[4]
mosdepthCov_path = args[5]
mosCOV50_path =args[6]
#pharma_marker_path = args[5]
pindel_path = args[7]
somaticseq_path = args[8]
#cnvkit_path = args[7]
#pharma_marker_path = args[8]

#csvfilenames=[filepath+sample+'.final.concat.csv',cava_path+sample+'.cava.csv',pindel_path,coverview_path,filepath+sample+'.artefacts.csv',cnvkit_path,pharma_marker_path]
#csvfilenames=[cava_path+sample+'.cava.csv',pindel_path,coverview_path,filepath+sample+'.artefacts.csv',cnvkit_path,pharma_marker_path]

csvfilenames = [sample +'.cava.csv']

mosdepth_df = pd.read_csv(mosdepth_path, sep='\t')
mosdepth_csv_dir = os.path.dirname(mosdepth_path)
os.makedirs(mosdepth_csv_dir, exist_ok=True)
mosdepth_csv_path = os.path.join(mosdepth_csv_dir, sample + '_mosdepth.csv')
mosdepth_df.to_csv(mosdepth_csv_path, index=False)
csvfilenames.append(mosdepth_csv_path)

moscovColumn = ['Chrom', 'Start' ,'End' ,'Median']
mosdepthCov_df = pd.read_csv(mosdepthCov_path, sep='\t',names=moscovColumn)
mosdepthCov_csv_dir = os.path.dirname(mosdepthCov_path)
os.makedirs(mosdepthCov_csv_dir, exist_ok=True)
mosdepthCov_csv_path = os.path.join(mosdepthCov_csv_dir, sample + '_mosdepthCov.csv')
mosdepthCov_df.to_csv(mosdepthCov_csv_path, index=False)
csvfilenames.append(mosdepthCov_csv_path)

mosCOV50_df = pd.read_csv(mosCOV50_path, sep='\t')
mosCOV50_csv_dir = os.path.dirname(mosCOV50_path)
os.makedirs(mosCOV50_csv_dir, exist_ok=True)
mosCOV50_csv_path = os.path.join(mosCOV50_csv_dir, sample + '_mosCOV50.csv')
mosCOV50_df.to_csv(mosCOV50_csv_path, index=False)
csvfilenames.append(mosCOV50_csv_path)

pindel_df = pd.read_csv(pindel_path, sep='\t')
pindel_csv_dir = os.path.dirname(pindel_path)
os.makedirs(pindel_csv_dir, exist_ok=True)
pindel_csv_path = os.path.join(pindel_csv_dir, sample + '_pindel.csv')
pindel_df.to_csv(pindel_csv_path, index=False)
csvfilenames.append(pindel_csv_path)

somaticseq_df = pd.read_csv(somaticseq_path, sep='\t')
somaticseq_csv_dir = os.path.dirname(somaticseq_path)
os.makedirs(somaticseq_csv_dir, exist_ok=True)
somaticseq_csv_path = os.path.join(somaticseq_csv_dir, sample + '_somaticseq.csv')
somaticseq_df.to_csv(somaticseq_csv_path, index=False)
csvfilenames.append(somaticseq_csv_path)

writer = pd.ExcelWriter(outfile)
for csvfilename in csvfilenames:
    sheetname = os.path.split(csvfilename)[1]
    df = pd.read_csv(csvfilename)
    print('process file:', csvfilename, 'shape:', df.shape)
    df.to_excel(writer, sheet_name=os.path.splitext(sheetname)[0], index=False)
writer.save()

