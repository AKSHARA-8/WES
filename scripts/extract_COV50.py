#! /usr/bin/env python3

import pandas as pd
import numpy as np
import os
import sys

inputFile = sys.argv[1] #coverage.bed
outputFile = sys.argv[2]

Column = ['Chrom', 'Start' ,'End' ,'Median']
df = pd.read_csv(inputFile, sep ='\t',names=Column)

medianData  = df[df.iloc[:,-1] < 50]

medianData.to_csv(outputFile, sep ='\t', header=True, index=False)
