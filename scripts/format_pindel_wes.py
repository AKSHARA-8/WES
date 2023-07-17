#! /usr/bin/env python3

import pandas as pd
import numpy as np
import os
import sys
import subprocess

args = sys.argv
bed = args[1]
pindel_vep = args[2]
outfile = args[3]

df=pd.read_csv(bed, sep="\t")
print(df)
if(os.stat(pindel_vep).st_size != 0):
	df.columns=["A","B","C","D","E"]
	df1=pd.read_csv(pindel_vep,sep="\t")
	df1.insert(loc=6, column="Ref_Count",value=0.0)
	df1.insert(loc=7, column="Alt_Count",value=0.0)
	df1.insert(loc =8,column = 'VAF',value=0.000)
	#df1.insert(loc =10,column = 'AR', value=0.000)
	#for i in range(len(df1)):
			#df1["Alt_Count"][i]= str(df1["Otherinfo1"][i]).split(",")[1]

	print(df1)
        #print(df,"llllllll")
        #print(df1["AF"])
	print(df1["AR"])
	print(df1["Ref_Count"])
	print(df1["Alt_Count"])
        #df1["values"].astype('float64', copy=False)
        #df["E"].astype('float64', copy=False)
        #df.astype({'E': 'float64'}).dtypes
        #df1.astype({'values': 'float64'}).dtypes
	print((df1.dtypes),(df.dtypes))


	for i in range(len(df)):
			for j in range(len(df1)):
					if((df1["POS"][j] in range(df["B"][i],df["C"][i])) and df1["CHROM"][j]==df["A"][i]):
							df1["Ref_Count"][j]=df["E"][i]+0.0
							#df1['AF'][j]=df1['Alt_Count'][j]/df1['Ref_Count'][j]
							#df1["AR"][j]=df1['AF'][j]/(1-df1['AF'][j])
							df1['VAF'][j]=df1['Alt_Count'][j]/df1['Ref_Count'][j]*100
        #df1['VAF'] = np.where(df1["Start"]==df["B"] & df1["Chr"]==df["A"],df1['values']/df['E'],0) #create a new column in df1 for price diff
	print(df1)
        #print(df1["AF"])
	print(df1["VAF"])
        #print(df1["AR"])
	print(df1["Ref_Count"])
	print(df1["Alt_Count"])
	df1.to_csv(outfile,sep = '\t', index=False)
else:
	subprocess.run(['touch', outfile], stdout=subprocess.DEVNULL)
