#!/usr/bin/env python3
## Goal Create a file filled with random gaussian covariate or phenotype compatible with Regenie
## requieres argparse pandas and random

import argparse
import pandas as pd
import random

parser = argparse.ArgumentParser(description='Create a Covariate or Phenotype file with random values')

parser.add_argument('--ncol', dest='ncol', type=int, help='number of column expected in the outfile')

parser.add_argument('--ID', dest='Input', type=str, help='file with individuals') ## FID and IID column of .psam or the like

parser.add_argument('--out', dest='output', type=str, help='output')

parser.add_argument('--prefix', dest='pre', type=str, help='prefix to each column')

args = parser.parse_args()


In=pd.read_table(args.Input)

# Create header
out= list()
header=["FID", "IID"]
for y in range(0,args.ncol):
	header.append(args.pre + str(y))
#fill each line
for i in range(0,len(In.index)):
    line = list()
    for a in range(0,args.ncol):
        rando=random.random()*5
        line.append(rando)
    out.append(line)
df=pd.DataFrame(data=out)
Q=pd.concat([In.iloc[: , :2], df], axis=1, ignore_index=True)
Q.columns = header
#write ouput
Q.to_csv(args.output, index=False, sep='\t',header=True)
