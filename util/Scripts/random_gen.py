#!/usr/bin/env python3

import argparse
import pandas as pd
import random

parser = argparse.ArgumentParser(description='Create a Covariate or Phenotype file with random values')

parser.add_argument('--ncol', dest='ncol', type=int, help='number of column expected in the outfile')

parser.add_argument('--ID', dest='Input', type=str, help='file with individuals')

parser.add_argument('--out', dest='output', type=str, help='output')

parser.add_argument('--prefix', dest='pre', type=str, help='prefix to column')

args = parser.parse_args()


In=pd.read_table(args.Input)


out= list()
header=["FID", "IID"]
for y in range(0,args.ncol):
	header.append(args.pre + str(y))
for i in range(0,len(In.index)):
    line = list()
    for a in range(0,args.ncol):
        rando=random.random()*5
        line.append(rando)
    out.append(line)
df=pd.DataFrame(data=out)
print(In.shape)
print(In.shape)
Q=pd.concat([In.iloc[: , :2], df], axis=1, ignore_index=True)
Q.columns = header

Q.to_csv(args.output, index=False, sep='\t',header=True)
