#!/usr/bin/env python3
## Goal Create a file filled with random gaussian covariate or phenotype compatible with Regenie
## requieres argparse pandas and random

import argparse
import pandas as pd
import random

parser = argparse.ArgumentParser(description='Create a Covariate or Phenotype file with random values')

parser.add_argument('--ncol', dest='ncol', type=int, help='number of column expected in the outfile')

parser.add_argument('--ID', dest='Input', type=str, help='file with individuals') ## FID and IID column of .psam or the like

parser.add_argument('--out', dest='output', default='output', type=str, help='output')

parser.add_argument('--prefix', dest='pre',default='', type=str, help='prefix to each column')

parser.add_argument('--binary', dest='binary', action='store_true', required = False, help = 'Flag to produce binary product instead of continuous')


args = parser.parse_args()

with open(args.Input,'rt') as f, open(args.output,'wt') as out:
    header=["FID", "IID"]
    for y in range(0,args.ncol):
        header.append(args.pre + str(y))
    out.write('\t'.join(header)+'\n')
    for line in f :
        record = line.strip()
        if ' ' in record :
            record = record.split(" ")[0:2]
        if '\t' in record:
            record = record.split('\t')[0:2]
        if record[0] == 'FID' or record[0] == '#FID':
            continue
        for a in range(0,args.ncol):
            rando=random.random()*5
            if args.binary:
                rando = '1' if rando >2.5 else '0'
            record.append(str(rando))
        out.write('\t'.join(record)+'\n')
