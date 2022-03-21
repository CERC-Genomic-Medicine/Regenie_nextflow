#!/usr/bin/env python3



from bioinfokit import analys, visuz
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import statsmodels.api as sm
import pylab as py

# load dataset as pandas dataframe
parser = argparse.ArgumentParser(description='manhattan plot')

parser.add_argument('--input', dest='input', type=str, help='input file for manhattan plot')


args = parser.parse_args()

df= pd.read_table(args.input,header=0,sep=" ",low_memory=False)

df['pvalue']=10**(-df['LOG10P'])
# create Manhattan plot with default parameters
visuz.marker.mhat(df=df, chr='CHROM',pv='pvalue')
# set parameter show=True, if you want view the image instead of saving
