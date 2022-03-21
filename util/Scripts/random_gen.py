#!/usr/bin/env python3

import pandas as pd
import sys
from random import seed
from random import gauss
import sys, getopt


def main(argv):
  inputfile = ''
  outputfile = ''
  ncolumn = ''
  try:
    opts, args = getopt.getopt(argv,"i:o:n",["ifile=","ofile=","nb_column="])
  except getopt.GetoptError:
    print('test.py -i <inputfile> -o <outputfile> -n <number_column> ')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      print('test.py -i <inputfile> -o <outputfile> -n <number_column>')
      sys.exit()
    elif opt in ("-i", "--ifile"):
      inputfile = arg
    elif opt in ("-o", "--ofile"):
      outputfile = arg
    elif opt in ("-n", "--nb_column"):
      ncolumn = arg
  input=pd.read_table(inputfile)
  for x in range(2):
    colname=x
    rando=list()
    for y in range(len(input)):
      rando.append(gauss(0, 1))
    input[colname] = rando
  input.to_csv(outputfile, sep='\t')

if __name__ == "__main__":
   main(sys.argv[1:])