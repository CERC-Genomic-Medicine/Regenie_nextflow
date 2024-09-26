#!/usr/bin/env python3

#
# Author: Vincent Chapdelaine <vincent.chapdelaine@mail.mcgill.ca>; Daniel Taliun <daniel.taliun@mcgill.ca>
# Version: 2.0
# Year: 2024
#

import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse


argparser = argparse.ArgumentParser(description = 'Script to plot either a Manhattan or Miami plot')
argparser.add_argument('-g', '--gwas', metavar = 'file', dest = 'in_gwas_files', type = str, nargs = '+', required = True, help = 'Compressed (gzip) GWAS result file(s) in Regenie format. Specify two files to create a Miami plot.')
argparser.add_argument('-t', '--title', metavar = 'name', dest = 'in_gwas_titles', type = str, nargs = '+', required = True, help = 'Title(s) for GWAS plots. Specify two titles if you are creating Miami plot e.g.: -t "Type 2 Diabetes in Males" "Type 2 Diabetes in Females"')
argparser.add_argument('-f', '--min-af', metavar = 'float', dest = 'min_af', type = float, required = False, default = 0.0, help = 'Threshold for the minimal alternate allele frequency (AF). Default: 0.0')
argparser.add_argument('-q', '--min-imp-quality', metavar = 'float', dest = 'min_info', type = float, required = False, default = 0.0, help = 'Threshold for the minimal imputation quality (INFO or Rsq). Defailt: 0.0') 
argparser.add_argument('-p', '--max-log-pvalue', metavar = 'float', dest = 'max_log10p', type = float, required = False, default = 100, help = 'Maximal allowed -log_10(pvalue). If this threshold is exceeded, then the -log_10(pvalue) is set to this value. This option is useful only for the visualization purposes to avoid squeezed Manhattan plots due to very extreme P-values. Default: 100.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_filename', type = str, required = True, help = 'Output file name. The ".png" suffix will be appended automotically.')


CHROMOSOME_CODES = dict(
        [(str(i), np.uint8(i)) for i in range(1, 24)] + 
        [(f'chr{i}', np.uint8(i)) for i in range(1, 23)] + 
        [('X', np.uint8(23)), ('chrX', np.uint8(23))]
    )


SIGNIFICANCE_THRESHOLD = -np.log10(0.05 * 1e-6)


def read_regenie_continuous(filename, min_af, min_info, max_log10p):
    required_columns = ['CHROM', 'GENPOS', 'A1FREQ', 'LOG10P']
    n_total_entries = 0
    n_good_entries = 0

    min_af = np.float32(min_af)
    min_info = np.float32(min_info)
    max_log10p = np.float16(max_log10p)
    
    with gzip.open(filename, 'rt') as ifile:
        print(f'Scanning {filename}...')
        header = ifile.readline().rstrip().split(' ')
        if any(c not in header for  c in required_columns):
            raise Exception('Required CHROM, GENPOS, A1FREQ, and LOG10P columns are missing from the GWAS header.')
        chrom_idx = header.index('CHROM')
        genpos_idx = header.index('GENPOS')
        a1freq_idx = header.index('A1FREQ')
        log10p_idx = header.index('LOG10P')
        if 'INFO' in header:
            info_idx = header.index('INFO')
        else:
            info_idx = None
        for n, line in enumerate(ifile, 1):
            if n % 1000000 == 0:
                print(f'\r{n} records scanned...', end = '', flush = True)
        print(f'\rDone. {n} records detected.\t\t')
        
        chromosomes = np.zeros(n, dtype = np.uint8)
        positions = np.zeros(n, dtype = np.uint32)
        minus_log10_pvalues = np.zeros(n, dtype = np.float16) # Use float16 because we do not need precision when visualizing.

    with gzip.open(filename, 'rt') as ifile:
        print(f'Loading {filename}...')
        i = 0
        ifile.readline() # skip header
        for line in ifile:
            fields = line.rstrip().split(' ')
            if fields[chrom_idx] not in CHROMOSOME_CODES:
                continue
            if np.float32(fields[a1freq_idx]) < min_af:
                continue
            if info_idx is not None and np.float32(fields[info_idx]) < min_info:
                continue
            chromosomes[i] = CHROMOSOME_CODES[fields[chrom_idx]]
            positions[i] = np.uint32(fields[genpos_idx])
            minus_log10_pvalue = np.float16(fields[log10p_idx])
            minus_log10_pvalues[i] = min(minus_log10_pvalue, max_log10p)
            i += 1
            if i % 1000000 == 0:
                print(f'\r{i} records loaded...', end = '', flush = True)
        print(f'\rDone. {i} records loaded (AF => {min_af}; INFO >= {min_info}).\t\t')
    if i < chromosomes.size:
        chromosomes.resize(i)
        positions.resize(i)
        minus_log10_pvalues.resize(i)
    return chromosomes, positions, minus_log10_pvalues


def manhattan(Fig, file, inverted, title, min_af, min_info):
    chromosomes, positions, minus_log10_pvalues = read_regenie_continuous(file, min_af, min_info, 100)

    df_gwas = pd.DataFrame(data = {'chromosome': chromosomes, 'position': positions, 'minus_log10_pvalue': minus_log10_pvalues}, copy = False)

    print(f'Plotting {file}... ', end = '', flush = True) 
    ax = Fig
    colors = ['black', 'gray']
    x_labels = []
    x_labels_pos = []
    position_offset = 0
    x_start = None
    x_stop = None
    for i in range(1, 24):
        df_gwas_chrom = df_gwas[df_gwas.chromosome == i]
        if len(df_gwas_chrom) == 0:
            continue
        if x_start is None:
            x_start = df_gwas_chrom.position.min()
        chrom_length = df_gwas_chrom.position.max()
        ax.scatter(x = df_gwas_chrom.position + position_offset, y = df_gwas_chrom.minus_log10_pvalue, color = colors[i % len(colors)], s = 9)
        x_labels.append(f'{i}' if i < 23 else 'X')
        x_labels_pos.append(position_offset + chrom_length / 2)
        position_offset += chrom_length
    x_stop = position_offset

    # set X-axis and Y-axis ticks and tick labels
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels, fontsize = 16)
    ax.tick_params(axis='y', labelsize = 16)

    # set axis limits
    ax.set_xlim([x_start, x_stop])
    ax.set_ylim(bottom = 0)

    # set and adjust axis labels depending on the plot type
    ax.set_ylabel('-'+r'$\log_{10}$' +'(P-value)', fontsize = 20)
    if inverted :
        ax.invert_yaxis()
        plt.tick_params(
            axis = 'x',          # changes apply to the x-axis
            which = 'both',      # both major and minor ticks are affected
            bottom = False,      # ticks along the bottom edge are off
            top = False,         # ticks along the top edge are off
            labelbottom = False)
        ax.set_title(title, fontsize = 24, y = 0, pad = -24)
    else :
        ax.set_xlabel('Chromosome', fontsize = 20)
        ax.set_title(title, fontsize = 24)
    
    # draw significance line
    ax.axhline(y=SIGNIFICANCE_THRESHOLD, color='r', linestyle='-')

    print('Done.')


if __name__ == '__main__':
    args = argparser.parse_args()
    
    if (len(args.in_gwas_files) != len(args.in_gwas_titles)):
        raise Exception('Number of titles must be equal to the number of input GWAS files')

    if (len(args.in_gwas_files) == 2):
        fig, axs = plt.subplots(2, 1, constrained_layout=True,figsize=(32,15))
        manhattan(axs[0], args.in_gwas_files[0], False, args.in_gwas_titles[0], args.min_af, args.min_info)
        manhattan(axs[1], args.in_gwas_files[1], True, args.in_gwas_titles[1], args.min_af, args.min_info)
        fig.savefig(f'{args.out_filename}.png', dpi=300)
    elif (len(args.in_gwas_files) == 1):
        fig, axs = plt.subplots(1, 1, constrained_layout=True,figsize=(32,8))
        manhattan(axs, args.in_gwas_files[0], False, args.in_gwas_titles[0], args.min_af, args.min_info)
        fig.savefig(f'{args.out_filename}.png', dpi=300)
    else:
        raise Exception("Maximum two GWAS files are allowed in input.")
