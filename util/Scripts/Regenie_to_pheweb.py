import gzip
import argparse


argparser = argparse.ArgumentParser(description = 'Script to plot either a Manhattan or Miami plot')
argparser.add_argument('-i', metavar = 'name', dest = 'input', type = str, required = True, help = 'input')
argparser.add_argument('-o', metavar = 'name', dest = 'output', type = str, required = True, help = 'output')




def process_and_write_compressed_file(input_file_path, output_file_path):
    with gzip.open(input_file_path, 'rt') as file:  # 'rt' mode for reading text
        # Open the output file
        with open(output_file_path, 'w') as outfile:
            # Read the header line and rename the columns
            header_line = file.readline().lower().strip()
            header = header_line.split()
            dic_header = {'chrom': 'chrom','genpos': 'pos', 'allele0': 'ref','allele1': 'alt','a1freq': 'af','n': 'num_samples','beta': 'beta', 'se': 'sebeta','log10p': 'pval', 'info':'r2'}
            new_header = [dic_header[item] if item in dic_header else item for item in header]
            pvalue_index=new_header.index('pval')
            # Write the new header to the output file
            outfile.write('\t'.join(new_header) + '\n')

            # Process each line
            for line in file:
                values = line.strip().split()
                # Convert the LOG10P value to p-value
                if values[pvalue_index]=='NA':
                    continue
                else :
                    values[pvalue_index] = str(10 ** -float(values[pvalue_index]))
                    if values[0] =='23':
                        values[0]='X'
                    # Write the new line to the output file
                    outfile.write('\t'.join(values) + '\n')
if __name__ == '__main__':
    args = argparser.parse_args()
    process_and_write_compressed_file(args.input, args.output)
