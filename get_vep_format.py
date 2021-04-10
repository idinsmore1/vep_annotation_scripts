import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input_file', help = 'Input file')
parser.add_argument('--output_file', help = 'Output file name')
parser.add_argument('--chr', default='CHROM', help = 'Chromosome Column name. Default CHROM')
parser.add_argument('--id', default= 'ID', help = 'ID column. Default ID')
parser.add_argument('--ref', default='ALLELE0', help = 'Ref allele column. Default ALLELE0')
parser.add_argument('--alt', default='ALLELE1', help = 'Alt allele column. Default ALLELE1')
parser.add_argument('--pos', default='GENPOS', help = 'Genomic position column. Default GENPOS')
args = parser.parse_args()

def create_allele_col(df, allele0, allele1):
    df['Allele'] = df[f'{allele0}'] + '/' + df[f'{allele1}']

def get_vep_bp_and_alleles(alleles):
    split_char = '/'
    num_bp = list()
    new_alleles = list()
    for allele in alleles:
        ins = allele.partition(split_char)
        if len(ins[0]) > len(ins[2]):
            num_bp.append(len(ins[0]))
            new_alleles.append(ins[0] + '/-')
        elif len(ins[2]) > len(ins[0]):
            num_bp.append(1)
            new_alleles.append('-/' + ins[2][0])
        else:
            num_bp.append(0)
            new_alleles.append(allele)
    return num_bp, new_alleles

def update_positions(IDs, pos1, pos2, indel_len):
    for i in range(len(IDs)):
        if 'I' in IDs[i]:
            pos1[i] = pos1[i] + indel_len[i]
        if 'D' in IDs[i]:
            pos2[i] = pos2[i] + indel_len[i]

def get_vep_format_file(filepath, outfile, chr, ids, pos, a1, a2):
    results = pd.read_csv(f'{filepath}', sep = '\s+')
    create_allele_col(results, a1, a2)
    alleles = list(results['Allele'])
    num_bp, new_alleles = get_vep_bp_and_alleles(alleles)
    IDs = list(results[f'{ids}'])
    POS1 = list(results[f'{pos}'])
    POS2 = list(results[f'{pos}'])
    update_positions(IDs, POS1, POS2, num_bp)
    results['POS1'] = POS1
    results['POS2'] = POS2
    results['new_allele'] = new_alleles
    results['strand'] = '+'
    vep_out = results[[f'{chr}', 'POS1', 'POS2', 'new_allele', f'{ids}']]
    vep_out.to_csv(f'{outfile}',sep = '\t', index = False, header=False)

get_vep_format_file(args.input_file, args.output_file, args.chr, args.id, args.pos, args.ref, args.alt)