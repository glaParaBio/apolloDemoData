#!/usr/bin/env python3

import argparse
import sys
import re
import pandas

def parse_cigar(cigar):
    cigar_list = []
    v = ''
    op = ''
    for x in cigar:
        if x.isdigit():
            v += x
        else:
            cigar_list.append((x, int(v)))
            v = ''
    return cigar_list

parser = argparse.ArgumentParser(description='Convert PAF file from miniprot from protein coordinates to genome coordinates')
parser.add_argument('--paf', '-p', type= str, help= 'PAF file typically from miniprot [%(default)s]', default='-')
parser.add_argument('--protein-locations', '-l', type= str, help= 'Tab delimited file giving the genomic start of each protein. Start in 1-based coordinates [%(default)s]', default='-')
parser.add_argument('--genome', '-g', type= str, help= 'Tab delimited file of chromosome sizes [%(default)s]', default='-')
parser.add_argument('--version', action='version', version='%(prog)s 0.1.0')

if __name__ == '__main__':
    args = parser.parse_args()

    paf = pandas.read_csv(args.paf, sep='\t', comment='#', header=None)
    plocs = pandas.read_csv(args.protein_locations, sep='\t')
    genome = pandas.read_csv(args.genome, sep='\t', header=None, usecols=[0, 1], names=['chrom', 'chrom_length'])

    paf.rename(columns={0: 'protein_id', 1: 'qlen', 2: 'qstart', 4: 'strand', 5: 'genome_chrom', 6: 'genome_length', 7: 'genome_start', 8: 'genome_end'}, inplace=True)

    gpaf = pandas.merge(paf, plocs, on='protein_id') 
    gpaf = pandas.merge(gpaf, genome, on='chrom')
    gpaf.rename(columns={'chrom': 'protein_chrom', 'chrom_length':'protein_chrom_length'}, inplace=True)
    gpaf['protein_start'] = gpaf['protein_start'] - 1

    for idx,row in gpaf.iterrows():
        cigar = None
        for x in row[12:]:
            if x.startswith('cg:Z:'):
                cg = re.sub('^cg:Z:', '', x)
                cigar = parse_cigar(cg)
                break
        if cigar is None:
            raise Exception('Cigar string not found in paf record:\n%s' % row)
        """
        From https://github.com/lh3/miniprot/blob/master/miniprot.1:
        M	Alignment match. Consuming n*3 nucleotides and n amino acids
        D	Delection. Consuming n*3 nucleotides
        F	Frameshift deletion. Consuming n nucleotides
        N	Phase-0 intron. Consuming n nucleotides
        G	Frameshift match. Consuming n nucleotides and 1 amino acid
        U	Phase-1 intron. Consuming n nucleotides and 1 amino acid
        V	Phase-2 intron. Consuming n nucleotides and 1 amino acid
        I	Insertion. Consuming n amino acids
        """
        protein_pos = row.protein_start
        genome_pos = row.genome_start
        for op,val in cigar:
            if op == 'M':
                protein_start = protein_pos
                if row.strand == '-':
                    protein_start += 3
                protein_end = protein_start + val * 3
                protein_pos = protein_end
                genome_start = genome_pos
                genome_end = genome_start + val * 3
                genome_pos = genome_end
                outpaf = [row.protein_chrom, row.protein_chrom_length, protein_start, protein_end, row.strand, row.genome_chrom, row.genome_length, genome_start, genome_end, genome_end - genome_start, genome_end - genome_start, 0, row.protein_id, cg]
                print('\t'.join([str(x) for x in outpaf]))
            elif op == 'D':
                genome_pos += val * 3
            elif op in ('F', 'N'):
                genome_pos += val
            elif op == 'G':
                protein_start = protein_pos
                if row.strand == '-':
                    protein_start += 3
                protein_end = protein_start + 3
                protein_pos = protein_end
                genome_start = genome_pos
                genome_end = genome_start + val
                genome_pos = genome_end
                outpaf = [row.protein_chrom, row.protein_chrom_length, protein_start, protein_end, row.strand, row.genome_chrom, row.genome_length, genome_start, genome_end, genome_end - genome_start, genome_end - genome_start, 0, row.protein_id, cg]
                print('\t'.join([str(x) for x in outpaf]))
            elif op in ('U', 'V'):
                protein_pos += 3
                genome_pos += val
            elif op == 'I':
                protein_pos += val * 3
            else:
                raise Exception('Unexpected cigar operator: %s' % op)
