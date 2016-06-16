__author__ = 'bruno'

import re
import sys
import os
import argparse as argp
import numpy as np
from Bio import SeqIO


def filterFasta(fastaIDs, fasta_in, fasta_out):
    """
    Filter selected set of sequences from given fasta

    Parameters
    ----------
    fastaIDs  : set containing desired IDs
    fasta_in  : path to input fasta
    fasta_out : path to output fasta

    """
    with open(fasta_out, 'w') as f_out:
        n = 0
        #         with gzip.open(fasta_in,'rt') as f_in:
        with open(fasta_in) as f_in:
            for seq in SeqIO.parse(f_in, 'fasta'):
                if seq.id in fastaIDs:
                    SeqIO.write(seq, f_out, 'fasta')
                    n += 1

    print('found {} sequences from {} requested ({} %)'.format(n, len(fastaIDs), round(100 * n / len(fastaIDs))))


def main():
    parser = argp.ArgumentParser(description="description")
    parser.add_argument('-i', '--seq_ids', help='file with desired seq_ids separated by line breaks or blank space')
    parser.add_argument('-fa_in', '--fasta_input', help='input fasta')
    parser.add_argument('-fa_out', '--fasta_output', help='output fasta')

    args = parser.parse_args()

    if args.fasta_input:
        if not os.path.exists(args.input):
            raise FileNotFoundError("ERROR: %s doesn't exist" % args.input)
        else:
            f_in = open(args.input, 'r')
    else:
        sys.exit('--fasta_input is required!')

    if args.fasta_output:
        if os.path.exists(args.output):
            raise FileExistsError('ERROR: %s already exists!!!' % args.output)
        else:
            f_out = open(args.output, 'w')
    else:
        sys.exit('--fasta_output is required!')


    f_in.close()
    f_out.close()


if __name__ == '__main__':
    main()
