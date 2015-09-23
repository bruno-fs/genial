#!/usr/bin/env python
__author__ = 'bruno'


import re
import sys
import os
import argparse as argp
import numpy as np
from gffparser import str2array



def main():
    parser = argp.ArgumentParser(description="cds / ex")
    parser.add_argument('-i','--input', help='input file (default is stdin)')
    parser.add_argument('-o', '--output', help='output file (default is stdout)')
    parser.add_argument('-c', '--column', help='coluna que contem os blockSizes')

    args = parser.parse_args()

    f_in = sys.stdin
    f_out = sys.stdout
    col = 2

    if args.input:
        if not os.path.exists(args.input):
            sys.exit('ERRO: gff %s não existe' % args.input)
        else:
            f_in = open(args.input, 'r')
    if args.output:
        if os.path.exists(args.output):
            sys.exit('ERRO: %s já existe!!!' % args.output)
        else:
            f_out = open(args.output, 'w')
    if args.column:
        col = args.column -1

    for line in f_in:
        line = line.strip('\n')
        f = line.split('\t')
        cds = str2array(f[col])

        mean = sum(cds)/len(cds)
        mean = round(mean, 2)
        print(f[0], f[1], f[2], mean, sep='\t', file=f_out)


    f_in.close()
    f_out.close()



if __name__ == '__main__':
    main()