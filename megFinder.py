#!/usr/bin/env python
__author__ = 'bruno'


import re
import sys
import os
import argparse as argp
import numpy as np
from gffparser import str2array



def main():
    parser = argp.ArgumentParser(description="busca genes com arquitetura de meg em arquivos bed e similares")
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
        field = line.split('\t')
        exons = str2array(field[col])
        if ismeg(exons):
            print(line, file=f_out)

    f_in.close()
    f_out.close()

def countMicExTandem(exons, MicEx=36, simmetric=True):
    """return the max sequence $\mu$Ex in tandem on the given array (of exons)
    if simmetric is set to True, the $\mu$Ex must also be simmetric"""
    c_max = c = 0

    def isSimmetric(x):
        if x%3 == 0:
            return True
        else:
            return False

    for Ex in exons:
        if Ex <= MicEx:
            c += 1
            if simmetric and not isSimmetric(Ex):
                c = 0
        else:
            c = 0
        if c > c_max:
            c_max = c
    return c_max

def ismeg(x):
    intEx = x[1:-1]
    if len(intEx) >=3:
        # if np.median(intEx) <= 36:
        if countMicExTandem(intEx, MicEx=42) >= 4:
            return True



if __name__ == '__main__':
    main()