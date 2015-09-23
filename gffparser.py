#!/usr/bin/env python
__author__ = 'bruno'

import re
import sys
import os
import argparse as argp
import numpy as np


def main():
    parser = argp.ArgumentParser(description="extrai os blockSizes de um gff3")
    parser.add_argument('-i', '--input', help='input file (default is stdin)')
    parser.add_argument('-o', '--output', help='output file (default is stdout)')
    # parser.add_argument('-f', '--input_format', help='gtf or gff3 (default: gff3)')

    args = parser.parse_args()

    f_in = sys.stdin
    f_out = sys.stdout
    # input_format = 'gff3'

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
    # if args.input_format:
    #     input_format = args.input_format

    coords = {}

    for line in f_in:
        line = line.strip('\n')
        parseGFF(line, coords)

    for trans in coords:
        strand = coords[trans]['strand']
        exons, introns = calcExInt(coords[trans])
        print(trans, strand, array2str(exons), array2str(introns), sep='\t', file=f_out)
    f_in.close()
    f_out.close()


array2str = lambda arr: ','.join(str(x) for x in arr)
str2array = lambda string: np.fromstring(string, sep=',', dtype=int)


def calcExInt(transcriptDict):
    starts = np.fromstring(transcriptDict['starts'], sep=',', dtype=int)
    ends = np.fromstring(transcriptDict['ends'], sep=',', dtype=int)
    exons = ends - starts + 1
    strand = transcriptDict['strand']
    introns = np.array([np.nan])
    if len(exons) > 1:
        if strand == '-':
            introns = starts[:-1] - ends[1:] - 1
        else:
            introns = starts[1:] - ends[:-1] - 1

    return exons, introns


def attributesParser(field9, file_format='gff3'):
    # pattern = re.compile(r'(\w+)=([^;]+)')
    if file_format == 'gff3':
        pattern = re.compile(r'(\w+)=([^;]+)')
    elif file_format == 'gtf':
        pattern = re.compile(r'[;"\s]+')
    else:
        sys.exit('ERROR: %s is not supported (only gtf and gff3 are valid)' % file_format)
    attributes = re.split(pattern, field9)
    # att = dict((k.upper(), v) for k, v in zip(attributes[1::2], attributes[2::2]))
    att = dict(zip(attributes[1::2], attributes[2::2]))

    return att

def parseGFF(line, coords):
    if not line.startswith('#'):
        fields = re.split(r'\t', line)
        if fields[2] == 'exon':
            att = attributesParser(fields[8])
            if 'Parent' in att:
                # transID = '{}::{}'.format(att['Parent'], fields[2])
                transID = att['Parent']
                if transID not in coords:
                    coords[transID] = {'starts': str(fields[3]),
                                       'ends': str(fields[4]),
                                       'strand': fields[6]
                                       }

                else:
                    coords[transID]['starts'] += ',%s' % str(fields[3])
                    coords[transID]['ends'] += ',%s' % str(fields[4])


        # else:
        #     # print(att)
        #     pass
        #     # sys.exit('ERRO: GFF mal formatado -> a seguinte linha nao tem atributo "Parent"\n%s' % line)


if __name__ == '__main__':
    main()

    # Deal with errno 32 (broken pipe)
    # http://stackoverflow.com/a/16865106
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

