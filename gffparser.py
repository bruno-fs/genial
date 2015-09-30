#!/usr/bin/env python
__author__ = 'bruno'

import re
import sys
import os
import argparse as argp
import numpy as np
import gzip

def main():
    parser = argp.ArgumentParser(description="extract blockSizes from a gff file")
    parser.add_argument('-i', '--input', help='input file (default is stdin)')
    parser.add_argument('-o', '--output', help='output file (default is stdout)')
    # parser.add_argument('-f', '--input_format', help='gtf or gff3 (default: gff3)')

    args = parser.parse_args()

    f_in = sys.stdin
    f_out = sys.stdout
    # input_format = 'gff3'

    if args.input:
        if not os.path.exists(args.input):
            sys.exit("ERROR: gff file %s doesn't exist" % args.input)
        # elif
        elif re.search(r'\.gz$', args.input):
            f_in = gzip.open(args.input, 'rt')
            # print(f_in.readline())
            # sys.exit('foi =)')
        else:
            f_in = open(args.input, 'r')
    if args.output:
        if os.path.exists(args.output):
            sys.exit('ERROR: %s already exists!!!' % args.output)
        else:
            f_out = open(args.output, 'w')
    # if args.input_format:
    #     input_format = args.input_format

    coords = {}

    for line in f_in:
        # stop parsing if fasta is found (assumes fasta will always be AFTER ALL gff fields
        if line.startswith('>'):
            break

        line = line.strip('\n')
        parseGFF(line, coords)

    for trans in coords:
        strand = coords[trans]['strand']
        exons, introns = calcExIntCds(coords[trans])
        print(trans, strand, len(exons), array2str(exons), array2str(introns), sep='\t', file=f_out)
    f_in.close()
    f_out.close()


array2str = lambda arr: ','.join(str(x) for x in arr)
str2array = lambda string: np.fromstring(string, sep=',', dtype=int)


def calcExIntCds(transcript_dict):
    starts = np.fromstring(transcript_dict['exon_starts'], sep=',', dtype=int)
    ends = np.fromstring(transcript_dict['exon_ends'], sep=',', dtype=int)
    strand = transcript_dict['strand']
    exons, introns = calcExInt(starts, ends, strand)



    return exons, introns


def calcExInt(starts, ends, strand):
    exons = ends - starts + 1
    introns = [np.nan]
    if len(exons) > 1:
        if strand == '-':
            introns = starts[:-1] - ends[1:] - 1
        else:
            introns = starts[1:] - ends[:-1] - 1
    return exons, introns



# def calcCDS(trans_dict):
#     if 'CDS_starts' not in trans_dict:
#         return np.nan, np.nan
#     else:




def countUTR(s, e, cdsS, cdsE):

    utr5 = 0
    utr3 = 0

    for x in e:
        if x > cdsS:
            break
        utr5 += 1

    for x in reversed(s):
        if x < cdsE:
            break
        utr3 += 1

    return utr5, utr3


def attributesParser(field9, file_format='gff3'):
    # pattern = re.compile(r';*(\w+)=([^;]+)')
    if file_format == 'gff3':
        pattern = re.compile(r'(\w+)=([^;]+)')
    elif file_format == 'gtf':
        pattern = re.compile(r'[;"\s]+')
    # else:
    #     sys.exit('ERROR: %s is not supported (only gtf and gff3 are valid)' % file_format)
    attributes = re.split(pattern, field9)
    att = dict((k, re.sub(r'^\w+:', '', v)) for k, v in (zip(attributes[1::3], attributes[2::3])))
    # att = dict(zip(attributes[1::3], attributes[2::3]))
    # print(attributes)
    return att

def parseGFF(line, coords, file_format='gff'):
    if not line.startswith('#'):
        fields = re.split(r'\t', line)
        # if file_format == 'gff':
        #     id = 'ID'
        #     parent = 'Parent'
        # if file_format == 'gtf':
        #     id = 'gene_id'
        #     parent = 'transcript_id'


        if re.match('transcript|mRNA', fields[2]):
            att = attributesParser(fields[8])
            # print(att)
            transID = att['ID']
            coords[transID] = {'gene': att['Parent'],
                               'strand': fields[6]}

        # if fields[2] == 'exon':
        if re.match('exon|CDS', fields[2]):
            att = attributesParser(fields[8])
            # if 'Parent' in att:
                # transID = '{}::{}'.format(att['Parent'], fields[2])
            transID = att['Parent']
            if transID in coords:
                starts = fields[2] + '_starts'
                ends = fields[2] + '_ends'
                if starts not in coords[transID]:
                    coords[transID][starts] = str(fields[3])
                    coords[transID][ends] = str(fields[4])
                    if fields[2] == 'CDS':
                        coords[transID]['phase'] = str(fields[7])

                else:
                    coords[transID][starts] += ',%s' % str(fields[3])
                    coords[transID][ends] += ',%s' % str(fields[4])
                    if fields[2] == 'CDS':
                        coords[transID]['phase'] += ',%s' % str(fields[7])


if __name__ == '__main__':

    # Deal with errno 32 (broken pipe)
    # http://stackoverflow.com/a/16865106
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

    main()

