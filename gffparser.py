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
    parser.add_argument('-f', '--input_format', help='gtf or gff3 (default: gff3)')

    args = parser.parse_args()

    f_in = sys.stdin
    f_out = sys.stdout
    input_format = 'gff3'

    if args.input_format:
        input_format = args.input_format
        if input_format not in ['gff3', 'gtf']:
            sys.exit('ERROR: %s extension is not supported' % input_format)

    if args.input:
        if not os.path.exists(args.input):
            sys.exit("ERROR: file %s doesn't exist" % args.input)
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

    coords = parseGFF(f_in, input_format)

    # print(len([x for x in coords]))
    coords2extb(coords, f_out)


    f_in.close()
    f_out.close()


array2str = lambda arr: ','.join(str(x) for x in arr)
# array2str = lambda arr: ','.join('%s' % x for x in arr)
str2array = lambda string: np.fromstring(string, sep=',', dtype=int)

def parseGFF(f_in, input_format):
    coords = {}
    for line in f_in:
        # stop parsing if fasta is found (assumes fasta will always be AFTER ALL gff fields
        if line.startswith('>'):
            break

        line = line.strip('\n')
        gff_line_parser(line, coords, file_format=input_format)
    return coords

def coords2extb(coords, f_out):
    for trans in coords:
        strand = coords[trans]['strand']
        exons, introns = calc_extb_features(coords[trans])
        print(trans, strand, len(exons), array2str(exons), array2str(introns), sep='\t', file=f_out)




def calc_extb_features(transcript_dict):
    starts = str2array(transcript_dict['exon_starts'])
    ends = str2array(transcript_dict['exon_ends'])
    strand = transcript_dict['strand']
    exons, introns = calcExInt(starts, ends, strand)

    gcoords = ';'.join('%d..%d' % (s, e) for s, e in zip(starts, ends))

    ## CDS
    if 'CDS_starts' in transcript_dict:
        phase = str2array(transcript_dict['phase'])
        cds_s = str2array(transcript_dict['CDS_starts'])
        cds_e = str2array(transcript_dict['CDS_ends'])

        cdsS, cdsE, cds = calcCDS(cds_s, cds_e, strand)
        cdsSE = '%d,%d' % (cdsS, cdsE)
    else:
        cdsSE = cds = phase = [np.nan]



    return exons, introns
    # return exons, cds

def calcExInt(starts, ends, strand):
    exons = ends - starts + 1
    introns = [np.nan]
    if len(exons) > 1:
        if strand == '-':
            introns = starts[:-1] - ends[1:] - 1
        else:
            introns = starts[1:] - ends[:-1] - 1
    return exons, introns


def calcCDS(cds_s, cds_e, strand):

    cds = cds_e - cds_s + 1
    if strand == '-':
        cds_s = cds_s[::-1]
        cds_e = cds_e[::-1]

    cdsS = cds_s[0]
    cdsE = cds_e[-1]

    # if strand == '-':
    #     utr3, utr5 = countUTR(cds_s, cds_e, cdsS, cdsE)
    # else:
    #     utr5, utr3 = countUTR(cds_s, cds_e, cdsS, cdsE)
    #
    #
    # cds = np.hstack([np.zeros(utr5), cds, np.zeros(utr3)])
    # cds = np.array(cds, dtype=int)
    return cdsS, cdsE, cds



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


# def attributesParser(field9, file_format='gff3'):
#     # pattern = re.compile(r';*(\w+)=([^;]+)')
#     if file_format == 'gff3':
#         pattern = re.compile(r'(\w+)=([^;]+)')
#         step = 3
#     elif file_format == 'gtf':
#         pattern = re.compile(r'[;\s"]*')
#         step = 2
#
#     attributes = re.split(pattern, field9)
#     att = dict((k, re.sub(r'^\w+:', '', v)) for k, v in (zip(attributes[1::step], attributes[2::step])))
#     # att = dict(zip(attributes[1::3], attributes[2::3]))
#     # print(attributes)
#     # print(att)
#     return att

def attributesParser(field9, file_format='gff3'):
    if file_format == 'gff3':
        pattern = re.compile(r'^\s*(\S+)\s*=\s*(.*)\s*$')
    if file_format == 'gtf':
        pattern = re.compile(r'^\s*(\S+)\s+\"([^\"]+)\"\s*')


    from html import unescape
    # html_parser = html.parser.HTMLParser()
    field9 = unescape(field9)

    attributes = {}
    atts = re.sub(';\s*$', '', field9)
    atts = atts.split(';')
    for att in atts:
        g = re.search(pattern, att)
        try:
            k, v = g.group(1,2)
            v = re.sub(r'^(transcript|gene):', '', v)
            attributes[k] = v
        except:
            sys.exit('PARSING ERROR: regex %s failed to capture attributes \nLINE: %s' % (str(pattern), field9))
    return attributes


def gff_line_parser(line, coords, file_format='gff3'):
    if not line.startswith('#'):
        fields = re.split(r'\t', line)

        if file_format == 'gff3' and re.match('transcript|mRNA', fields[2]):
        # if re.match('transcript|mRNA', fields[2]):
            att = attributesParser(fields[8]) #, file_format=file_format)
            # print(att)
            # if file_format == 'gff3':
            transc_id = 'ID'
            gene_id = 'Parent'
            # elif file_format == 'gtf':
            #     transc_id = 'transcript_id'
            #     gene_id = 'gene_id'

            transID = att[transc_id]
            coords[transID] = {'gene': att[gene_id],
                               'strand': fields[6]}

        # if fields[2] == 'exon':
        elif re.match('exon|CDS', fields[2]):
            if file_format == 'gff3':
                transc_id = 'Parent'
            elif file_format == 'gtf':
                transc_id = 'transcript_id'
                gene_id = 'gene_id'

            att = attributesParser(fields[8], file_format=file_format)

            transID = att[transc_id]

            # if transID not in coords:
            if file_format == 'gtf' and transID not in coords:
                coords[transID] = {'gene': att[gene_id],
                                 'strand': fields[6],}
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
