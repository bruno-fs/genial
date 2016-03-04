#!/usr/bin/env python3

__author__ = 'bruno'

import argparse as argp
import gzip
import os
import re
import sys


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

        elif re.search(r'\.gz$', args.input):
            f_in = gzip.open(args.input, 'rt')

        else:
            f_in = open(args.input, 'r')
    if args.output:
        if os.path.exists(args.output):
            sys.exit('ERROR: %s already exists!!!' % args.output)
        else:
            f_out = open(args.output, 'w')

    coords = GffDict(f_in, input_format)

    coords2extb(coords, f_out)

    f_in.close()
    f_out.close()


from classes import *


class GffDict(dict):
    def __init__(self, file_handle, file_format, **kwargs):
        # self.gff_dict = {}
        super().__init__(**kwargs)
        for line in file_handle:

            # stop reading the file when fasta begins
            # assumes fasta, if present, will always be after all GFF entries
            if line.startswith('>'):
                break
            # ignore comment lines
            elif line.startswith("#"):
                pass

            else:
                gff = GFF(line, file_format=file_format)

                if file_format == 'gff3' and re.match('transcript|mRNA', gff.feature):
                    transc_id = 'ID'
                    gene_id = 'Parent'
                    name = gff.attrib_dict[transc_id]
                    self[name] = {
                        'gene': gff.attrib_dict[gene_id],
                        'strand': gff.strand
                    }

                elif re.match('exon|CDS', gff.feature):
                    if file_format == 'gff3':
                        transc_id = 'Parent'

                    elif file_format == 'gtf':
                        transc_id = 'transcript_id'
                        gene_id = 'gene_id'

                    name = gff.attrib_dict[transc_id]

                    if file_format == 'gtf' and name not in self:
                        self[name] = {
                            'gene': gff.attrib_dict[gene_id],
                            'strand': gff.strand
                        }

                    if name in self:
                        starts = gff.feature + '_starts'
                        ends = gff.feature + '_ends'
                        if starts not in self[name]:
                            self[name][starts] = gff.start
                            self[name][ends] = gff.end

                            if gff.feature == 'CDS':
                                self[name]['frame'] = gff.frame

                        else:
                            self[name][starts] += ',%s' % gff.start
                            self[name][ends] += ',%s' % gff.end
                            if gff.feature == 'CDS':
                                self[name]['frame'] += ',%s' % gff.frame

#
def coords2extb(coords, f_out):
    for trans in coords:
        # strand = coords[trans]['strand']
        # exons, introns = calc_extb_features(coords[trans])
        # print(trans, strand, len(exons), array2str(exons), array2str(introns), sep='\t', file=f_out)
        print(coords[trans], file=f_out)


def attributes_parser(attributes, file_format='gff3'):
    if file_format == 'gff3':
        pattern = re.compile(r'^\s*(\S+)\s*=\s*(.*)\s*$')
    elif file_format == 'gtf':
        pattern = re.compile(r'^\s*(\S+)\s+\"([^\"]+)\"\s*')
    else:
        sys.exit('Unsupported format: %s')

    from html import unescape
    attributes = unescape(attributes)

    attrib_dict = {}
    atts = re.sub(';\s*$', '', attributes)
    atts = atts.split(';')
    for att in atts:
        g = re.search(pattern, att)
        try:
            k, v = g.group(1, 2)
            v = re.sub(r'^(transcript|gene):', '', v)
            attrib_dict[k] = v
        except:
            sys.exit('PARSING ERROR: regex %s failed to parse: "%s"' % (str(pattern), attributes))
    return attrib_dict



if __name__ == '__main__':
    # Deal with errno 32 (broken pipe)
    # http://stackoverflow.com/a/16865106
    from signal import signal, SIGPIPE, SIG_DFL

    signal(SIGPIPE, SIG_DFL)

    main()
