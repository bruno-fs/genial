#!/usr/bin/env python3

import argparse as argp
import gzip
import os
from .classes import *

__author__ = 'bruno'


def main():
    arg_parser = argp.ArgumentParser(description="extract blockSizes from a gff file")
    arg_parser.add_argument('-i', '--input', help='input file (default is stdin)')
    arg_parser.add_argument('-o', '--output', help='output file (default is stdout)')
    arg_parser.add_argument('-f', '--input_format', help='gtf or gff3 (default: gff3)')
    arg_parser.add_argument('--specie', help='specie name (default: parsed from filename)')

    args = arg_parser.parse_args()

    f_in = sys.stdin
    f_out = sys.stdout
    input_format = 'gff3'
    spec_name = 'Unknown'

    if args.input_format:
        input_format = args.input_format
        if input_format not in ['gff3', 'gtf']:
            sys.exit('ERROR: %s extension is not supported' % input_format)

    if args.input:
        if not os.path.exists(args.input):
            sys.exit("ERROR: file %s doesn't exist" % args.input)

        # ToDO: use libmagic to detect file type
        elif re.search(r'\.gz$', args.input):
            f_in = gzip.open(args.input, 'rt')

        else:
            f_in = open(args.input, 'r')
        if not args.specie:
            spec_name = args.input.split('.')[0]

    if args.output:
        if os.path.exists(args.output):
            sys.exit('ERROR: %s already exists!!!' % args.output)
        else:
            f_out = open(args.output, 'w')

    gff_dict = gff_parser(f_in, input_format)
    gff_dict.specie = spec_name

    gff_dict2extb(gff_dict, f_out)

    f_in.close()
    f_out.close()


def gff_parser(file_handle, file_format):
    gff_dict = GFF()

    for line in file_handle:

        # stop reading the file when fasta begins
        # assumes fasta, if present, will always be after all GFF entries
        if line.startswith('>'):
            break
        # ignore comment lines
        elif line.startswith("#"):
            pass

        else:
            gff_line = GffLine(line, file_format=file_format)

            if re.match('transcript|mRNA', gff_line.feature):

                rna_id = gff_line.id

                if rna_id not in gff_dict:
                    gff_dict[rna_id] = GffItem(gff_line)

            elif re.match('exon|CDS', gff_line.feature):

                rna_id = gff_line.id

                if rna_id not in gff_dict:
                    gff_dict[rna_id] = GffItem(gff_line)

                starts = gff_line.feature + '_starts'
                ends = gff_line.feature + '_ends'

                gff_dict[rna_id][starts] += '%s,' % gff_line.start
                gff_dict[rna_id][ends] += '%s,' % gff_line.end
                if gff_line.feature == 'CDS':
                    gff_dict[rna_id].frame += '%s,' % gff_line.frame

                if gff_dict.orientation == 'Unknown' and gff_line.strand == '-':
                    arr = str2array(gff_dict[rna_id][starts])
                    if len(arr) > 1:
                        dif = arr[-1] - arr[0]

                        if dif > 0:
                            gff_dict.orientation = 'genomic'
                        elif dif < 0:
                            gff_dict.orientation = 'transcript'

    return gff_dict


def gff_dict2extb(gff_dic: GFF, f_out):
    for rna in gff_dic:
        gff_item = gff_dic[rna]
        annotation = GenomicAnnotation(gff_item.exon_starts, gff_item.exon_ends,
                                       gff_item.strand, orientation=gff_dic.orientation)

        chr_str = '%s:%s-%s' % (gff_item.chrom, annotation.starts[0], annotation.ends[0])

        print(gff_item.id, gff_item.gene_id, chr_str, annotation.len, gff_item.strand,
              array2str(annotation.exons), array2str(annotation.introns), array2str(annotation.starts),
              file=f_out, sep='\t')


if __name__ == '__main__':
    # Deal with broken pipe error (errno 32)
    # http://stackoverflow.com/a/16865106
    from signal import signal, SIGPIPE, SIG_DFL

    signal(SIGPIPE, SIG_DFL)

    main()
