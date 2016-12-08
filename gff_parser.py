#!/usr/bin/env python3

import argparse as argp
import os
import sys

from GenIAL.gff import parse_to_dict
from GenIAL.utils import magic_open
from GenIAL.GenomicAnnotation import GeneAnnotation


def save_to_file(gff_dict, f_out: str, out_format='bed'):
    import io
    remember_to_close = False
    if not isinstance(f_out, io.TextIOWrapper):
        remember_to_close = True
        f_out = open(f_out, 'w')

    for rna in gff_dict:
        gff_item = gff_dict[rna]
        annotation = GeneAnnotation(
            starts=gff_item.exon_starts, ends=gff_item.exon_ends,
            strand=gff_item.strand, orientation=gff_dict.orientation,
            cds_starts=gff_item.CDS_starts, cds_ends=gff_item.CDS_ends,
            chrom=gff_item.chrom, transcript_id=gff_item.transcript_id,
            gene_id=gff_item.gene_id
        )

        # print(annotation.format(out_format), file=f_out, sep='\t')
        print(annotation.format(out_format), file=f_out, sep='\t')

    if remember_to_close:
        f_out.close()

def main():
    arg_parser = argp.ArgumentParser(description="extract blockSizes from a gff file")
    arg_parser.add_argument('input', nargs='?', help="input file", )    # default=sys.stdin)
    arg_parser.add_argument('output', nargs='?', help='output file', )  # default=sys.stdout)
    arg_parser.add_argument('-f', '--input_format', help='gtf or gff3 (default: gff3)')
    arg_parser.add_argument('--specie', help='specie name (default: parse_to_dict the first cha'
                                             'racters before a dot on input_file)')
    arg_parser.add_argument('-l', '--lazy_output_naming', action='store_true',
                            help='name output file based on specie name\nignored when'
                                 'output name is provided')
    arg_parser.add_argument('-s', '--shorten_spec_name', action='store_true',
                            help='ignored when used with --specie')

    args = arg_parser.parse_args()

    # default values
    # they will be replaced if their respective args are provided
    f_in = sys.stdin
    f_out = sys.stdout
    input_format = 'Unknown'
    species_name = 'Unknown'

    if args.input_format:
        input_format = args.input_format
        if input_format not in ['gff3', 'gtf']:
            sys.exit('ERROR: %s extension is not supported' % input_format)

    if args.input:
        if not os.path.exists(args.input):
            sys.exit("ERROR: file %s doesn't exist" % args.input)

        else:
            f_in = magic_open(args.input)

    if args.specie:
        species_name = args.specie
    elif args.input:
        species_name = args.input.split('/')[-1]
        species_name = species_name.split('.')[0]
        if args.shorten_spec_name:
            a, *b = species_name.split('_')
            species_name = a[0] + b[-1]

    if args.output:
        if args.output == sys.stdout:
            pass
        elif os.path.exists(args.output):
            sys.exit('ERROR: %s already exists!!!' % args.output)
        else:
            dirname = os.path.dirname(args.output)
            os.makedirs(dirname, exist_ok=True)
            print(dirname, file=sys.stderr)
            f_out = open(args.output, 'w')

    elif args.lazy_output_naming:
        if os.path.exists(species_name + '.extb'):
            sys.exit('ERROR: %s.extb already exists!!!' % species_name)
        else:
            f_out = open(species_name + '.extb', 'w')

    gff_dict = parse_to_dict(f_in)   # , input_format)
    gff_dict.specie = species_name

    save_to_file(gff_dict, f_out)

    f_in.close()
    f_out.close()

if __name__ == '__main__':
    # Deal with broken pipe error (errno 32)
    # http://stackoverflow.com/a/16865106
    from signal import signal, SIGPIPE, SIG_DFL, SIGINT

    # deal with broken pipe
    signal(SIGPIPE, SIG_DFL)
    # deal with KeyboardInterrupt
    signal(SIGINT, SIG_DFL)

    main()
