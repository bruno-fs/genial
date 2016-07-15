#!/usr/bin/env python3

import argparse as argp
import os
import sys

from extb.utils import magic_open
from extb import parse


def main():
    arg_parser = argp.ArgumentParser(description="extract blockSizes from a gff file")
    arg_parser.add_argument('--input', help="input file", )    # default=sys.stdin)
    arg_parser.add_argument('--output', help='output file', )  # default=sys.stdout)
    arg_parser.add_argument('-f', '--input_format', help='gtf or gff3 (default: gff3)')
    arg_parser.add_argument('-t', '--output_format', help='gtf or gff3 (default: gff3)')

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
        if input_format not in ['gff3', 'gtf', 'bed12']:
            sys.exit('ERROR: %s extension is not supported' % input_format)

    if args.output_format:
        output_format = args.output_format
        if output_format not in ['gff3', 'gtf', 'bed12', 'extb']:
            sys.exit('ERROR: %s extension is not supported' % output_format)

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

    for annotation in parse(f_in, input_format)
        print(annotation.format(output_format), file=f_out, sep='\t')

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
