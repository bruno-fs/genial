#!/usr/bin/env python3

import argparse as argp
import os
import sys

try:
    from extb import gff_parser
    from extb.utils import magic_open, str2array

except ImportError:
    # from .helper_classes import gff_parser, gff_dict2extb
    # from .utils import magic_open, str2array

    # get the REAL script dir, even if the script is a link
    *parent_dir, dir_name = os.path.dirname(os.path.realpath(__file__)).split('/')
    parent_dir = '/'.join(parent_dir)
    # add to path
    sys.path.append(parent_dir)

    # import regardless the module name
    import importlib
    extb = importlib.import_module(dir_name)

    gff_dict2extb = getattr(extb, 'gff_dict2extb')
    gff_parser = getattr(extb, 'gff_parser')

    utils = getattr(extb, 'utils')
    magic_open = getattr(utils, 'magic_open')
    str2array = getattr(utils, 'str2array')


def main():
    arg_parser = argp.ArgumentParser(description="extract blockSizes from a gff file")
    arg_parser.add_argument('-i', '--input', help='input file (default is stdin)')
    arg_parser.add_argument('-o', '--output', help='output file (default is stdout)')
    arg_parser.add_argument('-f', '--input_format', help='gtf or gff3 (default: gff3)')
    arg_parser.add_argument('--specie', help='specie name (default: parsed from filename)')

    args = arg_parser.parse_args()

    # default values
    # they will be replaced if their respective args are provided
    f_in = sys.stdin
    f_out = sys.stdout
    # input_format = 'gff3'
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

        if not args.specie:
            species_name = args.input.split('.')[0]

    if args.output:
        if os.path.exists(args.output):
            sys.exit('ERROR: %s already exists!!!' % args.output)
        else:
            f_out = open(args.output, 'w')

    gff_dict = gff_parser(f_in, input_format)
    gff_dict.specie = species_name

    gff_dict.to_exons_file(f_out)

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
