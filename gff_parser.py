#!/usr/bin/env python3

import argparse as argp
import os
import sys

try:
    from extb import gff_parser
    from extb.utils import magic_open, str2array

except ImportError:
    # from .helper_classes import gff_parser
    # from .utils import magic_open, str2array

    # get the REAL script dir, even if the script is a link
    *parent_dir, dir_name = os.path.dirname(os.path.realpath(__file__)).split('/')
    parent_dir = '/'.join(parent_dir)
    # add to path
    sys.path.append(parent_dir)

    # import regardless the module name
    import importlib
    extb = importlib.import_module(dir_name)

    gff_parser = getattr(extb, 'gff_parser')

    utils = getattr(extb, 'utils')
    magic_open = getattr(utils, 'magic_open')
    str2array = getattr(utils, 'str2array')


def main():
    arg_parser = argp.ArgumentParser(description="extract blockSizes from a gff file")
    arg_parser.add_argument('input', nargs='?', help="input file", ) #default=sys.stdin)
    arg_parser.add_argument('output', nargs='?', help='output file', ) #default=sys.stdout)
    # arg_parser.add_argument('-i', '--input', help='input file (default is stdin)')
    # arg_parser.add_argument('-o', '--output', help='output file (default is stdout)')
    arg_parser.add_argument('-f', '--input_format', help='gtf or gff3 (default: gff3)')
    arg_parser.add_argument('--specie', help='specie name (default: parse the first cha'
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
        species_name = args.input.split('.')[0]
        if args.shorten_spec_name:
            a, *b = species_name.split('_')
            species_name = a[0] + '_' + b[-1]

    if args.output:
        if args.output == sys.stdout:
            pass
        elif os.path.exists(args.output):
            sys.exit('ERROR: %s already exists!!!' % args.output)
        else:
            f_out = open(args.output, 'w')

    elif args.lazy_output_naming:
        if os.path.exists(species_name + '.extb'):
            sys.exit('ERROR: %s.extb already exists!!!' % species_name)
        else:
            f_out = open(species_name + '.extb', 'w')

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
