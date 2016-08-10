#!/usr/bin/env python3

import argparse as argp
import os
import sys

import numpy as np

from extb.utils import magic_open
from extb import parse, input_formats, output_formats


def main():
    ap = argp.ArgumentParser(description="parse annotation files")
    ap.add_argument('-i', '--input',
                    help="input file. to read from pipe, use the argument 'stdin'")
    ap.add_argument('-o', '--output',
                    help='output file', default=sys.stdout)
    ap.add_argument('-f', '--input_format',
                    help='input file format (supported formats: %s)' % ', '.join(input_formats),
                    default='bed')
    ap.add_argument('-t', '--output_format',
                    help='output file format (supported formats: %s)' % ', '.join(output_formats),
                    default='bed')
    ap.add_argument('-n', '--min_exon_count', type=int, default=1,
                    help='min number of exons (inclusive)')

    ap.add_argument('-igs', '--ignore_gaps_smaller_than', type=int, default=False)
    ap.add_argument('-igb', '--ignore_gaps_bigger_than', type=int, default=False)

    ap.add_argument('-v', '--invert_match', default=False, action='store_true',
                    help='only return false results (similar to grep -v)',)

    argv = ap.parse_args()

    if argv.input_format:
        input_format = argv.input_format
        if input_format not in input_formats:
            raise SystemExit('ERROR: %s extension is not supported' % input_format)

    if argv.output_format:
        output_format = argv.output_format
        if output_format not in output_formats:
            raise SystemExit('ERROR: %s extension is not supported' % output_format)

    if argv.input:
        if argv.input == 'stdin':
            f_in = sys.stdin
        elif not os.path.exists(argv.input):
            raise SystemExit("ERROR: input file %s doesn't exist" % argv.input)
        else:
            f_in = magic_open(argv.input)
    else:
        # ap.print_help()
        ap.print_usage()
        quit()

    if argv.output:
        if argv.output == sys.stdout:
            f_out = sys.stdout

        elif os.path.exists(argv.output):
            raise SystemExit('ERROR: %s already exists!!!' % argv.output)

        else:
            # dirname = os.path.dirname(argv.output)
            # os.makedirs(dirname, exist_ok=True)
            # print('created dir', dirname, file=sys.stderr)
            f_out = open(argv.output, 'w')

    def filter(annotation, min_exon_count, ignore_small_gaps, ignore_big_gaps):
        c = 0

        if annotation.blockCount() < min_exon_count:
            # if we already know the count is above the minimum, just ignore it
            return False

        elif annotation.blockCount() > 1:
            # ignore small introns
            if ignore_small_gaps:
                small_gap = ignore_small_gaps
                gaps = annotation.introns
                small_gaps_count = np.sum(gaps < small_gap)
                if small_gaps_count > 0:
                    c += 1

            # ignore HUGE introns
            if ignore_big_gaps:
                huge_gap = ignore_big_gaps
                gaps = annotation.introns
                huge_gaps_count = np.sum(gaps > huge_gap)
                if huge_gaps_count > 0:
                    c += 1

        if c > 0:
            return False
        else:
            return True

    # list of args to be used on filter function
    args = [argv.min_exon_count,
            argv.ignore_gaps_smaller_than,
            argv.ignore_gaps_bigger_than]

    # argv -v / --invert_match
    match = True
    if argv.invert_match:
        match = False


    for annotation in parse(f_in, input_format):

        if filter(annotation, *args) == match:
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
