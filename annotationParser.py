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
    ap.add_argument('-mg', '--merge_small_gaps', type=int, default=False)
    ap.add_argument('-igs', '--ignore_gaps_smaller_than', type=int, default=False)
    ap.add_argument('-igb', '--ignore_gaps_bigger_than', type=int, default=False)

    args = ap.parse_args()

    if args.input_format:
        input_format = args.input_format
        if input_format not in input_formats:
            raise SystemExit('ERROR: %s extension is not supported' % input_format)

    if args.output_format:
        output_format = args.output_format
        if output_format not in output_formats:
            raise SystemExit('ERROR: %s extension is not supported' % output_format)

    if args.input:
        if args.input == 'stdin':
            f_in = sys.stdin
        elif not os.path.exists(args.input):
            raise SystemExit("ERROR: input file %s doesn't exist" % args.input)
        else:
            f_in = magic_open(args.input)
    else:
        # ap.print_help()
        ap.print_usage()
        quit()

    if args.output:
        if args.output == sys.stdout:
            f_out = sys.stdout

        elif os.path.exists(args.output):
            raise SystemExit('ERROR: %s already exists!!!' % args.output)

        else:
            # dirname = os.path.dirname(args.output)
            # os.makedirs(dirname, exist_ok=True)
            # print('created dir', dirname, file=sys.stderr)
            f_out = open(args.output, 'w')

    for annotation in parse(f_in, input_format):
        if annotation.blockCount() > 1:
            if args.merge_small_gaps:
                small_gap = args.merge_small_gaps
                annotation = annotation.merge_small_gap(small_gap)

            # ignore small introns
            if args.ignore_gaps_smaller_than:
                small_gap = args.ignore_gaps_smaller_than
                gaps = annotation.introns
                small_gaps_count = np.sum(gaps <= small_gap)
                if small_gaps_count > 0:
                    continue

            # ignore HUGE introns
            if args.ignore_gaps_bigger_than:
                huge_gap = args.ignore_gaps_bigger_than
                gaps = annotation.introns
                huge_gaps_count = np.sum(gaps > huge_gap)
                if huge_gaps_count > 0:
                    continue
        # try:
        print(annotation.format(output_format), file=f_out, sep='\t')
        # except IndexError:
        #     print(annotation)

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
