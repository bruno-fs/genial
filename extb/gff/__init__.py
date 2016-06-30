import re

from .classes import GFF, GffLine
from .line_parser import get_format_file
from extb.utils import str2array


def parse(file_handle, ff='Unknown'):
    gff_dict = GFF()
    gff_dict.file_format = ff

    for line in file_handle:

        # stop reading the file when fasta begins
        # assumes fasta, if present, will always be after all GFF entries
        if line.startswith('>'):
            break
        # ignore comment lines
        elif line.startswith("#"):
            pass

        else:
            # detect format
            if gff_dict.file_format == 'Unknown':
                ff = get_format_file(line)
                gff_dict.file_format = ff
                import sys
                print("detected format {}".format(ff), file=sys.stderr)

            gff_line = GffLine(line, file_format=ff)

            if re.match('exon|CDS', gff_line.feature):
                add_exon(gff_dict, gff_line)

            else:
                gff_dict.add_kinship(gff_line)

    return gff_dict


def add_exon(gff_dict: GFF, gff_line: GffLine):
    from sys import intern
    starts = intern(gff_line.feature + '_starts')
    ends = intern(gff_line.feature + '_ends')

    # deal with gff with multiple parents_of_exon
    for rna_id in gff_line.parents_of_exon:
        if rna_id not in gff_dict:
            gff_dict[rna_id] = gff_line
            gff_dict.add_attribs(rna_id, gff_line)

        gff_dict[rna_id][starts] += '%s,' % gff_line.start
        gff_dict[rna_id][ends] += '%s,' % gff_line.end
        if gff_line.feature == 'CDS':
            gff_dict[rna_id].frame += '%s,' % gff_line.frame

        # detect orientation of gff
        # ToDo: create another function to do this
        if gff_dict.orientation == 'Unknown' and gff_line.strand == '-':
            arr = str2array(gff_dict[rna_id][starts])
            if len(arr) > 1:
                dif = arr[-1] - arr[0]

                if dif > 0:
                    gff_dict.orientation = 'genomic'
                elif dif < 0:
                    gff_dict.orientation = 'transcript'