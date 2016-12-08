import re

from .classes import GFF, GffLine
from .line_parser import guess_kind_of_gff, line_parser
from genial.utils import str2array


def parse_to_dict(file_handle, ff='Unknown'):
    gff_dict = GFF()
    gff_dict.file_format = ff

    for gff_line in line_parser(file_handle, gff_dict.file_format):
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
    # try:
    rna_ids = gff_line.parents_of_exon
    # except KeyError:
    #     pass
    # if len(list(rna_ids)) > 0:
    for rna_id in rna_ids:
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
