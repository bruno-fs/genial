import re

from .utils import get_format_file
from .helper_classes import GFF, GffItem
from gff.classes import GffLine, GffItem, GFF
from .utils import str2array


def gff_parser(file_handle, ff='Unknown'):
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
            if ff == 'Unknown':
                ff = get_format_file(line)

            gff_line = GffLine(line, file_format=ff)

            if re.match('transcript|mRNA', gff_line.feature):

                rna_id = gff_line.id

                if rna_id not in gff_dict:
                    gff_dict[rna_id] = GffItem(gff_line)

            elif re.match('exon|CDS', gff_line.feature):
                # if gff_line.has_multiple_parents:
                #     raise MultipleParentsGFF('fields with multiple Parents are currently not supported')

                add_exon_or_cds_to_gff(gff_line, gff_dict)

    return gff_dict


def add_exon_or_cds_to_gff(gff_line: GffLine, gff_dict: GFF):
    starts = gff_line.feature + '_starts'
    ends = gff_line.feature + '_ends'

    # support for gff with multiple parents
    rna_ids = gff_line.id.split(',')
    for rna_id in rna_ids:
        if rna_id not in gff_dict:
            gff_dict[rna_id] = GffItem(gff_line)

        gff_dict[rna_id][starts] += '%s,' % gff_line.start
        gff_dict[rna_id][ends] += '%s,' % gff_line.end
        if gff_line.feature == 'CDS':
            gff_dict[rna_id].frame += '%s,' % gff_line.frame

        # detect orientation of gff
        if gff_dict.orientation == 'Unknown' and gff_line.strand == '-':
            arr = str2array(gff_dict[rna_id][starts])
            if len(arr) > 1:
                dif = arr[-1] - arr[0]

                if dif > 0:
                    gff_dict.orientation = 'genomic'
                elif dif < 0:
                    gff_dict.orientation = 'transcript'