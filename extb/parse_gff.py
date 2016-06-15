import re

from .gff.line_parser import get_format_file
from .gff.classes import GffLine, GffItem, GFF
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
                import sys
                print(ff, file=sys.stderr)

            gff_line = GffLine(line, file_format=ff)

            if re.match('transcript|mRNA', gff_line.feature):

                rna_id = gff_line.transcript_id
                # gene_id = gff_line.attrib_dict['gene_id']

                if rna_id not in gff_dict:
                    gff_dict[rna_id] = GffItem(gff_line) #, gene_id=gene_id)

            elif re.match('exon|CDS', gff_line.feature):
                # if gff_line.has_multiple_parents:
                #     raise MultipleParentsGFF('fields with multiple Parents are currently not supported')

                add_exon_or_cds_to_gff(gff_line, gff_dict)

    return gff_dict


def add_exon_or_cds_to_gff(gff_line: GffLine, gff_dict: GFF):
    from sys import intern
    starts = intern(gff_line.feature + '_starts')
    ends = intern(gff_line.feature + '_ends')

    # deal with gff with multiple parents
    rna_ids = gff_line.transcript_id.split(',')
    for rna_id in rna_ids:
        if rna_id not in gff_dict:
            """
            Assuming the GFF file is ordered, all transcripts/mRNAs will be declared
            before their exons, so it is safe to skip it;

            On the other hand, for GTF files it is expected to find both gene_id and transcript_id
            on its attibute dict. If it has this information, lets use it.
            otherwise, skip it.

            """
            if gff_dict.file_format == 'gtf':
                if gff_line.gene_id:
                    gff_dict[rna_id] = GffItem(gff_line)

            continue

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