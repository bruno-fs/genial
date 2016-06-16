import re

from .gff.line_parser import get_format_file
from .gff.classes import GffLine, GFF, GffItem
from .utils import str2array


def gff_parser(file_handle, ff='Unknown'):
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

            if re.match('transcript|mRNA', gff_line.feature):

                rna_id = gff_line.transcript_id

                if rna_id not in gff_dict:
                    gff_dict[rna_id] = GffItem(gff_line)
                    # gff_dict[rna_id] = gff_line

            elif re.match('exon|CDS', gff_line.feature):
                # print('achei exon')
                add_exon_or_cds_to_gff(gff_line, gff_dict)

    return gff_dict


def add_exon_or_cds_to_gff(gff_line: GffLine, gff_dict: GFF):
    from sys import intern
    starts = intern(gff_line.feature + '_starts')
    ends = intern(gff_line.feature + '_ends')
    # print('dentro da func add')
    # deal with gff with multiple parents
    # ToDo: find a better way to support multi-parent exon/CDS features
    rna_ids = gff_line.transcript_id.split(',')
    for rna_id in rna_ids:
        # print('tentando add rna')
        if rna_id not in gff_dict:
            """
            Assuming the GFF file is ordered, all transcripts/mRNAs will be declared
            before their exons, so it is safe to skip it;

            On the other hand, some GTF files ONLY have exon and CDS entries, but have
            gene_id and transcript_id.

            Let's check If this information is available and use it.
            Otherwise, just skip.

            """
            # print('vamos checar se é gtf?')
            if gff_dict.file_format == 'gtf':
                # print('é gtf')
                if gff_line.gene_id:
                    # print('achei gene')
                    gff_dict[rna_id] = gff_line
                else:
                    continue  # skip gtf w/o gene
            else:
                # continue  # skip item for gff not included
                gff_dict[rna_id] = gff_line
        # print('vou add esse rna')

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
