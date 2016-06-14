import re

from exceptions import ParseError
from helper_classes import GFF, get_format_file, GffLine, GffItem
from utils import str2array


def attributes_parser(attributes, file_format='gff3'):
    def compile_pattern(ff):
        if ff == 'gff3':
            return re.compile(r'^\s*(\S+)\s*=\s*(.*)\s*$')
        elif ff == 'gtf':
            return re.compile(r'^\s*(\S+)\s+\"([^\"]+)\"\s*')
        else:
            raise Exception('Unsupported format: %s' % ff)

    from html import unescape
    # unescape HTML characters like &amp;
    attributes = unescape(attributes)
    # shamelessly copied from https://github.com/hammerlab/gtfparse/blob/master/gtfparse/line_parsing.py
    # Catch mistaken semicolons by replacing "xyz;" with "xyz"
    # Required to do this since the Ensembl GTF for Ensembl release 78 has
    # gene_name = "PRAMEF6;"
    # transcript_name = "PRAMEF6;-201"
    # attributes.replace(';\"', '\"').replace(";-", "-")

    attrib_dict = {}
    attribs = re.sub(';\s*$', '', attributes)
    attribs = attribs.split(';')
    pattern = compile_pattern(file_format)

    for att in attribs:
        g = re.search(pattern, att)

        try:
            k, v = g.group(1, 2)
            v = re.sub(r'^(transcript|gene):', '', v)
            attrib_dict[k] = v
        except AttributeError:
            # TODO: use/make more specific exceptions
            raise ParseError('regex %s failed to parse %s' % (str(pattern), attributes))
            # sys.exit('PARSING ERROR: regex %s failed to parse: \n%s' % (str(pattern), attributes))

    return attrib_dict


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