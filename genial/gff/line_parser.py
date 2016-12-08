from genial.exceptions import UnsupportedFile, ParseError
from .attrib_parser import attributes_parser
from .classes import GffLine


def guess_kind_of_gff(gff_line: str, ffs=['gff3', 'gtf']):
    """Detect whether given attribute field
    Parameters
    ----------
    gff_line
    ffs

    Returns
    -------

    """
    if type(gff_line) != GffLine:
        gff_line = GffLine(gff_line)

    attrib = gff_line.attributes
    ffs = ffs.copy()
    try:
        ff = ffs.pop()
    except IndexError:
        raise UnsupportedFile("This file doesn't look like gtf/gff")

    try:
        attributes_parser(attrib, ff)
    except ParseError:
        ff = guess_kind_of_gff(gff_line, ffs)

    return ff


def line_parser(file_handle, ff='Unknown'):
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
            if ff == 'Unknown':
                ff = guess_kind_of_gff(line)
                import sys

                print("detected format {}".format(ff), file=sys.stderr)

            yield GffLine(line, file_format=ff)