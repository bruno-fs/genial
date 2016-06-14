from extb.exceptions import UnsupportedFile, ParseError
from .attrib_parser import attributes_parser
from .classes import GffLine


def get_format_file(gff_line: str, ffs=['gff3', 'gtf']):
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
        ff = get_format_file(gff_line, ffs)

    return ff