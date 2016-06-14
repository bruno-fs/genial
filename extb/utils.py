import numpy as np

from .exceptions import UnsupportedFile, ParseError
from .gff.classes import GffLine
from .gff.attrib_parser import attributes_parser



def detect_mime(path_to_file, uncompress=False):
    import magic
    mime = magic.Magic(mime=True, uncompress=uncompress)
    mime_from_file = mime.from_file(path_to_file)
    # magic returns byte, not str
    # ToDo: remove decode when issue #98 is closed
    mime_from_file = mime_from_file.decode()
    return mime_from_file


def magic_open(path_to_file, mode='rt'):
    import gzip
    mime = detect_mime(path_to_file)

    if mime == 'text/plain':
        return open(path_to_file, mode=mode)

    elif mime == 'application/gzip':
        # if detect_mime(path_to_file, uncompress=True) == 'text/plain':
        return gzip.open(path_to_file, mode=mode)
    else:
        raise UnsupportedFile('File %s is type %s' % (path_to_file, mime))


def array2str(arr):
    return ','.join(str(x) for x in arr)


def str2array(string):
    return np.fromstring(string, sep=',', dtype=np.int64)


def rand_id():
    import string as s
    chars = s.ascii_lowercase + s.digits
    chars = [x for x in chars]
    _id = ''.join(np.random.choice(chars, 15))
    return "Random_ID_" + _id


class AttribDict(dict):
    # attributes are dict keys =D
    # [source](http://goodcode.io/articles/python-dict-object/)

    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("attribute %s doesn't exist" % name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("attribute %s doesn't exist" % name)


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