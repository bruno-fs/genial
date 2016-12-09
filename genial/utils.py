from sys import intern
import os

import numpy as np
import pandas as pd
import re
from .exceptions import UnsupportedFile


def sort_intervals(starts, ends):
    intervals = [(s, e) for s, e in zip(starts, ends)]
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    starts, ends = zip(*sorted_by_lower_bound)
    starts = np.array(starts, dtype=np.int64)
    ends = np.array(ends, dtype=np.int64)
    return starts, ends


def nice_sort(l):
    """ Sort given iterable in the way that humans expect.
    src: http://stackoverflow.com/a/2669120
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def detect_mime(path_to_file, uncompress=False):
    import magic
    mime = magic.Magic(mime=True, uncompress=uncompress)
    mime_from_file = mime.from_file(path_to_file)
    # magic returns byte, not str
    # mime_from_file = mime_from_file.decode()
    return mime_from_file


def magic_open(path_to_file, mode='rt'):
    import gzip
    # follow symlinks
    path_to_file = os.path.realpath(path_to_file)
    mime = detect_mime(path_to_file)

    if mime == 'text/plain':
        return open(path_to_file, mode=mode)

    elif mime == 'application/gzip':
        # if detect_mime(path_to_file, uncompress=True) == 'text/plain':
        return gzip.open(path_to_file, mode=mode)

    raise UnsupportedFile('File %s is type %s' % (path_to_file, mime))


def rand_id():
    import string as s
    chars = s.ascii_lowercase + s.digits
    chars = [x for x in chars]
    _id = ''.join(np.random.choice(chars, 15))
    return "Random_ID_" + _id


def array2str(arr):
    return ','.join(str(x) for x in arr)


def str2array(string):
    return np.fromstring(string, sep=',', dtype=np.int64)


def format_intervals(iterable_with_numbers):
    n_list = sorted(iterable_with_numbers)
    dist = 0
    curr = n_list[0]
    intervals = []

    for i, n in enumerate(n_list):
        if n - curr != dist:
            next_n = n_list[i - 1]

            if next_n - curr == 0:
                intervals.append(str(curr))

            elif next_n - curr == 1:
                intervals.append(str(curr))
                intervals.append(str(next_n))

            else:
                intervals.append('%d-%d' % (curr, next_n))

            curr = n
            dist = 1

        else:
            dist += 1

    if n - curr == 0:
        intervals.append(str(curr))

    elif n - curr == 1:
        intervals.append(str(curr))
        intervals.append(str(n))

    else:
        intervals.append('%d-%d' % (curr, n))

    return intervals


def stringfy(obj):
    """
    Transform an numpy array in a string (same as array2str) and any other kind of obj in a string.

    I was subclassing numpy.ndarray ONLY to override its print method.
    Even though it works, I'd have to ALWAYS import it and always enforce its use.

    All because I was to lazy to type str(obj) for some objects and array2str(obj)
    for others.

    A Better solution was to create a wrapper to do the job for me:
    stringfy!
    """

    if isinstance(obj, np.ndarray):
        return array2str(obj)
    return str(obj)


def read_extb(filepath):
    """
    returns a pandas dataframe with the content of the specified extb file.
    transcript_ID is used as index
    
    EXTB (EXon TaBle) columns:
    
      'organism'         # 0
      'assembly_version' # 1
      'coords'           # 2
      'strand'           # 3
      'biotype'          # 4
      'gene_name'        # 5
      'gene_id'          # 6
      'transcript_id'    # 7
      'protein_id'       # 8
      'exons'            # 9
      'introns'          #10
      'cds'              #11
      'phase'            #12
      'coords_array'     #13
      'ORF_coords'       #14
      
    """
    columns = ['organism',        # 0
              'assembly_version', # 1
              'coords',           # 2
              'strand',           # 3
              'biotype',          # 4
              'gene_name',        # 5
              'gene_id',          # 6
              'transcript_id',    # 7
              'protein_id',       # 8
              'exons',            # 9
              'introns',          #10
              'cds',              #11
              'phase',            #12
              'coords_array',     #13
              'ORF_coords',       #14
              ]
    extb = pd.read_csv(filepath, sep='\t', header=None)
    extb.columns = columns
    extb.index = extb.transcript_id
    
    return extb


def read_bed(filepath, cols=12):
    """
    returns a pandas dataframe with the content of the specified BED12 file. 
    transcript_ID is used as index
    
    Parameters
    ----------
    filepath: path to file
    cols: # columns in the bedfile    

    Columns
    -------
    
    chrom 
    chromStart 
    chromEnd 
    name 
    score 
    strand 
    thickStart 
    thickEnd 
    itemRgb 
    blockCount 
    blockSizes --> exons
    blockStarts 

    """
      
    columns = """
    chrom 
    chromStart 
    chromEnd 
    name 
    score 
    strand 
    thickStart 
    thickEnd 
    itemRgb 
    blockCount 
    exons 
    blockStart
    """.split()
    
    bed = pd.read_csv(filepath, sep='\t', header=None, names=columns[:cols])
    bed.index = bed.name
    return bed


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


class InternDict(dict):
    def __setitem__(self, key, value):
        if key is not None:
            key = intern(key)
        if isinstance(value, str):
            value = intern(value)

        super(InternDict, self).__setitem__(key, value)
