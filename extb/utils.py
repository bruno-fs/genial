from sys import intern
import os

import numpy as np
import pandas as pd
from .exceptions import UnsupportedFile


def detect_mime(path_to_file, uncompress=False):
    import magic
    mime = magic.Magic(mime=True, uncompress=uncompress)
    mime_from_file = mime.from_file(path_to_file)
    # magic returns byte, not str
    # ToDo: remove decode when issue #98 is closed
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


def stringfy(obj):
    """
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


def str2array(string):
    return np.fromstring(string, sep=',', dtype=np.int64)
    # import re
    # string = re.sub(r',$', '', string)
    # items = string.split(',')
    # size = len(items)
    # buffer = np.array([int(x) for x in items])
    # arr = Array(shape=(size,), buffer=buffer, dtype=np.int64)
    # return arr

# # noinspection PyClassHasNoInit
# class Array(np.ndarray):
#     """a subclass from numpy.ndarray
#
#     the ONLY difference is its string representation for 1D arrays
#     (with shape = (n, )) will be provided by array2str
#
#     """
#
#     def __str__(self):
#         tup = self.shape
#         if len(tup) == 1:
#             return array2str(self)
#         else:
#             return super(Array, self).__str__()


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
