import re
from collections import OrderedDict

import numpy as np

from parse_gff import attributes_parser
from .exceptions import *
from .utils import str2array, array2str, rand_id


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


class GffLine(object):
    def __init__(self, line: str, file_format='gff3'):

        assert type(line) is str, '%s not a string' % line

        self.field = line.strip().split('\t')
        if len(self.field) != 9:
            raise Exception('%s doesnt have 9 fields' % line)

        self.chrom = self.field[0]
        self.source = self.field[1]
        self.feature = self.field[2]
        self.start = self.field[3]
        self.end = self.field[4]
        self.score = self.field[5]
        self.strand = self.field[6]
        self.frame = self.field[7]
        self.attributes = self.field[8]

        self.file_format = file_format

    @property
    def has_multiple_parents(self):
        if self.file_format == 'gff3':
            if 'Parent' in self.attrib_dict:
                parent = self.attrib_dict['Parent'].split(",")
                if len(parent) > 1:
                    return True
        return False

    @property
    def attrib_dict(self):
        return attributes_parser(self.attributes, file_format=self.file_format)

    # ToDo: URGENT => fix this
    @property
    def id(self):

        if self.file_format == 'gff3':
            if re.match('transcript|mRNA', self.feature):
                return self.attrib_dict['ID']

            elif re.match('exon|CDS', self.feature) and self.file_format == 'gff3':
                return self.attrib_dict['Parent']

            else:
                return self.attrib_dict['ID']

        elif self.file_format == 'gtf':
            if re.match('exon|CDS|transcript|mRNA', self.feature):
                return self.attrib_dict['transcript_id']

            else:
                return self.attrib_dict['gene_id']

        else:
            return rand_id()

    @property
    def gene_id(self):
        # gff and gtf differ on how to get this
        if self.file_format == 'gff3' and re.match('transcript|mRNA', self.feature):
                gene_key = 'Parent'
        elif self.file_format == 'gff3' and self.feature == 'gene':
            gene_key = 'ID'

        # expected for gtf, but some GFFs have this
        else:
            gene_key = 'gene_id'

        return self.attrib_dict[gene_key]

    @property
    def transcript_id(self):
        if self.file_format == 'gff3' and re.match('exon|CDS', self.feature):
            rna_key = 'Parent'

        elif self.file_format == 'gff3' and re.match('transcript|mRNA', self.feature):
            rna_key = 'ID'
        else:
            rna_key = 'transcript_id'

        return self.attrib_dict[rna_key]


class GffItem(AttribDict):
    """
        An tem parsed from a GFF/GTF file

        the attributes of this object can be accessed as keys from a dict()

    """
    def __init__(self, gff_line: GffLine, **kwargs):
        super(GffItem, self).__init__(**kwargs)
        # super().__init__(**kwargs)

        if gff_line:
            if type(gff_line) == str:
                gff_line = GffLine(gff_line)

            assert type(gff_line) == GffLine

            self.file_format = gff_line.file_format
            self.chrom = gff_line.chrom
            self.strand = gff_line.strand
            self.source = gff_line.source
            self.attrib = gff_line.attrib_dict

        string_keys = {'gene_id', 'transcript_id', 'chrom',
                   'source', 'strand', }

        coord_keys = {'exon_starts',
                      'exon_ends',
                      'CDS_starts',
                      'CDS_ends',
                      'frame',  # gff3 uses phase, gff2/gtf uses frame
                      }

        # def keys: union of id + coord_keys
        def_keys = string_keys | coord_keys

        # create default keys/attributes
        for k in def_keys:
            if k in kwargs:
                self[k] = kwargs[k]
            elif k not in self:
                if k in string_keys:
                    self[k] = None
                elif k in coord_keys:
                    self[k] = ''
        if not self.attrib:
            self.attrib = {}

        # update attributes passed on init
        self.attrib.update((k, v) for (k, v) in kwargs if k not in def_keys)



class GFF(OrderedDict):
    orientation = 'Unknown'
    specie = 'Unknown'
    file_format = 'Unknown'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def to_exons_file(self, f_out: str):
        import io
        remember_to_close = False
        if not isinstance(f_out, io.TextIOWrapper):
            remember_to_close = True
            f_out = open(f_out, 'w')

        for rna in self:
            gff_item = self[rna]
            annotation = GenomicAnnotation(gff_item.exon_starts, gff_item.exon_ends,
                                           gff_item.strand, orientation=self.orientation)

            chr_str = '%s:%s-%s' % (gff_item.chrom, annotation.starts[0], annotation.ends[-1])

            print(gff_item.id,
                  gff_item.gene_id,
                  chr_str,
                  len(annotation),
                  sum(annotation.exons),
                  gff_item.strand,
                  array2str(annotation.exons),
                  array2str(annotation.introns),
                  array2str(annotation.starts),
                  file=f_out, sep='\t')

        if remember_to_close:
            f_out.close()


class GenomicAnnotation(object):
    cds_starts = cds_ends = np.array([np.nan])

    def __init__(self, starts, ends, strand, cds_starts=None, cds_ends=None,
                 starts_offset=1, orientation='Unknown', **kwargs):

        for k, v in kwargs.items():
            setattr(self, k, v)

        self.starts = str2array(starts)
        self.ends = str2array(ends)

        if not re.match('-|\+', strand):
            raise Exception('invalid strand value: %s' % strand)
        self.strand = strand

        assert len(self.starts) == len(self.ends)

        # make starts 0-based
        self.starts -= starts_offset

        if cds_starts and cds_ends:
            self.cds_starts = str2array(cds_starts)
            self.cds_ends = str2array(cds_ends)

            assert len(self.cds_starts) == len(self.cds_ends)

        self.__fix_orientation__(orientation)

    def __len__(self):
        return len(self.starts)

    @property
    def exons(self):
        exons = self.ends - self.starts
        return exons

    @property
    def introns(self):
        introns = np.array([np.nan])
        if len(self) > 1:
            introns = self.starts[1:] - self.ends[:-1]
        return introns
        # for intron in introns:
        #     yield intron

    def __fix_orientation__(self, orientation='Unknown'):

        if orientation != 'genomic' and self.strand == '-' and len(self) > 1:

            if orientation == 'transcript':
                self.__reverse__()

            elif orientation == 'Unknown':
                diff = self.starts[-1] - self.starts[0]
                if diff < 0:
                    self.__reverse__()

    def __reverse__(self):
        self.starts = self.starts[::-1]
        self.ends = self.ends[::-1]
        self.cds_starts = self.cds_starts[::-1]
        self.cds_ends = self.cds_ends[::-1]


