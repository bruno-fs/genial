import re
import sys

import numpy as np


array2str = lambda arr: ','.join(str(x) for x in arr)
str2array = lambda string: np.fromstring(string, sep=',', dtype=int)

randstr = lambda: 'RandStr_' + str(int(np.random.random()*10**10))


class GffLine:
    def __init__(self, line, file_format='gff3'):
        assert type(line) is str, '%s not a string' % line

        self.field = line.strip().split('\t')
        assert len(self.field) == 9, '%s doesnt have 9 fields' % line

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
    def attrib_dict(self):
        return attributes_parser(self.attributes, file_format=self.file_format)

    @property
    def id(self):
        if self.file_format == 'gff3' and re.match('transcript|mRNA', self.feature):
            return self.attrib_dict['ID']

        elif re.match('exon|CDS', self.feature) and self.file_format == 'gff3':
            return self.attrib_dict['Parent']

        elif re.match('exon|CDS', self.feature) and self.file_format == 'gtf':
            return self.attrib_dict['transcript_id']

        else:
            return randstr


class GffItem(dict):
    """
        An tem parsed from a GFF/GTF file

        the attributes of this object can be accessed as keys from a dict()

        """
    def __init__(self, gff_line: GffLine, **kwargs):
        super().__init__(**kwargs)

        if gff_line:
            if type(gff_line) == str:
                gff_line = GffLine(gff_line)

            assert type(gff_line) == GffLine

            file_format = gff_line.file_format

            if file_format == 'gff3' and re.match('transcript|mRNA', gff_line.feature):
                gene_key = 'Parent'

            elif re.match('exon|CDS', gff_line.feature):
                if file_format == 'gtf':
                    gene_key = 'gene_id'

            self.id = gff_line.id
            # if the feature is mRNA, transcript, exon or CDS, the id is the same
            # otherwise, the parse function will ignore the line
            self.transcript_id = self.id

            # workaround for gene_id.
            # ToDo: fix this properly
            try:
                self.gene_id = gff_line.attrib_dict[gene_key]
            except UnboundLocalError:
                self.gene_id = randstr

            self.chrom = gff_line.chrom
            self.strand = gff_line.strand
            self.source = gff_line.source
            self.attrib = gff_line.attrib_dict

        default_keys = {'gene_id',
                        'transcript_id',
                        'chrom',
                        'source',
                        'strand',
                        }

        coord_keys = {'exon_starts',
                      'exon_ends',
                      'CDS_starts',
                      'CDS_ends',
                      'frame',  # gff3 uses phase, gff2/gtf uses frame
                      }

        def_keys = default_keys | coord_keys

        # create default keys/attributes
        for k in def_keys:
            if k in kwargs:
                self[k] = kwargs[k]
            elif k not in self:
                if k in default_keys:
                    self[k] = None
                elif k in coord_keys:
                    self[k] = ''
        if not self.attrib:
            self.attrib = {}

        # update attributes passed on init

        self.attrib.update((k, v) for (k, v) in kwargs if k not in def_keys)

    # attributes are dict keys =D
    # [source](http://goodcode.io/articles/python-dict-object/)

    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError("No such attribute: " + name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError("No such attribute: " + name)


class GFF(dict):
    orientation = 'Unknown'
    specie = 'Unknown'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class GenomicAnnotation:
    cds_starts = cds_ends = np.array([np.nan])

    def __init__(self, starts, ends, strand, cds_starts=None, cds_ends=None,
                 starts_offset=1, orientation='guess', **kwargs):

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

    @property
    def len(self):
        return len(self.starts)

    @property
    def exons(self):
        exons = self.ends - self.starts
        return exons

    @property
    def introns(self):
        introns = np.array([np.nan])
        if self.len > 1:
            introns = self.starts[1:] - self.ends[:-1]
        return introns

    def __fix_orientation__(self, orientation='guess'):
        # valid_orientation = ['genomic', 'transcript', 'guess']

        if orientation != 'genomic' and self.strand == '-' and self.len > 1:

            if orientation == 'transcript':
                self.__reverse__()

            elif orientation == 'guess':
                diff = self.starts[-1] - self.starts[0]
                if diff < 0:
                    self.__reverse__()

    def __reverse__(self):
        self.starts = self.starts[::-1]
        self.ends = self.ends[::-1]
        self.cds_starts = self.cds_starts[::-1]
        self.cds_ends = self.cds_ends[::-1]


def attributes_parser(attributes, file_format='gff3'):
    if file_format == 'gff3':
        pattern = re.compile(r'^\s*(\S+)\s*=\s*(.*)\s*$')
    elif file_format == 'gtf':
        pattern = re.compile(r'^\s*(\S+)\s+\"([^\"]+)\"\s*')
    else:
        raise Exception('Unsupported format: %s')

    from html import unescape
    attributes = unescape(attributes)

    attrib_dict = {}
    atts = re.sub(';\s*$', '', attributes)
    atts = atts.split(';')
    for att in atts:
        g = re.search(pattern, att)
        try:
            k, v = g.group(1, 2)
            v = re.sub(r'^(transcript|gene):', '', v)
            attrib_dict[k] = v
        except AttributeError:
            sys.exit('PARSING ERROR: regex %s failed to parse: "%s"' % (str(pattern), attributes))
    return attrib_dict
