import re
from collections import OrderedDict

import numpy as np

from .utils import str2array, array2str, rand_id


class ParseError(Exception):
    pass


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
    ff = ffs.pop()

    try:
        attributes_parser(attrib, ff)
    except ParseError:
        ff = get_format_file(gff_line, ffs)

    # print('funcionou', c, ff)
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
    def attrib_dict(self):
        return attributes_parser(self.attributes, file_format=self.file_format)

    # ToDo: URGENT => fix this
    @property
    def id(self):
        if self.file_format == 'gff3' and re.match('transcript|mRNA', self.feature):
            return self.attrib_dict['ID']

        elif re.match('exon|CDS', self.feature) and self.file_format == 'gff3':
            return self.attrib_dict['Parent']

        elif re.match('exon|CDS', self.feature) and self.file_format == 'gtf':
            return self.attrib_dict['transcript_id']

        else:
            return rand_id


class GffItem(AttribDict):
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
                self.gene_id = rand_id

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


class GFF(OrderedDict):
    orientation = 'Unknown'
    specie = 'Unknown'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class GenomicAnnotation(object):
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
    def compile_pattern(ff):
        if ff == 'gff3':
            return re.compile(r'^\s*(\S+)\s*=\s*(.*)\s*$')
        elif ff == 'gtf':
            return re.compile(r'^\s*(\S+)\s+\"([^\"]+)\"\s*')
        else:
            raise Exception('Unsupported format: %s' % ff)

    from html import unescape
    attributes = unescape(attributes)

    attrib_dict = {}
    atts = re.sub(';\s*$', '', attributes)
    atts = atts.split(';')
    pattern = compile_pattern(file_format)

    for att in atts:
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



def gff_dict2extb(gff_dic: GFF, f_out):
    for rna in gff_dic:
        gff_item = gff_dic[rna]
        annotation = GenomicAnnotation(gff_item.exon_starts, gff_item.exon_ends,
                                       gff_item.strand, orientation=gff_dic.orientation)

        chr_str = '%s:%s-%s' % (gff_item.chrom, annotation.starts[0], annotation.ends[0])

        print(gff_item.id, gff_item.gene_id, chr_str, annotation.len, gff_item.strand,
              array2str(annotation.exons), array2str(annotation.introns), array2str(annotation.starts),
              file=f_out, sep='\t')


def gff_parser(file_handle, ff: str):
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

                rna_id = gff_line.id

                if rna_id not in gff_dict:
                    gff_dict[rna_id] = GffItem(gff_line)

                starts = gff_line.feature + '_starts'
                ends = gff_line.feature + '_ends'

                gff_dict[rna_id][starts] += '%s,' % gff_line.start
                gff_dict[rna_id][ends] += '%s,' % gff_line.end
                if gff_line.feature == 'CDS':
                    gff_dict[rna_id].frame += '%s,' % gff_line.frame

                if gff_dict.orientation == 'Unknown' and gff_line.strand == '-':
                    arr = str2array(gff_dict[rna_id][starts])
                    if len(arr) > 1:
                        dif = arr[-1] - arr[0]

                        if dif > 0:
                            gff_dict.orientation = 'genomic'
                        elif dif < 0:
                            gff_dict.orientation = 'transcript'

    return gff_dict
