import re
from collections import OrderedDict

from extb.GenomicAnnotation import GenomicAnnotation
from extb.gff.attrib_parser import attributes_parser
from extb.utils import AttribDict, array2str


def find_parent(
        child: str,
        parent_of: dict,
        recursive=True):

    """
    Find the parent or great-great-...-great-parent of a child

    Required Parameters
    -------------------
    child: str
        If this is already the greatest-parent, will return itself
        otherwise, raise KeyError

    parent_of: dict
        dictionary with key = child and value = parent, eg:
            parent_of = {}
            parent_of["child"] = "parent"

    Other Parameters
    ----------------
    recursive: bool (default: True)
        if True, look for greatest-parent of a child.

    Returns
    -------
    itself, the Parent or the greatest-parent
    """
    try:
        parent = parent_of[child]
    except KeyError:
        if child in parent_of.values():
            return child
        raise
    if recursive:
        return find_parent(parent, parent_of)
    else:
        return parent


class GffLine(object):
    def __init__(self, line: str, file_format: str = "Unknown"):

        assert type(line) is str, '%s not a string' % line

        self.field = line.strip().split('\t')
        if len(self.field) != 9:
            raise Exception('%s doesnt have 9 fields' % line)

        # this approach is more beautiful, but my IDE can't guess
        # completions for this
        # fieldNames = [
        #     'chrom',
        #     'source',
        #     'feature',
        #     'start',
        #     'end',
        #     'score',
        #     'strand',
        #     'frame',
        #     'attributes',
        # ]
        #
        # for i, name in enumerate(fieldNames):
        #     setattr(self, name, self.field[i])

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
            # if 'Parent' in self.attrib_dict:
            parents = self.attrib_dict['Parent'].split(",")
            if len(parents) > 1:
                return True
        return False

    @property
    def attrib_dict(self):
        return attributes_parser(self.attributes, file_format=self.file_format)

    @property
    def gene_id(self):
        # gff and gtf differ on how to get this
        print('dentro da linha, tentei achar o gene')
        try:
            return self.attrib_dict['gene_id']

        except KeyError:
            if self.file_format == 'gff3' and re.search(r'transcript|mRNA', self.feature):
                gene_key = 'Parent'

            elif self.file_format == 'gff3' and self.feature == 'gene':
                gene_key = 'ID'

            else:
                return None

            return self.attrib_dict[gene_key]



    @property
    def transcript_id(self):
        if self.file_format == 'gff3' and re.match(r'exon|CDS', self.feature):
            rna_key = 'Parent'

        elif self.file_format == 'gff3' and re.match('transcript|mRNA', self.feature):
            rna_key = 'ID'

        # usual on gtf, some GFFs (from ensembl) have the attribute transcript_id
        else:
            rna_key = 'transcript_id'

        return self.attrib_dict[rna_key]


class GffItem(AttribDict):
    """
        An tem parsed from a GFF/GTF file

        the attributes of this object can be accessed as keys from a dict()

    """
    def __init__(self, gff_line: GffLine = None, **kwargs):
        super(GffItem, self).__init__(**kwargs)

        if gff_line:
            if type(gff_line) == str:
                gff_line = GffLine(gff_line)

            assert type(gff_line) == GffLine

            self.file_format = gff_line.file_format
            self.chrom = gff_line.chrom
            self.strand = gff_line.strand
            self.source = gff_line.source
            self.attrib = gff_line.attrib_dict

            # self.gene_id = gff_line.gene_id
            # print(gff_line.gene_id)
            self.transcript_id = gff_line.transcript_id

        string_keys = {
            # 'gene_id',
            'transcript_id',
            'chrom',
            'source',
            'strand',
        }

        coord_keys = {
            'exon_starts',
            'exon_ends',
            'CDS_starts',
            'CDS_ends',
            'frame',  # gff3 uses phase, gff2/gtf uses frame
        }

        dict_attribs = {
            'attrib',
            'parent_of'
        }

        # def keys: union of id + coord_keys
        def_keys = string_keys | coord_keys | dict_attribs

        # create default keys/attributes
        for k in def_keys:
            if k in kwargs:
                self[k] = kwargs[k]
            elif k not in self:
                if k in string_keys:
                    self[k] = None

                elif k in coord_keys:
                    self[k] = ''
                elif k in dict_attribs:
                    if not hasattr(self, k):
                        setattr(self, k, {})

        # update attributes passed on init
        self.attrib.update((k, v) for (k, v) in kwargs.items() if k not in def_keys)

    @property
    def gene_id(self):
        try:
            gene_id = self.attrib['gene_id']
        except KeyError:
            gene_id = None

        if gene_id:
            return gene_id

        elif self.file_format == 'gff3':
            parent_of = self.parent_of
            if parent_of:
                try:
                    return find_parent(self.transcript_id, self.parent_of)
                except KeyError:
                    pass  # will return None
        return None


class GFF(OrderedDict):
    orientation = 'Unknown'
    specie = 'Unknown'
    file_format = 'Unknown'
    parent_of = {}
    attributes_of = {}

    # def add_item(self, gff_line: GffLine):
    #     rna_id = gff_line.transcript_id
    #     if rna_id not in self:
    #         self[rna_id] = GffItem(gff_line)

    # def __setitem__(self, *args, **kwargs):
    #     super(GFF, self).__setitem__(*args, **kwargs)

    def __setitem__(self, key, value, **kwargs):
        if isinstance(value, GffLine):
            value = GffItem(value, parent_of=self.parent_of)
        elif isinstance(value, str):
            GffItem(GffLine(value), parent_of=self.parent_of)
        elif isinstance(value, GffItem):
            value.parent_of = self.parent_of
        super(GFF, self).__setitem__(key, value, **kwargs)

    def __init__(self, *args, **kwargs):
        super(GFF, self).__init__(*args, **kwargs)

    def to_exons_file(self, f_out: str):
        import io
        remember_to_close = False
        if not isinstance(f_out, io.TextIOWrapper):
            remember_to_close = True
            f_out = open(f_out, 'w')

        for rna in self:
            gff_item = self[rna]
            annotation = GenomicAnnotation(
                starts=gff_item.exon_starts, ends=gff_item.exon_ends,
                strand=gff_item.strand, orientation=self.orientation,
                cds_starts=gff_item.CDS_starts, cds_ends=gff_item.CDS_ends,
                chrom=gff_item.chrom, transcript_id=gff_item.transcript_id,
                gene_id=gff_item.gene_id
            )

            # print(annotation.format('extb'), file=f_out, sep='\t')
            print('{}'.format(annotation), file=f_out, sep='\t')

        if remember_to_close:
            f_out.close()
