import re
from sys import intern
from collections import OrderedDict

from genial.gff.attrib_parser import attributes_parser
from genial.utils import AttribDict, InternDict


def find_parent(
        child: str,
        parent_of: dict,
        recursive=False):

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
    recursive: bool (default: False)
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
    field = ['']*9

    def __str__(self):
        return '\t'.join(self.field)

    def __init__(self, line: str, file_format: str = "Unknown"):

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
    def parents_of_exon(self):
        """for GFF, will return a generator for each parent of an item
            for GTF, return the transcript_id
        """
        if self.file_format == 'gff3':
            if 'Parent' in self.attrib_dict:
            # try:
                parents = self.attrib_dict['Parent'].split(",")
            # except KeyError:
            #     raise KeyError
            # else:
                for parent in parents:
                    yield parent
        else:
            yield self.transcript_id

    @property
    def attrib_dict(self):
        return attributes_parser(self.attributes, file_format=self.file_format)

    @property
    def gene_id(self):
        try:
            return self.attrib_dict['gene_id']

        except KeyError:
            return None



    @property
    def transcript_id(self):
        """
        For genes, return None. Otherwise, try to return a transcript_id.

        Raises KeyError if not found.

        For GFF format, is prefered to access ID and Parent features

        """
        if self.feature == 'gene':
            # return
            raise KeyError('genes don\'t have transcript_id')

        if self.file_format == 'gff3' and re.match(r'exon|CDS', self.feature):
            rna_key = intern('Parent')

        elif self.file_format == 'gff3' and re.match('transcript|mRNA', self.feature):
            rna_key = intern('ID')

        # usual on gtf, some GFFs (from ensembl) have the attribute transcript_id
        else:
            rna_key = intern('transcript_id')

        return self.attrib_dict[rna_key]


class GffItem(AttribDict):
    """
        An tem parsed from a GFF/GTF file

        the attributes of this object can be accessed as keys from a dict()

    """

    def __init__(self, gff_line: GffLine = None, **kwargs):
        super(GffItem, self).__init__(**kwargs)

        if gff_line:
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
            'parent_of',
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

        # update other attributes passed on init
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
            if self.parent_of:
                try:
                    return find_parent(self.transcript_id, self.parent_of)
                except KeyError:
                    pass  # will return None
        return None


class GFF(OrderedDict):
    orientation = 'Unknown'
    specie = 'Unknown'
    file_format = 'Unknown'

    def __setitem__(self, key, value, **kwargs):
        item = self._prepare_item(value)
        super(GFF, self).__setitem__(key, item, **kwargs)

    def _prepare_item(self, value):
        if isinstance(value, GffLine):
            value = GffItem(value, parent_of=self.parent_of)
        elif isinstance(value, str):
            GffItem(GffLine(value), parent_of=self.parent_of)
        elif isinstance(value, GffItem):
            value.parent_of = self.parent_of
        # value = GffItem(value, parent_of=self.parent_of)
        return value

    def __init__(self, *args, **kwargs):
        super(GFF, self).__init__(*args, **kwargs)
        self.parent_of = InternDict()
        self.attributes_of = InternDict()

    def add_attribs(self, key, item: GffLine):
        try:
            self.attributes_of[key].update(item.attrib_dict.copy())
        except KeyError:
            self.attributes_of[key] = item.attrib_dict

    def add_kinship(self, item: GffLine):
        if self.file_format == 'gtf':
            child_key = 'transcript_id'
            parent_key = 'gene_id'
        else:
            child_key = 'ID'
            parent_key = 'Parent'

        try:
            key = item.attrib_dict[child_key]
        except KeyError:
            # print(item) # ToDo: use a logger
            pass
        else:
            try:
                parent = item.attrib_dict[parent_key]
            except KeyError:
                pass  # so its a gene. skip
            else:  # gene wont reach here
                # deal with multiple parents?
                self.parent_of[key] = parent

            self.add_attribs(key, item)


