import re
from collections import OrderedDict

from extb.GenomicAnnotation import GenomicAnnotation
from extb.gff.attrib_parser import attributes_parser
from extb.utils import rand_id, AttribDict, array2str


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
    # replaced id for transcript_id and gene_id
    # @property
    # def id(self):
    #
    #     if self.file_format == 'gff3':
    #         if re.match('transcript|mRNA', self.feature):
    #             return self.attrib_dict['ID']
    #
    #         elif re.match('exon|CDS', self.feature) and self.file_format == 'gff3':
    #             return self.attrib_dict['Parent']
    #
    #         else:
    #             return self.attrib_dict['ID']
    #
    #     elif self.file_format == 'gtf':
    #         if re.match('exon|CDS|transcript|mRNA', self.feature):
    #             return self.attrib_dict['transcript_id']
    #
    #         else:
    #             return self.attrib_dict['gene_id']
    #
    #     else:
    #         return rand_id()

    @property
    def gene_id(self):
        # # gff and gtf differ on how to get this
        # if self.file_format == 'gff3' and re.match('transcript|mRNA', self.feature):
        #         gene_key = 'Parent'
        # elif self.file_format == 'gff3' and self.feature == 'gene':
        #     gene_key = 'ID'
        #
        # # usual on gtf, some GFFs (from ensembl) have the attribute transcript_id
        # else:
        #     gene_key = 'gene_id'

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

            self.gene_id = gff_line.gene_id
            self.transcript_id = gff_line.transcript_id

        string_keys = {
            'gene_id',
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
        self.attrib.update((k, v) for (k, v) in kwargs.items() if k not in def_keys)


class GFF(OrderedDict):
    orientation = 'Unknown'
    specie = 'Unknown'
    file_format = 'Unknown'

    def __init__(self, *args, **kwargs):
        # super().__init__(*args, **kwargs)
        super(GFF, self).__init__(*args, **kwargs)

    def to_exons_file(self, f_out: str):
        import io
        remember_to_close = False
        if not isinstance(f_out, io.TextIOWrapper):
            remember_to_close = True
            f_out = open(f_out, 'w')

        for rna in self:
            gff_item = self[rna]
            annotation = GenomicAnnotation(gff_item.exon_starts, gff_item.exon_ends,
                                           gff_item.strand, orientation=self.orientation,
                                           cds_starts=gff_item.CDS_starts, cds_ends=gff_item.CDS_ends)

            chr_str = '%s:%s-%s' % (gff_item.chrom, annotation.starts[0], annotation.ends[-1])

            print(gff_item.transcript_id,
                  gff_item.gene_id,
                  chr_str,
                  len(annotation),
                  sum(annotation.exons),
                  gff_item.strand,
                  array2str(annotation.exons),
                  # array2str(annotation.introns),
                  array2str(annotation.cds),
                  # array2str(annotation.starts),
                  file=f_out, sep='\t')

        if remember_to_close:
            f_out.close()
