

from gffdict import attributes_parser
import numpy as np


class GFF:
    def __init__(self, line, file_format='gff3'):
        """

        Parameters
        ----------
        line
        file_format

        Returns
        -------

        """
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


class GenomicAnnotation:
    cds_starts = cds_ends = np.array([np.nan])

    def __init__(self, starts, ends, cds_starts=None, cds_ends=None,
                 starts_offset=1, orientation='guess', **kargs):

        for k, v in kargs.items():
            setattr(self, k, v)

        self.starts = str2array(starts)
        self.ends = str2array(ends)

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
        introns = np.array(np.nan)
        if self.len > 1:
            introns = self.starts[1:] - self.ends[:-1]
        return introns

    def __fix_orientation__(self, orientation):
        # valid_orientation = ['genomic', 'transcript', 'guess']
        # assert orientation in valid_orientation, 'orientation "%s" is invalid (possible values ' \
        #                                          'are %s)' % (orientation, ', '.join(valid_orientation))

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


array2str = lambda arr: ','.join(str(x) for x in arr)
str2array = lambda string: np.fromstring(string, sep=',', dtype=int)