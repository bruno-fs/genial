import re

import numpy as np

from .utils import str2array


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
            cds_starts = str2array(cds_starts)
            cds_ends = str2array(cds_ends)

            assert len(cds_starts) == len(cds_ends)
            # make starts 0-based
            cds_starts -= 1
            self.cds_starts = cds_starts
            self.cds_ends = cds_ends

        # self.__fix_orientation__(orientation)

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


    @property
    def cds(self):
        return self.cds_ends - self.cds_starts



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