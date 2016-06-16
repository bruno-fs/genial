import re

import numpy as np

from .utils import str2array, array2str, stringfy


class GenomicAnnotation:
    cds_starts = np.array([np.nan])
    cds_ends = np.array([np.nan])
    chrom = None
    transcript_id = None
    gene_id = None

    def __init__(self, starts, ends, strand, cds_starts=None, cds_ends=None,
                 starts_offset=1, orientation='Unknown', **kwargs):
        """
        Parameters
        ----------
        starts
        ends
        strand
        cds_starts
        cds_ends
        starts_offset: 0 or 1 (default = 1)
            GTF/GFF => use 1
            BED => use 0
            extb => use 0
        orientation
        kwargs
        """
        # create all custom defined kwargs ;)
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
            cds_starts -= starts_offset
            self.cds_starts = cds_starts
            self.cds_ends = cds_ends
        else:
            self.cds_starts = np.array([np.nan])
            self.cds_ends = np.array([np.nan])

        self._fix_orientation(orientation)
        # self._reverse()

    def __len__(self):
        return len(self.starts)

    @property
    def exons(self):
        exons = self.ends - self.starts
        return exons

    @property
    def internal_exons(self):
        # first, *internal, last = self.exons
        # return internal
        return self.exons[1:-1]

    @property
    def introns(self):
        # introns = np.array([np.nan])
        # # introns = None
        if len(self) > 1:
            return self.starts[1:] - self.ends[:-1]
        return 'NA'

    @property
    def cds(self):
        if np.isnan(np.sum(self.cds_starts)):  # and self.cds_ends:
            return 'NA'
        else:
            return self.cds_ends - self.cds_starts

    @property
    def orf_blocks(self):
        if np.isnan(np.sum(self.cds_starts)):
            return 'NA'
        else:
            return self.cds_ends - self.cds_starts

    @property
    def orf_size(self):
        try:
            return sum(self.orf_blocks)
        except TypeError:
            return 'NA'

    @property
    def start(self):
        return self.starts[0]

    @property
    def end(self):
        return self.ends[-1]

    @end.setter
    def end(self, value):
        self.ends[-1] = value

    def __str__(self):
        return '{} -> {}'.format(self.transcript_id, array2str(self.exons))

    def format(self, format):

        if format == 'extb':
            chr_str = '%s:%s-%s' % (self.chrom, self.start + 1, self.end)
            extb = [
                self.transcript_id,
                # self.internal_exons,
                self.gene_id,
                chr_str,
                len(self),
                sum(self.exons),
                self.orf_size,
                self.strand,
                self.exons,
                self.introns,
                # self.cds,
                # self.starts,
                ]

            return '\t'.join(stringfy(x) for x in extb)

        elif format == 'bed':
            bed = [
                self.chrom,
                self.start,
                ]

            return '\t'.join(stringfy(x) for x in bed)
        else:
            super(GenomicAnnotation, self).__format__(format)

    def _fix_orientation(self, orientation='Unknown'):
        if orientation != 'genomic' and self.strand == '-' and len(self) > 1:

            if orientation == 'transcript':
                self._reverse()

            elif orientation == 'Unknown':
                diff = self.starts[-1] - self.starts[0]
                if diff < 0:
                    self._reverse()

    def _reverse(self):
        # print('reversing')
        self.starts = self.starts[::-1]
        self.ends = self.ends[::-1]

        if not np.isnan(np.sum(self.cds_starts)):
            self.cds_starts = self.cds_starts[::-1]
            self.cds_ends = self.cds_ends[::-1]

