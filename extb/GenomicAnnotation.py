import re

import numpy as np

from .utils import str2array, array2str


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

        self._fix_orientation(orientation)

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

    @property
    def start(self):
        return self.starts[0]

    @property
    def end(self, ei=False):
        if ei:
            return 'ei'
        return self.ends[-1]

    def __str__(self):
        # super(GenomicAnnotation, self).__str__()
        # return '%s:%s-%s\t%s' % (self.chrom, self.start + 1, self.end, array2str(self.exons))
        return '{} exons: {}'.format(self.transcript_id,
                                     # array2str(
                                         self.exons)#)

    def format(self, format):
        if format == 'extb':
            # return str(self)
            chr_str = '%s:%s-%s' % (self.chrom, self.start + 1, self.end)
            extb = [
                self.transcript_id,
                self.gene_id,
                chr_str,
                len(self),
                sum(self.exons),
                self.strand,
                array2str(self.exons),
                array2str(self.introns),
                array2str(self.cds),
                array2str(self.starts),
                ]
            return '\t'.join(str(x) for x in extb)

        elif format == 'bed':
            bed = [
                self.chrom,
                self.start,
            ]
            return '\t'.join(str(x) for x in bed)
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
        self.starts = self.starts[::-1]
        self.ends = self.ends[::-1]
        self.cds_starts = self.cds_starts[::-1]
        self.cds_ends = self.cds_ends[::-1]