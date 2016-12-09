import re
import numpy as np

from .utils import str2array, stringfy, sort_intervals


def _bed6_to_GeneAnnot(bed6):
    # print(bed6)
    lines = bed6.splitlines()
    g_annot = 0
    for line in lines:
        chrom, start, end, name, score, strand = line.split('\t')
        if g_annot == 0:
            g_annot = InteractiveAnnotation(start, end, strand,
                                            chrom=chrom, transcript_id=name, starts_offset=0)
        else:
            g_annot.starts = np.hstack([g_annot.starts, int(start)])
            g_annot.ends = np.hstack([g_annot.ends, int(end)])

    return g_annot


class InteractiveAnnotation:
    def __init__(self, starts, ends, strand, cds_starts=None, cds_ends=None,
                 starts_offset=1, orientation='Unknown', **kwargs):
        """

        A interactive and flexible genomic anotation.

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

        # self.starts = np.array([np.nan])
        # self.ends = np.array([np.nan])
        # self.strand = None
        # self.cds_starts = np.array([np.nan])
        # self.cds_ends = np.array([np.nan])
        self.chrom = None
        self.transcript_id = None
        self.gene_id = None
        self.thickStart = None
        self.thickEnd = None

        # create/update all custom defined kwargs
        for k, v in kwargs.items():
            setattr(self, k, v)

        if type(starts) == str:
            starts = str2array(starts)
            ends = str2array(ends)

        self.starts = starts.copy()
        self.ends = ends.copy()

        if not re.match('-|\+', strand):
            raise Exception('invalid strand value: %s' % strand)
        self.strand = strand

        # assert len(self.starts) == len(self.ends)

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

            self.thickStart = np.min(cds_starts)
            self.thickEnd = np.max(cds_ends)

        else:
            self.cds_starts = np.array([np.nan])
            self.cds_ends = np.array([np.nan])

        self._fix_orientation(orientation)
        # self._reverse()

    def __len__(self):
        return len(self.starts)

    def blockCount(self):
        return len(self.starts)

    def blockSizes(self):
        return self.ends - self.starts

    @property
    def exons(self):
        exons = self.ends - self.starts
        return exons

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
        return np.array([np.nan])

    @property
    def cds(self):
        if np.isnan(np.sum(self.cds_starts)):  # and self.cds_ends:
            return np.array([np.nan])
        else:
            return self.cds_ends - self.cds_starts

    @property
    def orf_blocks(self):
        if np.isnan(np.sum(self.cds_starts)):
            return np.array([np.nan])
        else:
            return self.cds_ends - self.cds_starts

    @property
    def exon_contrib_to_orf(self):
        # if self.cds == 'NA':
        if np.isnan(np.sum(self.cds_starts)):
            return np.zeros_like(self.exons)

        else:
            start = self.cds_starts[0]
            stop = self.cds_ends[-1]
            contrib = np.zeros_like(self.exons, dtype=np.float32)

            if len(self) == 1:
                coding = stop - start
                total = self.ends[-1] - self.starts[0]
                return np.array([coding/total])

            else:
                # find exon with start codon
                for i, g_start in enumerate(self.starts):
                    if start >= g_start:
                        startIndex = i
                        break

                # find exon with stop codon
                for i, g_end in reversed(list(enumerate(self.ends))):
                    if stop > g_end:
                        stopIndex = i + 1
                        break
                    elif stop == g_end:
                        stopIndex = i
                        break

            for i in range(self.blockCount()):
                if i < startIndex:
                    contrib[i] = 0

                elif i == startIndex:
                    contrib[i] = (self.ends[i] - start) / self.orf_size

                elif i < stopIndex:
                    contrib[i] = self.exons[i] / self.orf_size


                elif i == stopIndex:
                    # if stop == self.ends[i]:
                    #     contrib[i] = 1
                    # else:
                        # contrib[i] = (stop - self.starts[i]) / (self.ends[i] - self.starts[i])
                    contrib[i] = (stop - self.starts[i]) / self.orf_size


                else:
                    contrib[i] = 0

            return contrib

    def _find_orf_index(self):
        start = self.thickStart
        stop = self.thickEnd
        for i, g_start in enumerate(self.starts):
            if start >= g_start:
                startIndex = i
                break

        # find exon with stop codon
        for i, g_end in reversed(list(enumerate(self.ends))):
            if stop > g_end:
                stopIndex = i + 1
                break
            elif stop == g_end:
                stopIndex = i
                break

        return startIndex, stopIndex

    @property
    def orf_size(self):
        try:
            return sum(self.orf_blocks)
        except TypeError:
            return np.array([np.nan])

    @property
    def start(self):
        return self.starts[0]

    @property
    def end(self):
        return self.ends[-1]

    # @end.setter
    # def end(self, value):
    #     self.ends[-1] = value

    # def BedTool(self, format='bed'):
    #     try:
    #         from pybedtools import BedTool
    #     except ImportError:
    #         raise ImportError("pybedtools is required for this function")
    #
    #     return BedTool(self.format(format), from_string=True)

    def merge_small_gaps(self, gap=15, pybedtools=False):
        """
        Merge gaps smaller or equals the specified amount.
        Useful to remove gaps that should not be treated as introns.

        Reference: http://codereview.stackexchange.com/a/69249

        Parameters
        ----------
        gap:    int
            size of the gap
        pybedtools: Boolean
            use pybedtools (and bedtools). should be removed in future releases.

        Returns
        -------

        a new GeneAnnotation with gaps properly merged

        """

        if pybedtools:
            # if it has at least one small gap, call bedtools merge
            bed = self.format('bed6')
            try:
                from pybedtools import BedTool
            except ImportError:
                raise ImportError("pybedtools is required. "
                                  "Install it or run the function with pybedtools=False"
                                  "Tip:pip install pybedtools")

            bed_sorted = BedTool(bed, from_string=True).sort()
            # bed_merged = bed_sorted.intersect(bed_sorted).merge(d=gap, c='4,5,6', o='distinct')
            bed_merged = bed_sorted.merge(d=gap, c='4,5,6', o='distinct')

            return _bed6_to_GeneAnnot(str(bed_merged))

        else:
            #ToDO: add argument to allow this to be done to cds_starts / cds_ends?

            intervals = [(s, e) for s, e in zip(self.starts, self.ends)]
            # sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
            # hopefully we dont need to sort this data... we may activate this for some messed up data, though
            merged = []

            for higher in intervals:  # sorted_by_lower_bound: # intervals already sorted.
                if not merged:
                    merged.append(higher)
                else:
                    lower = merged[-1]
                    # test for intersection between lower and higher:
                    # we know via sorting that lower[0] <= higher[0]
                    if higher[0] <= lower[1] + gap:
                        upper_bound = max(lower[1], higher[1])
                        merged[-1] = (lower[0], upper_bound)  # replace by merged interval
                    else:
                        merged.append(higher)

            starts, ends = zip(*merged)
            # self.starts = np.array(starts, dtype=np.int64)
            # self.ends = np.array(ends, dtype=np.int64)
            setattr(self, 'starts', np.array(starts, dtype=np.int64))
            setattr(self, 'ends', np.array(ends, dtype=np.int64))

        return self

    def _fix_orientation(self, orientation='Unknown'):
        # if orientation == 'transcript' and len(self) > 1:
        #
        #     if orientation == 'transcript':
        #         self._reverse()
        #
        #     elif orientation == 'Unknown':
        #         diff = self.starts[-1] - self.starts[0]
        #         if diff < 0:
        #             self._reverse()
        #
        # if orientation == 'Unknown' and len(self) > 1:
        if orientation != 'genomic' and self.blockCount() > 1:
            self.starts, self.ends = sort_intervals(self.starts, self.ends)
            if not np.isnan(np.sum(self.cds_starts)):
                self.cds_starts, self.cds_ends = sort_intervals(self.cds_starts, self.cds_ends)

    def _reverse(self):
        # print('reversing')
        self.starts = self.starts[::-1]
        self.ends = self.ends[::-1]

        if not np.isnan(np.sum(self.cds_starts)):
            self.cds_starts = self.cds_starts[::-1]
            self.cds_ends = self.cds_ends[::-1]

    def __str__(self):
        return '{} {}'.format(self.transcript_id, stringfy(self.exons))

    def format(self, format):
        if format == 'extb':
            chr_str = '%s:%s-%s' % (self.chrom, self.start + 1, self.end)
            extb = [
                # self.internal_exons,
                chr_str,
                self.strand,
                self.transcript_id,
                self.gene_id,
                self.blockCount(),
                sum(self.exons),
                self.exons,
                self.introns,
                # self.cds,
                # self.starts,
                ]

            return '\t'.join(stringfy(x) for x in extb)

        elif format == 'bed':
            # if np.isnan(np.sum(self.cds_starts)):
            #     thickStart = self.start
            #     thickStop = self.start
            #
            # else:
            #     thickStart = self.cds_starts[0]
            #     thickStop = self.cds_ends[-1]
            #
            if hasattr(self, 'itemRgb'):
                item_rgb = self.itemRgb
            else:
                item_rgb = "200,155,55"

            if self.thickStart and self.thickEnd:
                thickStart = self.thickStart
                thickEnd = self.thickEnd
            else:
                thickStart = thickEnd = self.start

            bed = [
                self.chrom,
                self.start,
                self.end,
                self.transcript_id,
                "1000",
                self.strand,
                thickStart,
                thickEnd,
                item_rgb,
                self.blockCount(),
                self.exons,
                self.starts - self.start
                ]

            return '\t'.join(stringfy(x) for x in bed)

        elif format == 'bed6':
            block_count = self.blockCount()
            bed = ['']*block_count

            for block in range(block_count):
                line = [
                    self.chrom,
                    self.starts[block],
                    self.ends[block],
                    self.transcript_id,
                    "1000",
                    self.strand]

                bed[block] = '\t'.join(stringfy(x) for x in line)

            return '\n'.join(bed) + '\n'

        else:
            super(InteractiveAnnotation, self).__format__(format)
