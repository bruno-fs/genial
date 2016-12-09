import re

from .GenomeAnnotation import InteractiveAnnotation
from .utils import str2array

input_formats = {'gff3', 'gtf', 'bed'}
output_formats = {'bed', 'extb'}



def bed12_to_GeneAnnot(bed12):
    """
    01) chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
    02) chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
    03) chromEnd - The ending position of the feature in the chromosome or scaffold.
    04) name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
    05) score - A score between 0 and 1000.
    06) strand - Defines the strand - either '+' or '-'.
    07) thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
    08) thickEnd - The ending position at which the feature is drawn thickly (for example, the stop codon in gene displays).
    09) itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
    10) blockCount - The number of blocks (exons) in the BED line.
    11) blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
    12) blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.
    """

    bed12 = bed12.strip()
    bed_field = bed12.split('\t')
    starts = str2array(bed_field[11]) + int(bed_field[1])
    ends = str2array(bed_field[10]) + starts
    name = bed_field[3]
    chrom = bed_field[0]
    strand = bed_field[5]
    thickStart = bed_field[6]
    thickEnd = bed_field[7]
    itemRgb = bed_field[8]

    return InteractiveAnnotation(starts, ends, strand,
                                 chrom=chrom,
                                 transcript_id=name,
                                 thickStart=thickStart,
                                 thickEnd=thickEnd,
                                 starts_offset=0,
                                 itemRgb=itemRgb)


def parse(file_handle, format):
    """

    Parameters
    ----------
    file_handle
    format

    Returns
    -------

    """
    # --------------- gff3 / gtf -----------------------------
    if format in ['gff3', 'gtf']:
        from .gff import parse_to_dict
        gff = parse_to_dict(file_handle, format)

        for tranx in gff:
            # ToDo: support fivePrime and threePrime UTR elements
            if re.match(r'\s*$', gff[tranx].exon_starts):
                gff[tranx].exon_starts = gff[tranx].CDS_starts
                gff[tranx].exon_ends = gff[tranx].CDS_ends

            annotation = InteractiveAnnotation(
                starts=gff[tranx].exon_starts,
                ends=gff[tranx].exon_ends,
                strand=gff[tranx].strand,
                orientation=gff.orientation,
                cds_starts=gff[tranx].CDS_starts,
                cds_ends=gff[tranx].CDS_ends,
                chrom=gff[tranx].chrom,
                transcript_id=gff[tranx].transcript_id,
                gene_id=gff[tranx].gene_id)

            yield annotation

    # -------------------- bed -----------------------------
    elif format == 'bed':
        for line in file_handle:
            yield bed12_to_GeneAnnot(line)
