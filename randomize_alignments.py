import pysam
from collections import defaultdict
import numpy as np

def parse_bam_by_read(bam_filepath):
    """Parse multimapping reads from BAM file into an in-memory Python dictionary.

    """
    
    samfile = pysam.AlignmentFile(bam_filepath, "rb")
    reads = samfile.fetch()
    read_to_alignments = defaultdict(list)

    for read in reads:
        read_to_alignments[read.qname].append(read)

    return read_to_alignments
    
def permute_primary_flags(read_to_alignments):
    """Randomly permute primary flags

    """
    new_read_to_alignments = defaultdict(list)
    for name in read_to_alignments:
        alignments = read_to_alignments[name]
        primary_index = -1
        next_primary_index = -1

        for index in np.random.permutation(len(alignments)):
            alignment = alignments[index]
            is_primary = (bin(alignment.flag >> 8)[-1] == '0')

            if is_primary:
                primary_index = index
            elif (next_primary_index == -1):
                next_primary_index = index

        new_alignments = alignments.copy()
        if (primary_index != -1) and (next_primary_index != -1):
            # Permutation necessary
            new_alignments[primary_index].flag += 0x100
            new_alignments[next_primary_index].flag -= 0x100

        new_read_to_alignments[name] = new_alignments

    return new_read_to_alignments

def write_reads_to_bam(bam_filepath, read_to_alignments):
    """Write reads from multimapping dictionary to BAM file.

    """

    bamfile = pysam.AlignmentFile(bam_filepath, "wb")

    for read in read_to_alignments:
        for alignment in read_to_alignments[read]:
            bamfile.write(alignment)
