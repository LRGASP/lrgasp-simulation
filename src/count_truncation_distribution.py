import sys
import pysam
from traceback import print_exc
import numpy as np


def main():
    bamfile = pysam.AlignmentFile(sys.argv[1], "rb")
    left_truncations = []
    right_truncations = []
    for alignment in bamfile:
        if alignment.is_secondary or alignment.is_supplementary or alignment.is_unmapped:
              continue

        transcript_len = alignment.get_reference_length(alignment.reference_id)
        left_truncations.append(float(alignment.reference_start) / float(transcript_len))
        right_truncations.append(float(transcript_len - alignment.reference_end) / float(transcript_len))

    bins = [0.0] + [0.005 + 0.01 * i for i in range(0,100)] + [1.0]
    print(len(bins), bins)
    res_left = np.histogram(left_truncations, bins=bins)
    print(res_left[1])
    res_right = np.histogram(left_truncations, bins=bins)
    print(res_right[1])


if __name__ == "__main__":
   # stuff only to run when not called via 'import' here
    try:
        main()
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)