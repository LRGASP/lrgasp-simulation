#!/usr/bin/env python
#
# ############################################################################
# Copyright (c) 2020 LRGASP consortium
# ############################################################################

import sys
import os
from traceback import print_exc
import logging
import argparse
import subprocess
import pysam
import random
from collections import defaultdict

logger = logging.getLogger('LRGASP')

# list of "novel" transcript that will be added to the annotation / transcriptome for testing purposes
MIN_NOVEL_COUNT = 10
MAX_NOVEL_COUNT = 1000


def parse_args(args=None, namespace=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--output", "-o", help="output file with counts")
    parser.add_argument("--fastq", "-f", help="long reads in FASTQ format (PacBio/ONT)", type=str)
    parser.add_argument("--threads", "-t", help="number of CPU threads for minimap [16]", default=16, type=int)
    parser.add_argument("--reference_transcripts", "-r", help="reference transcripts in FASTA format", type=str)
    parser.add_argument("--mandatory", "-m", help="file with a list of mandatory transcripts to be included, "
                                                  "counts are assigned randomly", type=str)
    args = parser.parse_args(args, namespace)

    if not check_params(args):
        parser.print_usage()
        exit(-1)
    return args


def check_params(args):
    return args.fastq is not None and args.reference_transcripts is not None


def set_logger(args, logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def run_pipeline(args):
    logger.info(" === LRGASP quantification pipeline started === ")
    logger.info("Mapping reads with minimap2...")
    samfile_name = args.output + ".sam"
    result = subprocess.run(["minimap2", args.reference_transcripts, args.fastq, "-x", "map-ont",
                             "-t", str(args.threads), "-a", "--secondary=no", "-o", samfile_name])
    if result.returncode != 0:
        logger.error("minimap2 failed with code %d, make sure it is in you $PATH variable" % result.returncode)
        exit(-1)

    logger.info("Quantifying transcripts...")
    transcript_counts = defaultdict(int)
    with pysam.AlignmentFile(samfile_name, "r") as samfile_in:
        for alignment in samfile_in:
            transcript_id = alignment.reference_name
            if alignment.reference_id == -1 or alignment.is_supplementary or alignment.is_secondary:
                continue
            transcript_counts[transcript_id] += 1
    os.remove(samfile_name)

    mandatory_transcripts = []
    if args.mandatory is not None:
        for l in open(args.mandatory):
            if not l or l.startswith("#"):
                continue
            mandatory_transcripts.append(l.strip())

    # adding "novel" transcript that must be in the simulated data
    for transcript_id in mandatory_transcripts:
        # we have a few "novel" transcripts, so we don't really care if their abundance is uniformly distributed,
        # which is not biologically sound
        transcript_counts[transcript_id] = random.randint(MIN_NOVEL_COUNT, MAX_NOVEL_COUNT)

    count_sum = 0.0
    for count in transcript_counts.values():
        # we don't devide by gene length here as we assume each long read is a single transcripts (unlike in short reads)
        count_sum += count
    scaling_factor = count_sum / 1000000.0

    outf = open(args.output, "w")
    outf.write("#transcript_id\tcounts\ttpm\n")
    for transcript_id in transcript_counts.keys():
        outf.write("%s\t%d\t%0.6f\n" % (transcript_id, transcript_counts[transcript_id],
                                        transcript_counts[transcript_id] / scaling_factor))
    outf.close()
    logger.info("Done. Output counts are stored in " + args.output)
    logger.info(" === LRGASP quantification pipeline finished === ")


def main(args):
    args = parse_args(args)
    set_logger(args, logger)
    run_pipeline(args)


if __name__ == "__main__":
    # stuff only to run when not called via 'import' here
    try:
        main(sys.argv[1:])
    except SystemExit:
        raise
    except:
        print_exc()
        sys.exit(-1)
