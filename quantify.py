#!/usr/bin/env python
#
# ############################################################################
# Copyright (c) 2020 LRGASP consortium
# ############################################################################

import sys
from traceback import print_exc
import logging
import argparse

logger = logging.getLogger('LRGASP')


def parse_args(args=None, namespace=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--output", "-o", help="output file with counts")
    parser.add_argument("--fastq", "-f", help="long reads in FASTQ format (PacBio/ONT)", type=str)
    parser.add_argument("--reference_transcripts", "-r", help="reference transcripts in FASTA format", type=str)
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
