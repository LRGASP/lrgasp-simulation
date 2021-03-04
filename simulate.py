#!/usr/bin/env python
#
# ############################################################################
# Copyright (c) 2020 LRGASP consortium
# ############################################################################

import os
import sys
from traceback import print_exc
import logging
import argparse
import shutil
import random
import numpy as np


from src.simulate_pacbio import *

# TODO insert correct numbers
PACBIO_READ_COUNT = 1000


logger = logging.getLogger('LRGASP')


def parse_args(args=None, namespace=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--output", "-o", help="output folder, will be created automatically [default=lrgasp_simulation]",
                        type=str, default="lrgasp_simulation")
    parser.add_argument("--reference_prefix", "-r", help="a path to reference files",  type=str)
    parser.add_argument("--counts", "-c", help="transcript abundances in TSV format (output of quantify.py)", type=str)
    parser.add_argument("--threads", "-t", help="number of CPU threads for simulators [16]", default=16, type=int)
    parser.add_argument("--seed", "-s", help="randomizer seed [11]", default=11, type=int)
    parser.add_argument("--species", help="species to simulate from 'human' or 'mouse'", default="human", type=str)

    args = parser.parse_args(args, namespace)

    if os.path.exists(args.output):
        # logger is not defined yet
        print("WARNING! Output folder already exists, some files may be overwritten")
    else:
        os.makedirs(args.output)

    if args.species not in ['human', 'mouse']:
        raise ValueError("--species value must be 'human' or 'mouse'")

    args.tmp_dir = os.path.join(args.output, "tmp")
    if not os.path.exists(args.output):
        os.makedirs(args.tmp_dir)

    if not check_params(args):
        parser.print_usage()
        exit(-1)
    return args


def check_params(args):
    return args.counts is not None


def set_logger(args, logger_instance):
    logger_instance.setLevel(logging.INFO)
    log_file = os.path.join(args.output, "simulator.log")
    f = open(log_file, "w")
    f.write("CMD: " + ' '.join(sys.argv) + '\n')
    f.close()
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger_instance.addHandler(fh)
    logger_instance.addHandler(ch)


def run_pipeline(args):
    logger.info(" === LRGASP simulation pipeline started === ")
    random.seed(args.seed)
    np.random.seed(args.seed)
    # simulate short reads

    # simulate PacBio reads
    simulate_pacbio(args, PACBIO_READ_COUNT)
    # simulate ONT reads

    shutil.rmtree(args.tmp_dir)
    logger.info(" === LRGASP simulation pipeline finished === ")


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
