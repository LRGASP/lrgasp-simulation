# ############################################################################
# Copyright (c) 2020 LRGASP consortium
# ############################################################################

import logging
import subprocess
import os

logger = logging.getLogger('LRGASP')


def simulate_pacbio(args, read_count=1000):
    logger.info("Simulating PacBio reads with IsoSeqSim...")
    src_path = os.path.dirname(os.path.realpath(__file__))
    isoseqsim = os.path.join(src_path, "isoseqsim/bin/isoseqsim.py")
    param_dir = os.path.join(src_path, "isoseqsim/utilities/")
    ref_prefix = args.reference_prefix
    result = subprocess.run([isoseqsim, "--cpu", str(args.threads), "--tempdir", args.tmp_dir,
                             "--annotation", ref_prefix + ".annotation.gtf",
                             "--genome", ref_prefix + ".genome.fasta", "--expr", args.counts,
                             "--c5", os.path.join(param_dir, "5_end_completeness.PacBio-Sequel.tab"),
                             "--c3", os.path.join(param_dir, "3_end_completeness.PacBio-Sequel.tab"),
                             "--es", "0.004", "--ei", "0.006", "--ed", "0.006",
                             "--read_number", str(read_count / 1000000.0), "--polya",
                             "--transcript", os.path.join(args.output, "PacBio.simulated.tsv"),
                             "-o", os.path.join(args.output, "PacBio.simulated")])

    if result.returncode != 0:
        logger.error("IsoSeqSim failed, contact developers for support.")
        return

    logger.info("Done.")
