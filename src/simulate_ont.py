# ############################################################################
# Copyright (c) 2020 LRGASP consortium
# ############################################################################

import logging
import subprocess
import os

logger = logging.getLogger('LRGASP')


def simulate_ont(args, read_count=1000):
    logger.info("Simulating ONT reads...")
    src_path = os.path.dirname(os.path.realpath(__file__))
    nanosim = os.path.join(src_path, "NanoSim/src/simulator.py")
    model_dir = os.path.join(src_path, "NanoSim/pre-trained_models/")
    ref_prefix = os.path.join(args.reference_dir, args.reference_name)
    result = subprocess.run([nanosim, 'transcriptome',
                             "-rg", ref_prefix + ".genome.fasta",
                             "-rt", ref_prefix + ".transcripts.fasta",
                             "-e", args.counts,
                             "-t", str(args.threads),
                             "--seed", args.seed,
                             "-o", os.path.join(args.output, "ONT.simulated")])
                             "-n", read_count,
                             "--polya",
                             "-c", "TODO , whatever model",
                             

    if result.returncode != 0:
        logger.error("NanoSim failed, contact developers for support.")
        return

    logger.info("Done.")
