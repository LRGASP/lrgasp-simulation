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
    if args.species == 'human':
        model = model_dir+'human_NA12878_dRNA_Bham1_guppy.tar.gz'
        model_pref = model_dir+'human_NA12878_dRNA_Bham1_guppy/'
        molecule_type = 'dRNA'
    elif args.species == 'mouse':
        model = model_dir+'human_NA12878_cDNA_Bham1_guppy.tar.gz'
        model_pref = model_dir+'human_NA12878_cDNA_Bham1_guppy/'
        molecule_type = 'cDNA_1D2'

    # untar the model
    cmd = ['tar', '-xf', model]
    result = subprocess.run(cmd)
    if result.returncode != 0:
        logger.error("Untarring NanoSim training model failed.")

    ref_prefix = os.path.join(args.reference_dir, args.reference_name)
    result = subprocess.run([nanosim, 'transcriptome',
                             "-rg", ref_prefix + ".genome.fasta",
                             "-rt", ref_prefix + ".transcripts.fasta",
                             "-e", args.counts,
                             "-t", args.threads,
                             "--seed", args.seed,
                             "-b", "guppy",
                             "-o", os.path.join(args.output, "ONT.simulated")])
                             "-n", read_count,
                             "--polya",
                             "-c", "model_pref",
                             "-r", molecule_type
                             )

    # do we also want to simulate with U instead of T (--uracil option) when
    # simulating dRNA?
    # does the output need to be in fasta or fastq format?


    if result.returncode != 0:
        logger.error("NanoSim failed, contact developers for support.")
        return

    logger.info("Done.")
