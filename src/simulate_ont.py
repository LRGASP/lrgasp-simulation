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
        model_name = 'human_NA12878_dRNA_Bham1_guppy'
        model_pref = model_dir+model_name+'/'
        molecule_type = 'dRNA'
    elif args.species == 'mouse':
        model_name = 'human_NA12878_cDNA_Bham1_guppy'
        model_pref = model_dir + model_name + '/'
        molecule_type = 'cDNA_1D2'
    else:
        logger.error("Unexpected species value %s" % args.species)
        return

    # untar the model
    if not os.path.exists(model_pref):
        logger.info("Unpackking NanoSim model")
        cmd = ['tar', '-xzf', model_name+'.tar.gz']
        current_wd = os.getcwd()
        os.chdir(model_dir)
        result = subprocess.run(cmd)
        os.chdir(current_wd)
        if result.returncode != 0:
            logger.error("Untarring NanoSim training model failed.")

    logger.info("Simulating ONT reads with Trans-NanoSim...")
    ref_prefix = args.reference_prefix
    result = subprocess.run([nanosim, 'transcriptome',
                             "-rg", ref_prefix + ".genome.fasta",
                             "-rt", ref_prefix + ".transcripts.fasta",
                             "-e", str(args.counts),
                             "-t", str(args.threads),
                             "--seed", str(args.seed),
                             "-b", "guppy",
                             "-o", os.path.join(args.output, "ONT.simulated"),
                             "-n", str(read_count),
                             "-c", model_pref + 'training',
                             "-r", molecule_type, "--no_model_ir"])

    # do we also want to simulate with U instead of T (--uracil option) when
    # simulating dRNA?
    # does the output need to be in fasta or fastq format?


    if result.returncode != 0:
        logger.error("NanoSim failed, contact developers for support.")
        return

    logger.info("Done.")
