# ############################################################################
# Copyright (c) 2020 LRGASP consortium
# ############################################################################

import logging
import subprocess
import os

logger = logging.getLogger('LRGASP')


def simulate_ont(args, read_count=1000):
    src_path = os.path.dirname(os.path.realpath(__file__))
    nanosim = os.path.join(src_path, "NanoSim/src/simulator.py")
    model_dir = os.path.join(src_path, "NanoSim/pre-trained_models/")
    if args.ont_type == 'dRNA':
        model_name = 'human_NA12878_dRNA_Bham1_guppy'
        model_pref = model_dir+model_name+'/'
        molecule_type = 'dRNA'
        uracil = True
    elif args.ont_type == 'cDNA':
        model_name = 'human_NA12878_cDNA_Bham1_guppy'
        model_pref = model_dir + model_name + '/'
        molecule_type = 'cDNA_1D2'
        uracil = False
    else:
        logger.error("Unexpected ont_type value %s" % args.ont_type)
        return

    # untar the model
    if not os.path.exists(model_pref):
        logger.info("Unpacking NanoSim model")
        cmd = ['tar', '-xzf', model_name+'.tar.gz']
        current_wd = os.getcwd()
        os.chdir(model_dir)
        result = subprocess.run(cmd)
        os.chdir(current_wd)
        if result.returncode != 0:
            logger.error("Untarring NanoSim training model failed.")

    logger.info("Simulating ONT reads with Trans-NanoSim...")
    ref_prefix = args.reference_prefix
    cmd = [nanosim, 'transcriptome',
                             "-rg", ref_prefix + ".genome.fasta",
                             "-rt", ref_prefix + ".transcripts.fasta",
                             "-e", str(args.counts),
                             "-t", str(args.threads),
                             "--seed", str(args.seed),
                             "-b", "guppy",
                             "-o", os.path.join(args.output, "ONT.simulated"),
                             "-n", str(read_count),
                             "--fastq",
                             "-c", model_pref + 'training',
                             "-r", molecule_type, "--no_model_ir"]
    if uracil:
        cmd += ['--uracil']

    result = subprocess.run(cmd)

    if result.returncode != 0:
        logger.error("NanoSim failed, contact developers for support.")
        return

    # TODO - remove isoform IDs from output and rewrite fastas

    logger.info("Done.")
