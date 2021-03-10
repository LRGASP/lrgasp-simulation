# ############################################################################
# Copyright (c) 2020 LRGASP consortium
# ############################################################################

import logging
import subprocess
import os
import pandas as pd

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
    if not args.noise_reads:
        cmd.append("--aligned_only")

    result = subprocess.run(cmd)

    if result.returncode != 0:
        logger.error("NanoSim failed, contact developers for support.")
        return

    # get the real transcript IDs because NanoSim butchers them
    logger.info("Reading GTF annotation to restore original isoform names...")
    gtf = ref_prefix + ".annotation.gtf"
    df = pd.read_csv(gtf, sep='\t',  header=None, usecols=[2,8])
    df.columns = ['entry_type', 'fields']
    df = df.loc[df.entry_type == 'transcript']
    df['tid'] = df.fields.str.split(pat='transcript_id "', expand=True)[1]
    df['tid'] = df.tid.str.split(pat='";', n=1, expand=True)[0]
    df['nanosim_tid'] = df.tid.str.split(pat='.', expand=True)[0]
    df = df[['tid', 'nanosim_tid']]

    read_tid_map = {}
    fastqs = [os.path.join(args.output, "ONT.simulated_aligned_reads.fastq")]
    if args.noise_reads:
        fastqs.append(os.path.join(args.output, "ONT.simulated_unaligned_reads.fastq"))
    read_num = 0
    # dump everythin the same file
    fname = os.path.join(args.output, 'ONT.simulated.fastq')
    ofile = open(fname, 'w')

    logger.info("Renaming and counting ONT reads...")
    for f in fastqs:
        ifile = open(f, 'r')
        for line in ifile:
            if line.startswith('@'):

                # get transcript id from each simulated read
                tid = line[1:]
                tid = tid.split('_')[0]

                # convert to non-butchered tid and add to read map
                correct_tid = df.loc[df.nanosim_tid == tid, 'tid'].tolist()
                if not correct_tid:
                    logger.warning("%s was not found in the annotation, skipping" % tid)
                    correct_tid = tid
                else:
                    correct_tid = correct_tid[0]
                read_id = 'ONT_simulated_read_{}'.format(read_num)
                read_num += 1
                read_tid_map[read_id] = correct_tid

                ofile.write('@{}\n'.format(read_id))
            else:
                ofile.write(line)
        ifile.close()

    nanosim_aux_dir = os.path.join(args.output, "NanoSim_data")
    os.makedirs(nanosim_aux_dir)
    for aux_file_name in ["ONT.simulated_aligned_reads.fastq", "ONT.simulated_unaligned_reads.fastq",
                          "ONT.simulated_aligned_error_profile", "ONT.simulated_unaligned_error_profile"]:
        aux_path = os.path.join(args.output, aux_file_name)
        if os.path.exists(aux_path):
            os.rename(aux_path, os.path.join(nanosim_aux_dir, aux_file_name), )

    ofile.close()

    # write read to tid map
    logger.info("Saving de facto counts and read-to-isoform map")
    df = pd.DataFrame.from_dict(read_tid_map, orient='index')
    df.reset_index(inplace=True)
    df.columns = ['read_id', 'tid']
    fname = os.path.join(args.output, 'ONT.simulated.read_to_isoform.tsv')
    df.to_csv(fname, sep='\t', header=None, index=False)

    # compute / write counts per isoform
    df = df.groupby('tid').count()
    df.reset_index(inplace=True)
    df.rename({'read_id': 'counts'}, axis=1, inplace=True)
    fname = os.path.join(args.output, 'ONT.simulated.isoform_counts.tsv')
    df.to_csv(fname, sep='\t', header=None, index=False)

    logger.info("Done.")
