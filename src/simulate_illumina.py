# ############################################################################
# Copyright (c) 2020 LRGASP consortium
# ############################################################################

import logging
import subprocess
import os
from Bio import SeqIO
from collections import defaultdict

logger = logging.getLogger('LRGASP')


def prepare_rsem_counts(ref_transcripts, counts_file, out_file):
    tpm_dict = defaultdict(float)
    for l in open(counts_file):
        t = l.strip().split()
        tpm_dict[t[0]] = float(t[1])

    outf = open(out_file, "w")
    outf.write("transcript_id\tgene_id\tlength\teffective_length\texpected_count\tTPM\tFPKM\tIsoPct\n")
    for transcirpt in SeqIO.index(ref_transcripts, "fasta"):
        seq_len = len(transcirpt)
        eff_len = float(max(1, seq_len - 220))
        tpm = tpm_dict[transcirpt.id]
        iso_pct = 100.0 if tpm > 0.0 else 0.0
        outf.write("%s\t%s\t%d\t%.2f\t0.00\t%.2f\t0.00\t%.2f\n" %
                   (transcirpt.id, transcirpt.id, seq_len, eff_len, tpm, iso_pct))
    outf.close()


def simulate_illumina(args, read_count=100000):
    logger.info("Simulating Illumina reads...")
    src_path = os.path.dirname(os.path.realpath(__file__))
    rsem_dir = os.path.join(src_path, "RSEM")

    if not os.path.exists(os.path.join(rsem_dir, "rsem-simulate-reads")):
        # compile rsem
        current_wd = os.getcwd()
        os.chdir(rsem_dir)
        result = subprocess.run(["make"])
        os.chdir(current_wd)
        if result.returncode != 0:
            logger.error("Compilation of RSEM simulator failed, contact developers for support.")
            return
    assert os.path.exists(os.path.join(rsem_dir, "rsem-simulate-reads"))

    ref_prefix = args.reference_prefix
    rsem_data_dir = os.path.join(args.output, "rsem_data")
    os.makedirs(rsem_data_dir)
    prepare_reference = os.path.join(rsem_dir, "rsem-prepare-reference")
    rsem_ref_path = os.path.join(rsem_data_dir, "RSEM.reference")
    ref_transcripts = ref_prefix + ".transcripts.fasta"
    result = subprocess.run([prepare_reference, ref_transcripts, rsem_ref_path])
    if result.returncode != 0:
        logger.error("RSEM reference data preparation failed, contact developers for support.")
        return

    rsem_count_file = os.path.join(rsem_data_dir, "RSEM.counts.tsv")
    prepare_rsem_counts(ref_transcripts, args.counts, rsem_count_file)

    rsem_simulate = os.path.join(rsem_dir, "rsem-simulate-reads")
    model_file = os.path.join(rsem_dir, "models/Illumina150.RNA.model")
    result = subprocess.run([rsem_simulate, rsem_ref_path, model_file, rsem_count_file, "0.0", str(read_count),
                             os.path.join(args.output, "Illumina.simulated"), "--seed", str(args.seed)])
    if result.returncode != 0:
        logger.error("RSEM simulator failed, contact developers for support.")
        return

    logger.info("Done.")
