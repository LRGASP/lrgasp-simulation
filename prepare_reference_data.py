#!/usr/bin/env python
#
# ############################################################################
# Copyright (c) 2020 LRGASP consortium
# Prepares reference data in needed format and infers artifical novel isoforms
# ############################################################################

import sys
import os
from traceback import print_exc
import logging
import shutil
import random
import gffutils
from Bio import SeqIO
from Bio import Seq
import argparse
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger('LRGASP')

POLYA_LEN = 100


def parse_args(args=None, namespace=None):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--output", "-o", help="output prefix")
    parser.add_argument("--reference_annotation", "-r", help="reference annotation (GTF/.db)", type=str)
    parser.add_argument("--reference_transcripts", "-t", help="reference transcripts in FASTA format", type=str)
    parser.add_argument("--reference_genome", "-g", help="reference genome in FASTA format", type=str)
    parser.add_argument("--sqanti_prefix", "-q", help="path to SQANTI output "
                                                      "(_classification.txt and _corrected.gtf are needed)", type=str)
    parser.add_argument("--n_random_isoforms", "-n", help="infer this number of random isoforms into the annotation",
                        type=int)
    parser.add_argument("--isoform_list", "-l", help="infer only isoforms from a given file", type=str)
    parser.add_argument("--seed", "-s", help="randomizer seed [11]", default=11, type=int)

    args = parser.parse_args(args, namespace)

    if not check_params(args):
        parser.print_usage()
        exit(-1)
    return args


def check_params(args):
    if args.n_random_isoforms is not None and args.isoform_list is not None:
        logger.warning("Both --n_random_isoforms and --isoform_list are provided, only ones from the list will be used")
    return args.reference_annotation is not None and args.reference_transcripts is not None and args.sqanti_prefix is not None


def replace_gene_id(l, new_gene_id):
    gene_id_pos = l.find("gene_id")
    if gene_id_pos == -1:
        return l
    end_pos = l.find(";", gene_id_pos)

    return l[:gene_id_pos + len("gene_id")] + ' "' + new_gene_id + '"' + l[end_pos:]


def select_sqanti_isoforms(args):
    total_transcripts = 0
    # isofrom_id -> gene_ind
    novel_isoforms = {}
    logger.info("Loading SQANTI output")
    for l in open(args.sqanti_prefix + "_classification.txt"):
        total_transcripts += 1
        tokens = l.strip().split()
        isoform_id = tokens[0]
        isoform_type = tokens[5]
        gene_id = tokens[6]
        is_canonical = tokens[16] == 'canonical'

        if is_canonical and isoform_type in ['novel_not_in_catalog', 'novel_in_catalog']:
            novel_isoforms[isoform_id] = gene_id
    logger.info("Total isoforms read: %d, nic and nnic selected: %d" % (total_transcripts-1, len(novel_isoforms)))

    # selecting isoforms
    logger.info("Selecting novel isoforms")
    if args.isoform_list is not None:
        selected_isoform_set = set()
        for l in open(args.isoform_list):
            isoform_id = l.strip().split()[0]
            if isoform_id in novel_isoforms:
                selected_isoform_set.add(isoform_id)
            else:
                logger.warning("Specified isoform id %s is not found among canonical novel isoforms in SQANTI output" % isoform_id)
    elif args.n_random_isoforms is not None and args.n_random_isoforms < len(novel_isoforms):
        selected_isoform_set = set()
        all_novel_isoforms = list(novel_isoforms.keys())
        while len(selected_isoform_set) < args.n_random_isoforms:
            isoform_to_add = random.choice(all_novel_isoforms)
            selected_isoform_set.add(isoform_to_add)
    else:
        selected_isoform_set = set(novel_isoforms.keys())
    logger.info("Total isoform selected: %d" % len(selected_isoform_set))

    novel_isoform_file = args.output + ".novel_isoforms.tsv"
    novel_isoform_handler = open(novel_isoform_file, "w")
    for t_id in sorted(selected_isoform_set):
        novel_isoform_handler.write(t_id + "\n")
    novel_isoform_handler.close()

    # gene_id -> [(isoform_id, GTF lines)]
    novel_annotation = defaultdict(list)
    logger.info("Loading SQANTI annotation file")
    current_transcript_entry = ""
    current_transcrtip_id = ""

    for l in open(args.sqanti_prefix + "_corrected.gtf"):
        tokens = l.strip().split()
        enrty_type = tokens[2]

        if enrty_type == 'transcript':
            if current_transcrtip_id:
                # previous transcript was selected, adding to the annotation and clearing current entries
                gene_id = novel_isoforms[current_transcrtip_id]
                novel_annotation[gene_id].append((current_transcrtip_id, current_transcript_entry))
                current_transcript_entry = ""
                current_transcrtip_id = ""

            isoform_id = tokens[9].replace('"','').replace(';','')
            if isoform_id in selected_isoform_set:
                current_transcript_entry = replace_gene_id(l, novel_isoforms[isoform_id])
                current_transcrtip_id = isoform_id

        elif enrty_type == 'exon' and current_transcrtip_id:
            current_transcript_entry += replace_gene_id(l, novel_isoforms[current_transcrtip_id])

    if current_transcrtip_id:
        # last transcript was selected, adding to the annotation
        gene_id = novel_isoforms[current_transcrtip_id]
        novel_annotation[gene_id].append((current_transcrtip_id, current_transcript_entry))

    return novel_annotation


def load_gene_db(args):
    if args.reference_annotation.endswith(".db"):
        logger.info("Loading annotation from database")
        return gffutils.FeatureDB(args.reference_annotation, keep_order=True)

    logger.info("Converting gene annotation file to .db format (takes a while)...")
    db_file = args.output + ".reference_annotation.db"
    gffutils.create_db(args.reference_annotation, db_file, force=True, keep_order=True, merge_strategy='merge',
                       sort_attribute_values=True, disable_infer_transcripts=True, disable_infer_genes=True)
    logger.info("Gene database written to " + db_file)
    logger.info("Provide this database next time to avoid excessive conversion")

    logger.info("Loading annotation from database")
    return gffutils.FeatureDB(db_file, keep_order=True)


def infer_novel_transcripts_to_gtf(args, genedb, novel_annotation):
    extended_annotation_file = args.output + ".annotation.gtf"
    logger.info("Saving extended gene annotation (takes a while)...")
    with open(extended_annotation_file, "w") as f:
        for record in genedb.all_features():
            f.write(str(record) + '\n')
            if record.featuretype == 'gene' and record.id in novel_annotation:
                # add novel isoforms
                for t in novel_annotation[record.id]:
                    f.write(t[1])
    logger.info("Extended annotation saved to %s" % extended_annotation_file)


def parse_transcript(gtf_fragment):
    lines = gtf_fragment.split('\n')
    tokens = lines[0].split()
    assert tokens[2] == 'transcript'

    chr_id = tokens[0]
    strand = tokens[6]

    exons = []
    for l in lines[1:]:
        if not l:
            continue
        tokens = l.split()
        exons.append((int(tokens[3]), int(tokens[4])))

    return chr_id, strand, sorted(exons)


def extract_transcript_from_fasta(chr_seq, strand, exons):
    transcript = ""
    for e in exons:
        transcript += str(chr_seq[e[0]-1:e[1]])

    if strand == '-':
        transcript = str(Seq(transcript).reverse_complement()).upper()

    transcript += "A" * POLYA_LEN
    return transcript


def infer_novel_transcripts_to_fasta(args, novel_annotation):
    logger.info("Loading reference genome from " + args.reference_genome)
    genome_records = SeqIO.index(args.reference_genome, "fasta")
    # chrid -> [(t_id, strand, exons)]
    genomic_regions = defaultdict(list)
    for gene_id in novel_annotation:
        for t in novel_annotation[gene_id]:
            chr_id, strand, exons = parse_transcript(t[1])
            genomic_regions[chr_id].append((t[0], strand, exons))

    new_transcripts = []
    for chr_id in genomic_regions:
        chr_seq = genome_records[chr_id]
        for transcript_tuple in genomic_regions[chr_id]:
            transcript_seq = extract_transcript_from_fasta(chr_seq, transcript_tuple[1], transcript_tuple[2])
            record = SeqRecord(Seq(transcript_seq), id=transcript_tuple[0])
            new_transcripts.append(record)

    for record in SeqIO.parse(args.reference_transcripts, "fasta"):
        record.seq = Seq(str(record.seq) + "A" * POLYA_LEN)
        new_transcripts.append(record)

    SeqIO.write(new_transcripts, args.output + ".transcripts.fasta", 'fasta')


def set_logger(logger_instance):
    logger_instance.setLevel(logging.INFO)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    logger_instance.addHandler(ch)


def run_pipeline(args):
    logger.info(" === LRGASP reference data preparation started === ")
    random.seed(args.seed)
    novel_annotation = select_sqanti_isoforms(args)
    gene_db = load_gene_db(args)
    infer_novel_transcripts_to_gtf(args, gene_db, novel_annotation)
    infer_novel_transcripts_to_fasta(args, novel_annotation)
    shutil.copyfile(args.reference_genome, args.output + ".transcripts.fasta")
    logger.info(" === LRGASP reference data preparation finished === ")


def main(args):
    args = parse_args(args)
    set_logger(logger)
    out_dir = os.path.dirname(args.output)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
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
