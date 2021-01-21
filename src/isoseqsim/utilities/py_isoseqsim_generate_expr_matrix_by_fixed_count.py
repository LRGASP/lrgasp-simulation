#!/usr/bin/env python3

# NOTE: this file was created by the LRGASP consortium and is not included the original IsoSeqSim source code
# Contributor: Andrey Prjibelski

import sys
import time
import argparse
from collections import defaultdict

def main(args):
    sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
    sys.stdout.flush()
    generate_expr_matrix(args.input, args.input_expr, args.output)
    sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
    sys.stdout.flush()


def generate_expr_matrix(input_gpd_fl, input_txt_fl, output_expr_mtx):
    # parse txt
    dic_iso_expr = defaultdict(int)
    for line in input_txt_fl:
        if line.startswith("#"):
            continue
        vals = line.strip().split("\t")
        iso_id = vals[0]
        expr_v = vals[2]  # expr_v is TPM now
        # we aim to generate 10 million reads per run
        dic_iso_expr[iso_id] = int(round(10 * float(expr_v)))

    for line in input_gpd_fl:
        iso_id = line.strip().split("\t")[1]
        print(line.strip() + "\t" + str(dic_iso_expr[iso_id]), file=output_expr_mtx)

    input_txt_fl.close()
    input_gpd_fl.close()
    output_expr_mtx.close()


def do_inputs():
    parser = argparse.ArgumentParser(
        description="Randomly generate read count for each isoform based on negative binomial (NB) distribution. Read count is shown in last column of output file",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=argparse.FileType('r'), required=True, help="Input: gpd file")
    parser.add_argument('-e', '--input_expr', type=argparse.FileType('r'), required=True,
                        help="Input: expression txt file (first colunm is isoform ID, second colunmn is expression value)")
    parser.add_argument('-o', '--output', type=argparse.FileType('w'), required=True,
                        help="Output: gpd + read count file")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = do_inputs()
    main(args)
