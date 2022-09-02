#!/usr/bin/env python3

"""Module containing the function to extract HSPs from a list of proteins

Copyright (c) 2022, François Charih
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

from loguru import logger
import Bio.SeqIO

import rsprint

def main(args):
    logger.info("Loading the protein sequences...")
    proteins = [(p.id, str(p.seq)) for p in Bio.SeqIO.parse(args.input, "fasta")]

    logger.info(f"Extracting the HSPs from the {len(proteins)} provided protein sequences...")
    hsps = rsprint.extract_hsps(proteins)

    logger.info(f"Processing the HSPs to account for overrepresented residues...")
    processed_hsps = rsprint.process_hsps(proteins, hsps, kmer_size=args.kmer_size, t_count=args.t_count)

    logger.info(f"Writing to a file...")
    with open(args.output, "w") as output_file:
        output_file.write("\n".join([f"{x[0]} {x[1]} {x[2]} {x[3]} {x[4]}" for x in processed_hsps]))

    logger.info(f"Done!")