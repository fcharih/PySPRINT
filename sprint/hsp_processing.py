#!/usr/bin/env python3

"""Module containing the function to process a set of HSPs

Copyright (c) 2022, François Charih
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""
from typing import Tuple

from loguru import logger
import Bio.SeqIO

import sprint

def parse_hsp(line) -> Tuple[str, str, int, int, int]:
    split = line.split()
    return (split[0], split[1], int(split[2]), int(split[3]), int(split[4]))

def main(args):
    logger.info("Loading the protein sequences...")
    proteins = [(p.id, str(p.seq)) for p in Bio.SeqIO.parse(args.input, "fasta")]

    logger.info(f"Loading the HSPs...")
    hsps = set([parse_hsp(line) for line in open(args.hsps).read().splitlines()])

    logger.info(f"Processing the HSPs to account for overrepresented residues...")
    processed_hsps = sprint.process_hsps(proteins, hsps, kmer_size=args.kmer_size, t_count=args.t_count)

    logger.info(f"Writing to a file...")
    with open(args.output, "w") as output_file:
        output_file.write("\n".join([f"{x[0]} {x[1]} {x[2]} {x[3]} {x[4]}" for x in processed_hsps]))

    logger.info(f"Done!")
