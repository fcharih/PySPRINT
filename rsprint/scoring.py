#!/usr/bin/env python3
"""Module containing the function to score interactions

Copyright (c) 2022, François Charih
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

from typing import Set, Tuple
import Bio.SeqIO
from loguru import logger 

import rsprint

from .predictions import Predictions

def load_hsps(filename: str) -> Set[Tuple[str, str, int, int, int]]:
    hsps = set()
    for hsp in open(filename).read().split("\n"):
        split = hsp.split()
        p1 = split[0]
        p2 = split[1]
        l1 = int(split[2])
        l2 = int(split[3])
        length = int(split[4])
        hsps.add((p1, p2, l1, l2, length))
    return hsps

def score_all_to_all(args):
    logger.info("Loading the protein sequences...")
    proteins = [(p.id, str(p.seq)) for p in Bio.SeqIO.parse(args.input, "fasta")]

    logger.info(f"Loading the HSPs...")
    hsps = load_hsps(args.hsps)

    logger.info(f"Loading the training pairs...")
    training_pairs = [tuple(x.rstrip("\n").split()) for x in open(args.training_pairs)]

    logger.info(f"Scoring the interactions...")
    scores = rsprint.score_interactions(proteins, hsps, training_pairs, kmer_size=args.kmer_size)
    predictions = Predictions(scores, [p[0] for p in proteins])

    logger.info(f"Saving the scores to a {args.output}...")
    predictions.save(args.output)

    logger.info(f"Done!")

def score_peptides(args):
    logger.info("Loading the protein sequences...")
    proteins = [(p.id, str(p.seq)) for p in Bio.SeqIO.parse(args.input, "fasta")]

    logger.info("Loading the peptide sequences...")
    peptides = [(p.id, str(p.seq)) for p in Bio.SeqIO.parse(args.peptides, "fasta")]

    logger.info(f"Loading the HSPs...")
    hsps = load_hsps(args.hsps)

    logger.info(f"Loading the training pairs...")
    training_pairs = [tuple(x.rstrip("\n").split()) for x in open(args.training_pairs)]

    logger.info(f"Extracting the peptide hsps...")
    peptide_hsps = rsprint.extract_peptide_hsps(proteins, peptides, t_smer=args.t_smer, t_hsp=args.t_hsp, kmer_size=args.kmer_size)
    hsps = hsps.union(peptide_hsps)

    logger.info(f"Scoring the interactions...")
    scores = rsprint.score_interactions(proteins, hsps, training_pairs, kmer_size=args.kmer_size)
    predictions = Predictions(scores, [p[0] for p in proteins], peptide_names=[p[0] for p in peptides])

    logger.info(f"Saving the scores to a {args.output}...")
    predictions.save(args.output)

    logger.info(f"Done!")