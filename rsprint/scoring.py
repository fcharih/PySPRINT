#!/usr/bin/env python3
"""Module containing the function to score interactions

Copyright (c) 2022, François Charih
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

import Bio.SeqIO
from loguru import logger 

import rsprint

from .predictions import Predictions

def score_all_to_all(args):
    logger.info("Loading the protein sequences...")
    proteins = [(p.id, str(p.seq)) for p in Bio.SeqIO.parse(args.input, "fasta")]

    logger.info(f"Loading the HSPs...")
    hsps = set([tuple(x.split()) for x in open(args.hsps).read().split("\n")])

    logger.info(f"Loading the training pairs...")
    training_pairs = [tuple(x.rstrip("\n").split()) for x in open(args.training_pairs)]

    logger.info(f"Scoring the interactions...")
    scores = rsprint.score_interactions(proteins, hsps, training_pairs, kmer_size=args.kmer_size)
    predictions = Predictions(scores, proteins)

    logger.info(f"Saving the scores to a {args.output}...")
    predictions.save(args.output)

    logger.info(f"Done!")

def score_peptides(args):
    logger.info("Loading the protein sequences...")
    proteins = [(p.id, str(p.seq)) for p in Bio.SeqIO.parse(args.input, "fasta")]

    logger.info("Loading the peptide sequences...")
    peptides = [(p.id, str(p.seq)) for p in Bio.SeqIO.parse(args.peptides, "fasta")]

    logger.info(f"Loading the protein HSPs...")
    hsps = set([tuple(x.split()) for x in open(args.hsps).read().split("\n")])

    logger.info(f"Loading the training pairs...")
    training_pairs = [tuple(x.rstrip("\n").split()) for x in open(args.training_pairs)]

    logger.info(f"Extracting the peptide hsps...")
    peptide_hsps = rsprint.extract_peptide_hsps(proteins, peptides, t_smer=args.t_smer, t_hsp=args.t_hsp, kmer_size=args.kmer_size)
    hsps = hsps.union(peptide_hsps)

    logger.info(f"Scoring the interactions...")
    scores = rsprint.score_interactions(proteins, hsps, training_pairs, kmer_size=args.kmer_size)
    predictions = Predictions(scores, proteins)

    logger.info(f"Saving the scores to a {args.output}...")
    predictions.save(args.output)

    logger.info(f"Done!")