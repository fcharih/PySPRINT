#!/usr/bin/env python3

"""Extracts High-Scoring Segment Pairs from proteins in a FASTA file.

Copyright (c) 2022, François Charih
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""
import rsprint
import argparse
import mpi4py

from .hsp_extraction import main as extract_hsps
from .hsp_processing import main as process_hsps
from .scoring import score_all_to_all, score_peptides

def main():
    parser = argparse.ArgumentParser(description=("rSPRINT PPI suite."))
    subparsers = parser.add_subparsers(dest='command', help='Command to run', required=True)

    # HSP extraction
    hsp_extraction = subparsers.add_parser("extract_hsps", help="Extract HSPs from a set of protein sequences.")
    hsp_extraction.add_argument("-i", "--input", type=str, required=True,
        help=("FASTA-formatted file containing the protein sequences."))
    hsp_extraction.add_argument("-o", "--output", type=str, required=True,
        help=("Path to the raw HSP file."))
    hsp_extraction.add_argument("-s", "--t_sim", type=int, required=False, default=15,
        help=("Threshold at which two s-mers are considered similar."))
    hsp_extraction.add_argument("-p", "--t_hsp", type=int, required=False, default=35,
        help=("Threshold at which two HSPs are considered similar."))
    hsp_extraction.add_argument("-k", "--kmer_size", type=int, required=False, default=20,
        help=("Minimum length of an HSP."))
    hsp_extraction.add_argument("-c", "--t_count", type=int, required=False, default=40,
        help=("Threshold on the max number of times a residue can be involved in a HSP before being removed."))
    hsp_extraction.set_defaults(func=extract_hsps)

    # HSP processing
    hsp_processing = subparsers.add_parser("process_hsps", help="Extract HSPs from a set of peptide sequences.")
    hsp_processing.add_argument("-i", "--input", type=str, required=True,
        help=("FASTA-formatted file containing the protein sequences."))
    hsp_processing.add_argument("-o", "--output", type=str, required=True,
        help=("Path to the output file."))
    hsp_processing.add_argument("-s", "--hsps", type=str, required=True,
        help=("Path to the raw HSP file."))
    hsp_processing.add_argument("-c", "--count_threshold", type=int, required=False, default=40,
        help=("Threshold at which two HSPs are considered similar."))
    hsp_processing.add_argument("-k", "--kmer_size", type=int, required=False, default=20,
        help=("Minimum length of an HSP."))
    hsp_processing.set_defaults(func=process_hsps)

    # All-to-all scoring
    all_to_all = subparsers.add_parser("extract_peptide_hsps", help="Score all PPIs.")
    all_to_all.add_argument("-i", "--input", type=str, required=True,
        help=("FASTA-formatted file containing the protein sequences."))
    all_to_all.add_argument("-o", "--output", type=str, required=True,
        help=("Path to the output file."))
    all_to_all.add_argument("-s", "--hsps", type=str, required=True,
        help=("Path to the processed HSP file."))
    all_to_all.add_argument("-t", "--training_pairs", type=str, required=True,
        help=("Path to the training pairs (one per line separated by a space or tab)."))
    all_to_all.add_argument("-k", "--kmer_size", type=int, required=False, default=20,
        help=("Minimum length of an HSP."))
    all_to_all.set_defaults(func=score_all_to_all)

    # Peptide scoring
    peptide_scoring = subparsers.add_parser("extract_peptide_hsps", help="Score peptides against a target.")
    peptide_scoring.add_argument("-i", "--input", type=str, required=True,
        help=("FASTA-formatted file containing the protein sequences."))
    peptide_scoring.add_argument("-p", "--peptides", type=str, required=True,
        help=("FASTA-formatted file containing the peptide sequences."))
    peptide_scoring.add_argument("-o", "--output", type=str, required=True,
        help=("Path to the output file."))
    peptide_scoring.add_argument("-s", "--hsps", type=str, required=True,
        help=("Path to the processed HSP file."))
    peptide_scoring.add_argument("-t", "--training_pairs", type=str, required=True,
        help=("Path to the training pairs (one per line separated by a space or tab)."))
    peptide_scoring.add_argument("-k", "--kmer_size", type=int, required=False, default=20,
        help=("Minimum length of an HSP."))
    peptide_scoring.set_defaults(func=score_peptides)

    args = parser.parse_args()
    args.func(args)