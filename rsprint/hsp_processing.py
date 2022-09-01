#!/usr/bin/env python3

"""Module containing the function to process a set of HSPs

Copyright (c) 2022, François Charih
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""
import Bio.SeqIO
import rsprint

def main(args):
    print("processing hsps")
    proteins = [(p.id, str(p.seq)) for p in Bio.SeqIO.parse(args.input, "fasta")]
    # TODO

