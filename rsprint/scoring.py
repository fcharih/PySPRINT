#!/usr/bin/env python3
"""Module containing the function to score interactions

Copyright (c) 2022, François Charih
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

import Bio.SeqIO
import rsprint

def score_all_to_all(args):
    print("Scoring all-to-all")
    # scores = rsprint.score_interactions()

def score_peptides(args):
    print("Scoring peptides")
    # scores = rsprint.score_interactions()
    #rsprint.score_peptides
    #proteins = [(p.id, str(p.seq)) for p in Bio.SeqIO.parse(args.input, "fasta")]
    #hsps = rsprint.extract_hsps(proteins)
    #with open(args.output, "w") as output_file:
    #    output_file.write("\n".join([f"{x[0]} {x[1]} {x[2]} {x[3]} {x[4]}" for x in hsps]))

