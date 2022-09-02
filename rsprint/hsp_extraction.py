#!/usr/bin/env python3

"""Module containing the function to extract HSPs from a list of proteins

Copyright (c) 2022, François Charih
This source code is licensed under the MIT license found in the
LICENSE file in the root directory of this source tree.
"""

import Bio.SeqIO
import rsprint

def main(args):
    proteins = [(p.id, str(p.seq)) for p in Bio.SeqIO.parse(args.input, "fasta")]
    hsps = rsprint.extract_hsps(proteins)
    # TODO process hsps here
    with open(args.output, "w") as output_file:
        output_file.write("\n".join([f"{x[0]} {x[2]} {x[1]} {x[3]} {x[4]}" for x in hsps]))
