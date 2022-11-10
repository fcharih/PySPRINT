# rSPRINT

Author: François Charih <francois@charih.ca>

The SPRINT algorithm was designed by Li and Ilie
([paper](https://doi.org/10.1186/s12859-017-1871-x), [repo](https://github.com/lucian-ilie/SPRINT)). This repo contains a faster reimplementation in Rust for usage with Python.

## Installation

This package was only tested on macOS and Linux (Ubuntu).

**OpenMPI must be installed on your system**. This can be done
with Homebrew (macOS) or your package manager on Linux platforms.

### Building from source

If you have a recent version of [Rust](https://www.rust-lang.org/learn/get-started) on
your system, building from scratch is easy:

```
$ pip install git+https://github.com/fcharih/rSPRINT.git
```

**Boom! Done!**

### Installation with pip

This is in progress...

## Using rSPRINT

### Extracting HSPs

The SPRINT algorithm requires one to first extract the regions of similarity
or HSPs (high-scoring segment pairs) between all proteins to be scored
and the proteins in the training set (protein pairs known to interact).
They are tuple that describe the location of the similarity region in *Protein A*
and *Protein B*: *(ID Prot A, ID Prot B, Position in A, Position in B, Length)*.

This can be achieved in CLI mode:

```
$ rsprint extract_hsps <params>
```

See the help file for the parameters.

The HSPs can also be extracted with a Python call:

```python
from rsprint import extract_hsps

# Proteins from which HSPs should be extracted (DO NOT INCLUDE PEPTIDE HERE)
proteins = [("protein1", "SEQUENCE"), ...]

# Extract the HSPs (result is a set of tuples representing the HSPs)
hsps = extract_hsps(proteins, t_smer=15, t_hsp=35, kmer_size=20)
```

The HSPs can be used or saved to a file (*e.g.* in plain text format or
in binary format with `pickle`).

### Processing HSPs

### Scoring all interactions

### Scoring new proteins (peptides)