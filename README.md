# rSPRINT

Author: François Charih <francois@charih.ca>

The SPRINT algorithm was designed by Li and Ilie
([paper](https://doi.org/10.1186/s12859-017-1871-x), [repo](https://github.com/lucian-ilie/SPRINT)). This repo contains a faster reimplementation in Rust for usage with Python.

## Installation

This package was only tested on macOS and Linux (Ubuntu).

**OpenMPI must be installed on your system**. This can be done
with [Homebrew (macOS)](https://formulae.brew.sh/formula/open-mpi)
or your package manager on Linux platforms.

### Building from source

If you have a recent version of [Rust](https://www.rust-lang.org/learn/get-started) on
your system, building from scratch is easy:

```
$ pip install git+https://github.com/fcharih/rSPRINT.git
```

**Boom! Done!**

### Installation with pip

TODO
## Using rSPRINT

I suggest using default settings unless you understand how the SPRINT algorithm
works.

### Extracting HSPs

The SPRINT algorithm requires one to first extract the regions of similarity
or HSPs (high-scoring segment pairs) between all proteins to be scored
and the proteins in the training set (protein pairs known to interact).
They are tuple that describe the location of the similarity region in *Protein A*
and *Protein B*: *(ID Prot A, ID Prot B, Position in A, Position in B, Length)*.

**Note:**
The proteins from which we should extract hsps are:
1. Proteins that we wish to compute interactions for (all proteins in all-to-all runs).
2. Proteins that are in the training data.

Normally, this means the human proteome.

This can be achieved in CLI mode:

```
$ rsprint extract_hsps -h
```

Arguments listed when running with `-h`.

The HSPs can also be extracted with a Python call:

```python
from rsprint import extract_hsps

# Proteins from which HSPs should be extracted (DO NOT INCLUDE PEPTIDE HERE)
proteins = [("protein1", "SEQUENCE"), ...]

# Extract the HSPs (result is a set of tuples representing the HSPs)
hsps = extract_hsps(proteins, t_smer=15, t_hsp=35, kmer_size=20)

# Use or save to a file
# 
# in this e.g. saved to plain text file where each line has the format
# protein 1,protein 2,position a, position b,length
with open("my_file", "w") as output_file:
    file_content = "\n".join([f"{h[0]},{h[1]},{h[2]},{h[3]},{h[4]}" for h in hsps])
    output_file.save(file_content)
```

The HSPs can be used or saved to a file (*e.g.* in plain text format or
in binary format with `pickle`).

### Processing HSPs

In order to account for the occurence of regions that occur frequently in
proteins, but that account little for actual protein interactions, (*e.g.*
signal peptides) the SPRINT algorithm alter HSPs that come from such regions.

```
$ rsprint process_hsps -h
```

Arguments listed when running with `-h`.

The HSPs can also be processed with a Python call:

```python
from rsprint import extract_hsps, process_hsps

# Proteins from which HSPs should be extracted (DO NOT INCLUDE PEPTIDE HERE)
proteins = [("protein1", "SEQUENCE"), ...]

# Extract the HSPs (result is a set of tuples representing the HSPs)
#
# Note: in this example, we extract them, but they could also be
#       loaded from a file.
hsps = extract_hsps(proteins, t_smer=15, t_hsp=35, kmer_size=20)

# Process the HSPs
processed_hsps = process_hsps(hsps)

# Save to a file or use
...
```

**HSPs extracted from peptide binders SHOULD NOT be processed.**

### Scoring all interactions

TODO

### Scoring new proteins (peptides)

This is what you want to use to score peptide binders! Note that you must
have already extracted and processed the HSPs from proteins that you want
to score the peptides against AND from proteins in the training set
(training/interacting pairs).

```python
from rsprint import score_peptides, Predictions

proteins = [("Protein 1", "SEQUENCE"), ...] # tuples
peptides_to_score = [("Peptide 1", "SEQUENCE"), ...] # tuples
hsps = ... # extracted or retrieved from file
training_pairs = [("Protein 1", "Protein 5"), ("Protein 3", "Protein 10"), ...] # tuples

prediction_matrix = score_peptides(proteins, peptides, hsps, training_pairs, kmer_size=20) # numpy array`

# Convert the matrix to a Predictions object
protein_names = [p[0] for p in proteins]
peptide_names = [p[1] for p in peptides]
predictions = Predictions(prediction_matrix, protein_names, peptide_names)

# Now you can easily retrieve scores for your peptides
scores = predictions.get_scores("Peptide 1") # returns a dict of scores for all proteins/peptides that interact with Peptide 1

# You can get a specific score
score = predictions.get_score("Peptide 1", "Protein 10") # returns a float (score)

# You can save the predictions to a file
predictions.save("some_file.mat")

# You can load the predictions from a file
Predictions.from_file("some_file.mat")
```