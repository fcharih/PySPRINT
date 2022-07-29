from types import SimpleNamespace
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO

def main(args: SimpleNamespace):
    df = pd.read_csv(args.input)

    ## Find the number of peptides and the target protein length
    grouped_peptides = df.groupby("peptide")
    num_peptides = len(grouped_peptides.size())
    target_length = grouped_peptides.size()[0]

    ## Format the data into a matrix form with peptide along rows
    ## and position in target protein along columns
    pivoted = df.pivot(index="peptide", columns="target_position", values="contribution")
    hmap = pivoted.values

    if args.sequence:
        records = [rec for rec in SeqIO.parse(args.sequence, "fasta")]
        xticklabels = list(str(records[0].seq)) # There should only be 1 seq (target)
    else:
        xticklabels = list(range(1, target_length + 1))

    ## Setup the plot
    fig, ax = plt.subplots(figsize=(100, 30))
    im = ax.imshow(hmap)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=12)
    ax.set_xticks(range(0, target_length), labels=xticklabels)
    ax.set_yticks(range(0, num_peptides), labels=pivoted.index.values)
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
         rotation_mode="default")

    # Save or show plot
    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Plots the contribution of amino acids in the target protein to the interaction score.")
    parser.add_argument("-i", "--input", type=str, help="Path to the CSV file containing the contributions of the amino acids.")
    parser.add_argument("-s", "--sequence", type=str, help="Path to FASTA file containing the sequence of the targeted protein.")
    parser.add_argument("-o", "--output", type=str, required=False, help="Path to where the plot should be saved (if not provided, plot is displayed).")
    args = parser.parse_args()

    main(args)