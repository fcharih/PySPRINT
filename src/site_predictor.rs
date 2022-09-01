use clap::Parser;

use rsprint::rsprint::{fileio::{load_fasta, load_hsps, load_pairs, save_contributions}, proteinset::ProteinSet, extraction::{extract_hsps}, sites::{compute_contributions}};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
pub struct SitePredictionArgs {
    #[clap(value_parser, short='i', long="sequences")]
    pub sequences_path: String,

    #[clap(value_parser, short='p', long="peptides")]
    pub peptides_path: String,

    #[clap(value_parser, short='s', long="hsps")]
    pub hsps_path: String,

    #[clap(value_parser, short='t', long="target_name")]
    pub target_name: String,

    #[clap(value_parser, short='r', long="training_pairs")]
    pub training_pairs_path: String,

    #[clap(value_parser, long="t_sim", default_value="15")]
    pub t_sim: i16,

    #[clap(value_parser, long="t_hsp", default_value="35")]
    pub t_hsp: i16,

    #[clap(value_parser, short='o', long="output")]
    pub output_path: String,

    #[clap(value_parser, long="kmer_size", default_value="20")]
    pub kmer_size: usize,
}

fn main() {

    let args = SitePredictionArgs::parse();

    // Load the sequences
    let mut protein_set = ProteinSet::from_file(&args.sequences_path).unwrap();

    // Load the peptides and add to the protein set
    let peptides = load_fasta(&args.peptides_path, true).unwrap();
    protein_set.add_new(peptides, true);

    // Load the processed HSPs
    let mut hsps = load_hsps(&args.hsps_path, &protein_set);

    // Add the HSPs from the peptide sequences to these
    let new_hsps = extract_hsps(
        &protein_set,
        args.kmer_size,
        args.t_sim,
        args.t_hsp,
        0,
        1,
        true, // Only new HSPs (those involving peptides)
        true,
        true
    );

    hsps.extend(new_hsps);

    // Load the training pairs used for scoring
    let training_pairs = load_pairs(&args.training_pairs_path);

    // Compute the contributions of residues within the target to the interaction score
    // for the peptides of interest (new)
    let contributions = compute_contributions(
        &args.target_name, &protein_set, &hsps, &training_pairs, args.kmer_size, 0, 1, true);

    // Save the scores to a file
    save_contributions(&contributions, &protein_set, &args.output_path).unwrap();
}