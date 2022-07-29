use clap::Parser;

use rsprint::rsprint::{
    proteinset::ProteinSet,
    fileio::{load_hsps, load_pairs, save_peptide_scores},
    extraction::extract_hsps,
    prediction::score_interactions
};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
pub struct PeptideScoringArgs {
    #[clap(value_parser, short='i', long="sequences")]
    pub sequences_path: String,

    #[clap(value_parser, short='p', long="peptides")]
    pub peptides_path: String,

    #[clap(value_parser, short='s', long="hsps")]
    pub hsps_path: String,

    #[clap(value_parser, short='r', long="training_pairs")]
    pub training_pairs_path: String,

    #[clap(value_parser, short='o', long="output")]
    pub output_path: String,

    #[clap(value_parser, long="t_sim", default_value="15")]
    pub t_sim: i16,

    #[clap(value_parser, long="t_hsp", default_value="35")]
    pub t_hsp: i16,

    #[clap(value_parser, short='k', long="kmer_size", default_value="20")]
    pub kmer_size: usize,
}

pub fn main() {

    let args = PeptideScoringArgs::parse();

    // Training proteins
    let mut protein_set = ProteinSet::from_file(&args.sequences_path).unwrap();

    // Add the peptides
    protein_set.add_from_file(&args.peptides_path, true);

    // Load training HSPs
    let mut hsps = load_hsps(&args.hsps_path, &protein_set);

    // Load the training pairs
    let training_pairs = load_pairs(&args.training_pairs_path);

    // Compute the peptide HSPs
    let peptide_hsps = extract_hsps(&protein_set, args.kmer_size, args.t_sim,
        args.t_hsp, 0, 1, true, true, false);

    // Merge the HSPs
    hsps.extend(peptide_hsps);

    // Score the interactions
    let score_matrix = score_interactions(&protein_set, &hsps,
        &training_pairs, args.kmer_size, 0, 1, false);

    // Save the scores
    save_peptide_scores(&score_matrix, &protein_set, &args.output_path).unwrap();
}