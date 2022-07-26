use clap::Parser;

use rsprint::rsprint::{fileio::{load_hsps, load_pairs, save_scores}, proteinset::ProteinSet, prediction::{score_interactions}};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
pub struct PredictionArgs {
    #[clap(value_parser, short='i', long="sequences")]
    pub sequences_path: String,

    #[clap(value_parser, short='s', long="hsps")]
    pub hsps_path: String,

    #[clap(value_parser, long="training_pairs")]
    pub training_pairs_path: String,

    #[clap(value_parser, long="output")]
    pub output_path: String,

    #[clap(value_parser, long="kmer_size", default_value="20")]
    pub kmer_size: usize,
}
fn main() {

    let args = PredictionArgs::parse();

    // Load the sequences
    let protein_set = ProteinSet::from_file(&args.sequences_path).unwrap();

    // Load the processed HSPs
    let hsps = load_hsps(&args.hsps_path, &protein_set);

    // Load the training pairs used for scoring
    let training_pairs = load_pairs(&args.training_pairs_path);

    // Score the interactions
    let scores = score_interactions(
        &protein_set, &hsps, &training_pairs, args.kmer_size, 0, 1, true);

    // Save the scores to a file
    save_scores(&scores, &protein_set, &args.output_path).unwrap();
}