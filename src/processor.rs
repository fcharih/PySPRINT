use clap::Parser;
use std::path::Path;

use rsprint::rsprint::{fileio::{load_hsps, save_hsps}, proteinset::ProteinSet, processing::{process_hsps}};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
pub struct ProcessingArgs {
    #[clap(value_parser, short='i', long="sequences")]
    pub sequences_path: String,

    #[clap(value_parser, short='s', long="hsps")]
    pub hsps_path: String,

    #[clap(value_parser, short='o', long="output")]
    pub output_path: String,

    #[clap(value_parser, short='k', long="kmer_size", default_value="20")]
    pub kmer_size: usize,

    #[clap(value_parser, short='c', long="count_threshold", default_value="40")]
    pub count_threshold: u16,
}

fn main() {
    let args = ProcessingArgs::parse();

    // Ensure that the output has a .phsp extension (processed hsp)
    let output_path = Path::new(&args.output_path);
    let output_filepath = output_path.with_extension("phsp");

    // Load the sequences
    let protein_set = ProteinSet::from_file(&args.sequences_path).unwrap();

    // Load the unprocessed HSPs
    let hsps = load_hsps(&args.hsps_path, &protein_set);

    // Process the HSPs
    let processed = process_hsps(
        &protein_set, hsps, args.kmer_size, args.count_threshold, true);

    // Save the processed HSPs to a file
    save_hsps(processed, &protein_set, &output_filepath.to_str().unwrap()).unwrap();
}
