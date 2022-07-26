use clap::Parser;

use rsprint::rsprint::proteinset::ProteinSet;
use rsprint::rsprint::extraction::extract_hsps;
use rsprint::rsprint::fileio::save_hsps;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
pub struct ExtractionArgs {
    #[clap(value_parser, short='i', long="input")]
    pub input_path: String,

    #[clap(value_parser, short='o', long="output")]
    pub output_path: String,

    #[clap(value_parser, short='k', long="kmer_size", default_value="20")]
    pub kmer_size: usize,

    #[clap(value_parser, short='s', long="t_sim", default_value="15")]
    pub t_sim: i16,

    #[clap(value_parser, short='t', long="t_hsp", default_value="35")]
    pub t_hsp: i16,
}

fn main() {
    let args = ExtractionArgs::parse();

    // Load the sequences
    let set = ProteinSet::from_file(&args.input_path).unwrap();

    // Extract the HSPs
    let hsps = extract_hsps(
        &set,
        args.kmer_size,
        args.t_sim,
        args.t_hsp,
        0,
        1,
        false,
        true,
        true
    );

    // Save the HSPs to a file
    save_hsps(hsps, &set, &args.output_path).unwrap();
}