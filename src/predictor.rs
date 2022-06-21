use clap::{App, Arg};

use rsprint::rsprint::{fileio::{load_hsps, load_pairs, save_scores}, proteinset::ProteinSet, prediction::{score_interactions}};

fn main() {
    let matches = App::new("SPRINT predictor")
    .author("Francois Charih, francois@charih.ca")
    .version("1.0.0")
    .about("Scores the interaction using the training pairs and HSPs provided.")
    .arg(
        Arg::new("input")
        .short('i')
        .long("input")
        .required(true)
        .value_name("FILE")
        .help("Path to the FASTA file containing the sequences from which HSPs should be extracted.")
    )
    .arg(
        Arg::new("hsps")
        .short('s')
        .long("hsps")
        .required(true)
        .value_name("FILE")
        .help("Path to the file containing the HSPs to process.")
    )
    .arg(
        Arg::new("training_pairs")
        .short('t')
        .long("training_pairs")
        .required(true)
        .value_name("FILE")
        .help("Path to the file containing the training pairs.")
        )
    .arg(
        Arg::new("output")
        .short('o')
        .long("output")
        .required(true)
        .value_name("FILE")
        .help("Path to the output file where HSPs should be saved.")
    )
    .arg(
        Arg::new("kmer_size")
        .short('k')
        .long("kmer_size")
        .required(false)
        .default_value("20")
        .value_name("INTEGER")
        .help("Minimum length of an hsp.")
    )
    .get_matches();

    let input = matches.value_of("input")
        .expect("You did not provide an appropriate FASTA file with sequences.")
        .to_string();
    let hsp_file = matches.value_of("hsps")
        .expect("You did not provide an appropriate file where to dump the HSPs.")
        .to_string();
    let training_pairs_file = matches.value_of("training_pairs")
        .expect("You did not provide an appropriate file with training pairs to use.")
        .to_string();
    let output = matches.value_of("output")
        .expect("You did not provide an appropriate file where to dump the HSPs.")
        .to_string();
    let kmer_size = matches.value_of("kmer_size")
        .unwrap()
        .parse::<usize>()
        .expect("You did not provide an integer value for kmer_size!");
        
    let protein_set = ProteinSet::from_file(&input).unwrap();
    let hsps = load_hsps(&hsp_file, &protein_set);
    let training_pairs = load_pairs(&training_pairs_file);
    let scores = score_interactions(&protein_set, hsps, training_pairs, kmer_size, 0, 1, true);
    save_scores(&scores, &protein_set, &output).unwrap();
}