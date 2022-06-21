use clap::{App, Arg};

use rsprint::rsprint::{fileio::{load_hsps, save_hsps}, proteinset::ProteinSet, processing::{process_hsps}};

fn main() {
    let matches = App::new("HSP processing")
    .author("Francois Charih, francois@charih.ca")
    .version("1.0.0")
    .about("Processes the HSPs as done in the SPRINT paper.")
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
    .arg(
        Arg::new("count_threshold")
        .short('t')
        .long("count_threshold")
        .required(false)
        .default_value("40")
        .value_name("INTEGER")
        .help("Frequency count above which an HSP is split.")
    )
    .get_matches();

    let input = matches.value_of("input")
        .expect("You did not provide an appropriate FASTA file with sequences.")
        .to_string();
    let hsp_file = matches.value_of("hsps")
        .expect("You did not provide an appropriate file where to dump the HSPs.")
        .to_string();
    let output = matches.value_of("output")
        .expect("You did not provide an appropriate file where to dump the HSPs.")
        .to_string();
    let kmer_size = matches.value_of("kmer_size")
        .unwrap()
        .parse::<usize>()
        .expect("You did not provide an integer value for kmer_size!");
    let t_count = matches.value_of("count_threshold")
        .unwrap()
        .parse::<u16>()
        .expect("You did not provide an integer value for count_threshold.");
        
    let protein_set = ProteinSet::from_file(&input).unwrap();
    let hsps = load_hsps(&hsp_file, &protein_set);
    let processed = process_hsps(&protein_set, hsps, kmer_size, t_count, true);
    save_hsps(processed, &protein_set, &output).unwrap();
}