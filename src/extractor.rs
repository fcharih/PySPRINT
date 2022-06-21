use clap::{App, Arg};

use rsprint::rsprint::proteinset::ProteinSet;
use rsprint::rsprint::extraction::extract_hsps;
use rsprint::rsprint::fileio::save_hsps;

fn main() {

    let matches = App::new("HSP extractor")
    .author("Francois Charih, francois@charih.ca")
    .version("1.0.0")
    .about("Extracts the HSPs from ")
    .arg(
        Arg::new("input")
        .short('i')
        .long("input")
        .required(true)
        .value_name("FILE")
        .help("Path to the FASTA file containing the sequences from which HSPs should be extracted.")
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
        .help("Minimum length of an HSP.")
    )
    .arg(
        Arg::new("smer_threshold")
        .short('s')
        .long("smer_threshold")
        .required(false)
        .default_value("15")
        .value_name("INTEGER")
        .help("Score above which two s-mers are considered similar.")
    )
    .arg(
        Arg::new("hsp_threshold")
        .short('t')
        .long("hsp_threshold")
        .required(false)
        .default_value("35")
        .value_name("INTEGER")
        .help("Score above which two regions are considered HSPs.")
    )
    .get_matches();

    let input = matches.value_of("input")
        .expect("You did not provide an appropriate FASTA file with sequences.")
        .to_string();
    let output = matches.value_of("output")
        .expect("You did not provide an appropriate file where to dump the HSPs.")
        .to_string();
    let kmer_size = matches.value_of("kmer_size")
        .unwrap()
        .parse::<usize>()
        .expect("You did not provide an integer value for kmer_size!");
    let t_smer = matches.value_of("smer_threshold")
        .unwrap()
        .parse::<i16>()
        .expect("You did not provide an integer value for smer_threshold!");
    let t_hsp = matches.value_of("hsp_threshold")
        .unwrap()
        .parse::<i16>()
        .expect("You did not provide an integer value for hsp_threshold!");

    let set = ProteinSet::from_file(&input).unwrap();
    let hsps = extract_hsps(
        &set,
        kmer_size,
        t_smer,
        t_hsp,
        0,
        1,
        false,
        true,
        true
    );
    save_hsps(hsps, &set, &output).unwrap();
}