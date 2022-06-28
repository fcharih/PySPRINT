//use clap::{App, Arg};
//
//use rsprint::rsprint::{
//    fileio::{load_hsps, load_pairs, save_scores},
//    protein::Protein,
//    proteinset::ProteinSet,
//    //pymodules::SPRINT,
//};
//
//fn main() {
//    let matches = App::new("SPRINT peptide scorer")
//    .author("Francois Charih, francois@charih.ca")
//    .version("1.0.0")
//    .about("Scores the interactions involving the peptide using the training pairs and HSPs provided.")
//    .arg(
//        Arg::new("input")
//        .short('i')
//        .long("input")
//        .required(true)
//        .value_name("FILE")
//        .help("Path to the FASTA file containing the sequences from which HSPs should be extracted.")
//    )
//        .arg(
//        Arg::new("peptides")
//        .short('p')
//        .long("peptides")
//        .required(true)
//        .value_name("FILE")
//        .help("Path to the FASTA file containing the sequences of the peptides to score.")
//    )
//    .arg(
//        Arg::new("hsps")
//        .short('s')
//        .long("hsps")
//        .required(true)
//        .value_name("FILE")
//        .help("Path to the file containing the HSPs to process.")
//    )
//    .arg(
//        Arg::new("training_pairs")
//        .short('t')
//        .long("training_pairs")
//        .required(true)
//        .value_name("FILE")
//        .help("Path to the file containing the training pairs.")
//        )
//    .arg(
//        Arg::new("output")
//        .short('o')
//        .long("output")
//        .required(true)
//        .value_name("FILE")
//        .help("Path to the output file where HSPs should be saved.")
//    )
//    .arg(
//        Arg::new("kmer_size")
//        .short('k')
//        .long("kmer_size")
//        .required(false)
//        .default_value("20")
//        .value_name("INTEGER")
//        .help("Minimum length of an hsp.")
//    )
//    .arg(
//        Arg::new("smer_threshold")
//        .short('e')
//        .long("smer_threshold")
//        .required(false)
//        .default_value("15")
//        .value_name("INTEGER")
//        .help("Score above which two s-mers are considered similar.")
//    )
//    .arg(
//        Arg::new("hsp_threshold")
//        .short('t')
//        .long("hsp_threshold")
//        .required(false)
//        .default_value("35")
//        .value_name("INTEGER")
//        .help("Score above which two regions are considered HSPs.")
//    )
//
//    .get_matches();
//
//    let input = matches
//        .value_of("input")
//        .expect("You did not provide an appropriate FASTA file with sequences.")
//        .to_string();
//    let peptide_file = matches
//        .value_of("peptides")
//        .expect("You did not provide an appropriate FASTA file with sequences.")
//        .to_string();
//    let hsp_file = matches
//        .value_of("hsps")
//        .expect("You did not provide an appropriate file where to dump the HSPs.")
//        .to_string();
//    let training_pairs_file = matches
//        .value_of("training_pairs")
//        .expect("You did not provide an appropriate file with training pairs to use.")
//        .to_string();
//    let output = matches
//        .value_of("output")
//        .expect("You did not provide an appropriate file where to dump the HSPs.")
//        .to_string();
//    let kmer_size = matches
//        .value_of("kmer_size")
//        .unwrap()
//        .parse::<usize>()
//        .expect("You did not provide an integer value for kmer_size!");
//    let t_smer = matches
//        .value_of("smer_threshold")
//        .unwrap()
//        .parse::<i16>()
//        .expect("You did not provide an integer value for smer_threshold!");
//    let t_hsp = matches
//        .value_of("hsp_threshold")
//        .unwrap()
//        .parse::<i16>()
//        .expect("You did not provide an integer value for hsp_threshold!");
//
//    // Training proteins
//    let mut protein_set = ProteinSet::from_file(&input).unwrap();
//
//    // Peptides
//    let peptides_set = ProteinSet::from_file(&peptide_file).unwrap();
//
//    // Merged
//    let peptides_as_proteins: Vec<Protein> = peptides_set.iter().cloned().collect();
//    protein_set.add_new(peptides_as_proteins);
//
//    // TODO This is really hacky but meh...
//    let scorer = SPRINT::new(0, 1, input, hsp_file, training_pairs_file);
//    let peptides: Vec<(String, String)> =
//        peptides_set.iter().map(|p| (p.name(), p.seq())).collect();
//    let peptides_copy: Vec<(String, String)> = peptides.iter().cloned().collect();
//
//    let peptide_hsps = scorer
//        .extract_peptide_hsps(peptides, kmer_size, t_smer, t_hsp)
//        .unwrap();
//    let scores = scorer
//        .score_interactions(peptides_copy, peptide_hsps, kmer_size)
//        .unwrap();
//
//    save_scores(&scores, &protein_set, &output).unwrap();
//}
