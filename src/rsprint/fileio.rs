use std::{io::Error, collections::HashSet};
use std::fs::read_to_string;

use rayon::prelude::*;
use bio::io::fasta::{Reader};

use super::protein::Protein;
use super::hsp::HSP;
use super::proteinset::ProteinSet;

/// Open a fasta file and returns a vector of Protein
pub fn load_fasta(filename: &str, training: bool) -> Result<Vec<Protein>, Error> {
    let reader = Reader::from_file(filename).unwrap();
    let mut proteins = vec![];

    for (i, protein) in reader.records().enumerate() {
        let record = protein?;
        let protein_id = String::from(record.id());
        let sequence_string: String = 
            std::str::from_utf8(record.seq()).unwrap().to_string();
        let new_protein: Protein = Protein::new(i, protein_id.clone(), sequence_string, training);
        proteins.push(new_protein);
    }

    Ok(proteins)
}

///// Saves a HashSet of HSPs to a file with the locations 
///// for a protein pair sorted by position
pub fn save_hsps(hsps: HashSet<HSP>,
    protein_set: &ProteinSet,
    filename: &str) -> Result<(), Error> {

    let output = hsps.into_iter()
    .map(|hsp| hsp.as_string(protein_set))
    .collect::<Vec<String>>()
    .join("\n");

    std::fs::write(filename, output)?;
    Ok(())
}

/// Loads the HSPs from a file into a HashSet of HSPs
pub fn load_hsps(filename: &str, protein_set: &ProteinSet) -> HashSet<HSP> {
    let file_contents = read_to_string(filename)
        .expect(&format!("The file {} could not be read.", filename));
    let mut hsps: HashSet<HSP> = HashSet::new();

    for line in file_contents.lines() {
        hsps.insert(HSP::from(line.to_owned(), protein_set));
    }

    hsps
}

/// Loads protein pairs
pub fn load_pairs(filename: &str) -> Vec<(String, String)> {
    let file_contents: String = read_to_string(filename)
        .expect(&format!("The file {} could not be read.", filename));

    let mut parsed: HashSet<(String, String)> = HashSet::new();
    let mut pairs = Vec::new();

    for (i, line) in file_contents.lines().enumerate() {
        let mut delimiter: &str = " "; //default
        if line.contains(" ") {
            delimiter = " ";
        } else if line.contains(",") {
            delimiter = ",";
        } else if line.contains("\t") {
            delimiter = "\t";
        } else {
            println!("WARNING: The pairs in a protein pair file must all be on their own
                   line and separated by a space, a comma or a tab. ({}, line {})", filename, i + 1);
        }
        let proteins: Vec<&str> = line.split(delimiter).collect();

        // Push them to the pairs vector, removing all whitespaces
        let pair = (proteins[0].split_whitespace().collect::<String>(), 
                                    proteins[1].split_whitespace().collect::<String>());
        let inverted_pair = (pair.1.clone(), pair.0.clone());

        // Only add new pairs, not existing pairs
        // The order of the proteins in the pair does not matter
        if !parsed.contains(&pair) {
            parsed.insert(pair.clone());
            parsed.insert(inverted_pair);
            pairs.push(pair);
        }
    }

    pairs
}

/// Saves scores
pub fn save_scores(scores: &Vec<f32>, 
                   protein_set: &ProteinSet,
                 filename: &str) -> std::io::Result<()> {

    // Serial
    //let mut file_contents = String::from("");
    //for i in 0..sequence_set.len() {
    //    let num_sequences_in_row = i + 1;
    //    for j in 0..num_sequences_in_row {
    //        let protein1 = sequence_set.get_protein_name(i);
    //        let protein2 = sequence_set.get_protein_name(j);
    //        file_contents.push_str(&format!("{} {} {}\n", protein1, protein2, scores[(i * (i + 1))/2 + j]));
    //    }
    //}
    let file_contents: String = (0..protein_set.len())
        .collect::<Vec<usize>>()
        .par_iter()
        .map(|&row_index| {
            let mut row_scores: Vec<String> = vec![];
            let num_sequences_in_row = row_index + 1;
            for j in 0..num_sequences_in_row {
                let protein1 = protein_set.get_protein_by_id(row_index).name();
                let protein2 = protein_set.get_protein_by_id(j).name();
                row_scores.push(format!("{} {} {}", protein1, protein2, scores[(row_index * (row_index + 1))/2 + j]));
            }
            row_scores.join("\n")
        })
    .collect::<Vec<String>>()
    .join("\n");

    //Parallel
    std::fs::write(filename, file_contents)?;
    Ok(())
}