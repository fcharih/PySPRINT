use crate::rsprint::protein::Protein;
use crate::rsprint::constants::{PAM120};

pub fn score_sequences(protein1: &Protein, protein2: &Protein, start1: usize, start2: usize, length: usize) -> i16 {
    let mut score = 0i16;
    for i in 0..length {
        score += PAM120[protein1.residue(i + start1)][protein2.residue(i + start2)];
    }
    score
}

pub fn score_position(protein1: &Protein, protein2: &Protein, position1: usize, position2: usize) -> i16 {
    PAM120[protein1.residue(position1)][protein2.residue(position2)]
}

pub fn score_hsp(protein1: &Protein, protein2: &Protein, start1: usize, start2: usize, length: usize, kmer_size: usize) -> i16 {
    let mut score: i16 = 0;
    for i in 0..length - kmer_size + 1 {
        let mut inner_score: i16 = 0;
        for j in 0..kmer_size {
            inner_score += score_position(protein1, protein2, i + j + start1, i + j + start2);
        }
        score += inner_score
    }
    score
}

