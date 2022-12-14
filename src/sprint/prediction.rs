use crate::sprint::hsp::HSP;
use crate::sprint::proteinset::ProteinSet;
use crate::sprint::scoring::score_hsp;
use ndarray::{Array2, Axis};
use rayon::prelude::*;
use std::cell::UnsafeCell;
use std::collections::HashSet;
use std::time::Instant;

pub struct PredictionMatrix {
    pub scores: UnsafeCell<Array2<f32>>,
}

impl PredictionMatrix {
    pub fn new(protein_set_size: usize) -> PredictionMatrix {
        PredictionMatrix {
            scores: UnsafeCell::new(Array2::zeros((protein_set_size, protein_set_size))),
        }
    }
}

unsafe impl Sync for PredictionMatrix {}
unsafe impl Send for PredictionMatrix {}

pub fn score_interactions(
    protein_set: &ProteinSet,
    hsps: &HashSet<HSP>,
    training_pairs: &Vec<(String, String)>,
    kmer_size: usize,
    process_rank: usize,
    world_size: usize,
    verbose: bool,
) -> Array2<f32> {
    let mapped_training_pairs: Vec<(usize, usize)> = training_pairs
        .iter()
        .filter(|pair| protein_set.contains(&pair.0) && protein_set.contains(&pair.1))
        .map(|pair| {
            (
                protein_set.get_protein_by_name(&pair.0).index(),
                protein_set.get_protein_by_name(&pair.1).index(),
            )
        })
        .collect();

    if verbose {
        println!("Process {} - Identifying interactors ", process_rank);
    }

    let mut interactors: HashSet<usize> = HashSet::new();
    mapped_training_pairs.iter().for_each(|pair| {
        interactors.insert(pair.0);
        interactors.insert(pair.1);
    });

    if verbose {
        println!(
            "Process {} - Identification of relevant pairs.",
            process_rank
        );
    }
    // Prepare the batch of training pairs to use in this process
    let mut training_pairs_to_process: Vec<(usize, usize)> = vec![];
    for i in 0..mapped_training_pairs.len() {
        if (i + process_rank) % world_size == 0 {
            training_pairs_to_process.push(mapped_training_pairs[i].clone());
        }
    }

    if verbose {
        println!("Process {} -: Building the HSP table", process_rank);
    }
    let hsp_table = build_hsp_table(hsps, &protein_set, &interactors, kmer_size);

    unsafe {
        if verbose {
            println!("Process {} - Initializing the score matrix.", process_rank);
        }
        let matrix = PredictionMatrix::new(protein_set.len());

        if verbose {
            println!("Process {} - Scoring the interactions...", process_rank);
        }
        let start = Instant::now();

        training_pairs_to_process
            .par_iter()
            .for_each(|pair| fill_matrix(pair, &hsp_table, kmer_size as f32, &matrix));

        // Divide the diagonal by two because it is double scored
        let matrix_ptr = matrix.scores.get();
        for i in 0..(*matrix_ptr).len_of(Axis(0)) {
            (*matrix_ptr)[[i, i]] /= 2.0;
        }

        if verbose {
            println!(
                "Process {} - Scored the interactions in {}...",
                process_rank,
                start.elapsed().as_secs()
            );
        }

        (*matrix.scores.get()).clone()
    }
}

pub unsafe fn initialize_score_matrix(matrix: &mut Vec<f32>, protein_set: &ProteinSet) {
    let sequence_set_size = protein_set.len();
    matrix.clear();
    for _ in 0..sequence_set_size * (sequence_set_size + 1) / 2 {
        matrix.push(0f32);
    }
}

/// Creates a table of HSPs used in scoring
///
/// It lists the HSPs in which a given protein is involved
/// The format of the table is a Vec<Vec<ScoredHSP>> as follows:
///
/// <--------------- HSPs (vector) --------------->
/// ^
/// |
/// |
/// Protein index
/// |
/// |
/// v
///
/// The HSPs are sorted in order of the partner protein's index (entry 0 of the scored HSP tuple)
/// Tuples are (partner_index, partner protein length, length of HSP, HSP score)
pub fn build_hsp_table(
    hsps: &HashSet<HSP>,
    protein_set: &ProteinSet,
    interactors: &HashSet<usize>,
    kmer_size: usize,
) -> Vec<Vec<(usize, f32, f32, f32)>> {
    let mut table: Vec<Vec<(usize, f32, f32, f32)>> = Vec::new();

    // Initialize the table
    for _ in 0..protein_set.len() {
        table.push(vec![]);
    }

    hsps.iter().for_each(|hsp| {
        let protein1 = protein_set.get_protein_by_id(hsp.location(0).index());
        let protein2 = protein_set.get_protein_by_id(hsp.location(1).index());
        let hsp_score = score_hsp(
            protein1,
            protein2,
            hsp.location(0).position(),
            hsp.location(1).position(),
            hsp.len(),
            kmer_size,
        );

        if interactors.contains(&hsp.location(0).index())
            || interactors.contains(&hsp.location(1).index())
        {
            table[hsp.location(0).index()].push((
                hsp.location(1).index(),
                protein2.len() as f32,
                hsp.len() as f32,
                hsp_score as f32,
            ));

            // If the HSP is not within the same protein, add the reciprocal HSP
            // in the other row of the table
            if hsp.location(0).index() != hsp.location(1).index() {
                table[hsp.location(1).index()].push((
                    hsp.location(0).index(),
                    protein1.len() as f32,
                    hsp.len() as f32,
                    hsp_score as f32,
                ));
            }
        }
    });

    for row in &mut table {
        row.sort_by_key(|partner| partner.0);
    }

    table
}

#[inline(always)]
pub fn get_1d_index(position1: usize, position2: usize) -> usize {
    let smallest = std::cmp::min(position1, position2);
    let largest = std::cmp::max(position1, position2);
    let row_start = largest * (largest + 1) / 2;
    return row_start + smallest;
}

/// Fill the score matrix using the similarity-to-interacting-pair principle
///
///               Interactor 1 --------------- Interactor 2
///                    |                            |
///               HSP1 |                       HSP2 |
///                    |                            |
///                 Query 1                      Query 2
///
pub unsafe fn fill_matrix(
    interacting_pair: &(usize, usize),
    hsps: &Vec<Vec<(usize, f32, f32, f32)>>,
    kmer_size: f32,
    prediction_matrix: &PredictionMatrix,
) {
    let matrix_ptr = prediction_matrix.scores.get();

    // If the training pair used for scoring involved the same protein (oligomerization),
    // then, only use the HSP pairs once!
    if interacting_pair.0 == interacting_pair.1 {
        for i in 0..hsps[interacting_pair.0].len() {
            for j in i..hsps[interacting_pair.1].len() {
                let hsp1 = &hsps[interacting_pair.0][i];
                let hsp2 = &hsps[interacting_pair.1][j];

                let term1 = hsp1.3 * (hsp2.2 - kmer_size + 1f32);
                let term2 = hsp2.3 * (hsp1.2 - kmer_size + 1f32);
                let contribution = (term1 + term2) / (hsp1.1 * hsp2.1);
                (*matrix_ptr)[[hsp1.0, hsp2.0]] += contribution;
                (*matrix_ptr)[[hsp2.0, hsp1.0]] += contribution;
            }
        }

    // We are not concerned with this if the interacting proteins are different.
    } else {
        let hsps1 = &hsps[interacting_pair.0];

        for hsp1 in hsps1 {
            let hsps2 = &hsps[interacting_pair.1];

            for hsp2 in hsps2 {
                let term1 = hsp1.3 * (hsp2.2 - kmer_size + 1f32);
                let term2 = hsp2.3 * (hsp1.2 - kmer_size + 1f32);
                let contribution = (term1 + term2) / (hsp1.1 * hsp2.1);
                (*matrix_ptr)[[hsp1.0, hsp2.0]] += contribution;
                (*matrix_ptr)[[hsp2.0, hsp1.0]] += contribution;
            }
        }
    }
}
