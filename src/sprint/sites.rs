use crate::sprint::hsp::HSP;
use crate::sprint::proteinset::ProteinSet;
use crate::sprint::scoring::score_hsp;
use rayon::prelude::*;
use std::cell::UnsafeCell;
use std::collections::{HashSet, HashMap};
use std::time::Instant;

pub struct ContributionsDict {
    pub dict: UnsafeCell<HashMap<usize, Vec<f32>>>
}

impl ContributionsDict {
    pub fn new(target_index: usize, protein_set: &ProteinSet) -> ContributionsDict {
        let target_length = protein_set.get_protein_by_id(target_index).len();
        let mut contributions_dict = HashMap::new();

        for p in protein_set.iter() {
            if p.is_new() {
                contributions_dict.insert(p.index(), vec![0f32; target_length]);
            }
        }

        ContributionsDict { dict: UnsafeCell::new(contributions_dict) }
    }
}

unsafe impl Sync for ContributionsDict {}
unsafe impl Send for ContributionsDict {}


pub fn compute_contributions(
    target_name: &String,
    protein_set: &ProteinSet,
    hsps: &HashSet<HSP>,
    training_pairs: &Vec<(String, String)>,
    kmer_size: usize,
    process_rank: usize,
    world_size: usize,
    verbose: bool,
) -> HashMap<usize, Vec<f32>> {
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

        let target_index = protein_set.get_protein_by_name(target_name).index();
        let dict = ContributionsDict::new(target_index, protein_set);

        if verbose {
            println!("Process {} - Scoring the interactions...", process_rank);
        }
        let start = Instant::now();

        training_pairs_to_process
            .par_iter()
            .for_each(|pair| fill_dict(target_index, pair, &hsp_table, kmer_size as f32, protein_set, &dict));

        if verbose {
            println!(
                "Process {} - Scored the interactions in {}...",
                process_rank,
                start.elapsed().as_secs()
            );
        }

        return (*dict.dict.get()).clone();
    }
}

pub unsafe fn initialize_score_matrix(matrix: &mut Vec<f32>, protein_set: &ProteinSet) {
    let sequence_set_size = protein_set.len();
    matrix.clear();
    for _ in 0..sequence_set_size * (sequence_set_size + 1) / 2 {
        matrix.push(0f32);
    }
}

pub fn build_hsp_table(
    hsps: &HashSet<HSP>,
    protein_set: &ProteinSet,
    interactors: &HashSet<usize>,
    kmer_size: usize,
) -> Vec<Vec<(usize, f32, f32, f32, usize, usize)>> {
    let mut table: Vec<Vec<(usize, f32, f32, f32, usize, usize)>> = Vec::new();

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
                hsp.location(0).position(),
                hsp.location(1).position()
            ));

            // If not within the same pro
            if hsp.location(0).index() != hsp.location(1).index() {
                table[hsp.location(1).index()].push((
                    hsp.location(0).index(),
                    protein1.len() as f32,
                    hsp.len() as f32,
                    hsp_score as f32,
                    hsp.location(1).position(),
                    hsp.location(0).position()
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

pub unsafe fn fill_dict(
    target_index: usize,
    interacting_pair: &(usize, usize),
    hsps: &Vec<Vec<(usize, f32, f32, f32, usize, usize)>>,
    kmer_size: f32,
    protein_set: &ProteinSet,
    contributions_dict: &ContributionsDict,
) {
    let dict_ptr = contributions_dict.dict.get();

    if interacting_pair.0 == interacting_pair.1 {
        for i in 0..hsps[interacting_pair.0].len() {
            for j in i..hsps[interacting_pair.1].len() {
                let hsp1 = &hsps[interacting_pair.0][i];
                let hsp2 = &hsps[interacting_pair.1][j];

                let term1 = hsp1.3 * (hsp2.2 - kmer_size + 1f32);
                let term2 = hsp2.3 * (hsp1.2 - kmer_size + 1f32);
                let contribution = (term1 + term2) / (hsp1.1 * hsp2.1); // TODO divide at the end

                let p1 = protein_set.get_protein_by_id(hsp1.0);
                let p2 = protein_set.get_protein_by_id(hsp2.0);
                let p1_is_target = p1.index() == target_index;

                if (p1.is_new() && p2.index() == target_index) || (p2.is_new() && p1.index() == target_index) {
                    if p1_is_target {
                        let vec = (*dict_ptr).get_mut(&p2.index()).unwrap();
                        for i in 0..hsp1.2 as usize {
                            vec[hsp1.5 + i] += contribution / hsp1.2; // Contribution is distributed over similarity region
                        }
                    } else {
                        let vec = (*dict_ptr).get_mut(&p1.index()).unwrap();
                        for i in 0..hsp2.2 as usize {
                            vec[hsp2.5 + i] += contribution / hsp2.2; 
                        }
                    }
                }
            }
        }
    } else {
        let hsps1 = &hsps[interacting_pair.0];

        for hsp1 in hsps1 {
            let hsps2 = &hsps[interacting_pair.1];

            for hsp2 in hsps2 {
                let term1 = hsp1.3 * (hsp2.2 - kmer_size + 1f32);
                let term2 = hsp2.3 * (hsp1.2 - kmer_size + 1f32);
                let contribution = (term1 + term2) / (hsp1.1 * hsp2.1); // TODO divide at the end


                let p1 = protein_set.get_protein_by_id(hsp1.0);
                let p2 = protein_set.get_protein_by_id(hsp2.0);
                let p1_is_target = p1.index() == target_index;

                if (p1.is_new() && p2.index() == target_index) || (p2.is_new() && p1.index() == target_index) {
                    if p1_is_target {
                        let vec = (*dict_ptr).get_mut(&p2.index()).unwrap();
                        for i in 0..hsp1.2 as usize {
                            vec[hsp1.5 + i] += contribution / hsp1.2;
                        }
                    } else {
                        let vec = (*dict_ptr).get_mut(&p1.index()).unwrap();
                        for i in 0..hsp2.2 as usize {
                            vec[hsp2.5 + i] += contribution / hsp2.2;
                        }
                    }
                }
            }
        }
    }
}

