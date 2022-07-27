use std::collections::{HashMap, HashSet};
use std::sync::Mutex;
use rayon::prelude::*;
use std::time::Instant;
use crate::rsprint::proteinset::ProteinSet;
use crate::rsprint::scoring::{score_sequences, score_position};
use crate::rsprint::seed::Seed; 
use crate::rsprint::smer::{Smer, SmerCollection}; 
use crate::rsprint::protein::Protein; 
use crate::rsprint::hsp::HSP;
use crate::rsprint::location::Location;
use crate::rsprint::constants::{RESIDUE_CODES, SEEDS};
use crate::rsprint::similarity::{compute_similar_smers};

/// Given a set of sequences (ProteinSet struct) and a Seed, extracts
/// all the s-mers from the sequences
pub fn extract_all_smers(sequences: &ProteinSet, seed: &Seed) -> Vec<SmerCollection> {

    // Extract the smers
    let smers: Vec<Smer>;
    smers = sequences.iter()
        .map(|prot| extract_smers(&prot, &seed))
        .flatten()
        .collect();

    // Group the smers by smer value into a hashmap that maps the
    // u64 value of the smer to its locations within the sequence_set
    let mut grouped_smers: HashMap<u64, HashSet<Location>>;
    grouped_smers = HashMap::new();
    smers.iter().for_each(|smer| {
       grouped_smers.entry(smer.value())
           .or_insert(HashSet::new())
           .insert(smer.location().clone());
    });

    // Create a vector of Smers (sorted by value)
    let mut smers_vec: Vec<SmerCollection> = grouped_smers.into_iter()
        .map(|(value, locations)| SmerCollection::new(value, locations))
        .collect();
    smers_vec.sort();

    smers_vec
}

/// Given a protein and a seed, extracts all the s-mers from that protein
pub fn extract_smers(protein: &Protein, seed: &Seed) -> Vec<Smer> {
    let mut smers = vec![];
    let seed_length = seed.len();
    let protein_length = protein.len();
    let num_smers = protein_length - seed_length + 1;

    if protein_length < seed_length {
        return smers;
    }

    for position in 0..num_smers {
        let start = position as usize;
        let end = position + seed.len() as usize - 1;
        let subsequence: &str = &protein.sub(start, end);
        let smer_value: u64 = compute_smer(subsequence, &seed);
        let location = Location::new(protein.index(), position);
        smers.push(Smer::new(smer_value, location));
    }
    smers
}

/// Computes the u64 value of an smer
pub fn compute_smer(sequence: &str, seed: &Seed) -> u64 {
    let mut u64_sequence: u64 = 0;

    for amino_acid in sequence.chars() {
        u64_sequence = u64_sequence << 5;
        u64_sequence = u64_sequence | *RESIDUE_CODES.get(&amino_acid).unwrap() as u64;
    }

    u64_sequence & seed.value()
}

pub fn extract_hsps(
    protein_set: &ProteinSet,
    kmer_size: usize,
    t_sim: i16,
    t_hit: i16,
    process_rank: usize,
    world_size: usize,
    new_only: bool,
    trivial_hsps: bool,
    verbose: bool
) -> HashSet<HSP> {

    let hsps: Mutex<HashSet<HSP>> = Mutex::new(HashSet::new());

    for (i, seed_string) in SEEDS.iter().enumerate() {

        // Create the seed struct
        let seed = Seed::new(seed_string);

        // Extract the seed's smers
        let smers = extract_all_smers(&protein_set, &seed);
        let smer_map: HashMap<u64, usize> = smers.iter().enumerate().map(|(i, smer)| (smer.value(), i)).collect();

        // Identify s-mers for this rank
        let mut smers_to_process = vec![];
        for i in 0..smers.len() {
            if (i + process_rank) % world_size == 0 {
                smers_to_process.push(i);
            }
        }

        // Compute the HSPs
        if verbose {
            println!("Process {} - Working with seed {}/{} ({}).", process_rank, i + 1, SEEDS.len(), seed_string);
        }
        let start = Instant::now();
        
        smers_to_process
            .into_par_iter()
            .for_each(|index| unsafe {
                // Compute hsps
                let results = compute_hsps_for_smer(index, &smers, &smer_map, &protein_set, &seed, kmer_size, t_sim, t_hit, new_only);

                // Return results
                let mut guard = hsps.lock().unwrap();
                (*guard).extend(results);
            });

        if verbose {
            println!("Process {} - Completed table {} in {}s...", process_rank, i + 1, start.elapsed().as_secs());
        }
    }

    // Add the HSPs of proteins with themselves
    if trivial_hsps {
        if !new_only {
            for protein in protein_set.iter() {
                let mut guard = hsps.lock().unwrap();
                let location = Location::new(protein.index(), 0);
                let hsp = HSP::new(location.clone(), location.clone(), protein.len());
                (*guard).insert(hsp);
            }
        } else {
            for protein in protein_set.iter() {
                if !protein.is_new() {
                    continue;
                }
                let mut guard = hsps.lock().unwrap();
                let location = Location::new(protein.index(), 0);
                let hsp = HSP::new(location.clone(), location.clone(), protein.len());
                (*guard).insert(hsp);
            }
        }
    }
   
    let guard = hsps.lock().unwrap();
    (*guard).iter()
        .map(|hsp| {
            HSP::new(hsp.location(0).clone(), hsp.location(1).clone(), hsp.len())
        })
    .collect::<HashSet<HSP>>()
}

/// Computes the HSPs that arise from the smer
pub unsafe fn compute_hsps_for_smer(
    smer_index: usize, 
     smer_list: &Vec<SmerCollection>,
     smer_map: &HashMap<u64, usize>,
     protein_set: &ProteinSet,
     seed: &Seed,
     kmer_size: usize,
     t_sim: i16,
     t_hit: i16,
     new_only: bool
 ) -> HashSet<HSP> {
    let mut hsps = HashSet::new();
    let smer = &smer_list[smer_index];

    let raw_similar = compute_similar_smers(smer.value(), &seed, t_sim);
    let similar_indices: Vec<usize> = raw_similar.into_iter()
        .filter(|x| x >= &smer.value() && smer_map.contains_key(&x))
        .map(|value| *smer_map.get(&value).unwrap())
        .collect();
    
    for index in similar_indices {
        let similar_smer = &smer_list[index];
        if smer != similar_smer {
            for smer1_location in smer.locations() {
                for smer2_location in similar_smer.locations() {
                    if new_only && !protein_set.is_new(smer1_location.index()) && !protein_set.is_new(smer2_location.index()) {
                        continue
                    }
                    let hits = find_hits(smer1_location, smer2_location, protein_set, kmer_size, t_hit);
                    let smer_hsps: Vec<HSP> = hits.iter().map(|hit| {
                        let hit = extend_hit(hit, protein_set, kmer_size, t_hit);
                        return hit;
                    }).collect();
                    hsps.extend(smer_hsps);
                }
            }
        } else {
            let locations = smer.locations().iter().cloned().collect::<Vec<Location>>();
            for i in 0..smer.locations().len() {
                for j in i+1..smer.locations().len() { 
                    if new_only && !protein_set.is_new(locations[i].index()) && !protein_set.is_new(locations[j].index()) {
                        continue
                    }
                    let hits = find_hits(&locations[i], &locations[j], protein_set, kmer_size, t_hit);
                    let smer_hsps: Vec<HSP> = hits.iter().map(|hit| {
                        let hit = extend_hit(hit, protein_set, kmer_size, t_hit);
                        return hit;
                    }).collect();
                    hsps.extend(smer_hsps);
                }
            }
        }
    }

    hsps
}

/// Retrieves hits around similar smers
pub unsafe fn find_hits(location1: &Location, location2: &Location, protein_set: &ProteinSet, kmer_size: usize, t_hit: i16) -> Vec<(usize, usize, usize, usize, i16)> {
    //println!("start hit");
    let mut hits = vec![];
    
    // Do not allow a hit for the two exact same locations, as these
    // will automatically extend to both ends of the proteins in both
    // direction and give the trivial HSP (P1 P1 0 0 LENGTH_OF_PROTEIN).
    if location1 == location2 {
        return hits;
    }

    let protein1 = protein_set.get_protein_by_id(location1.index());
    let length1 = protein1.len();
    let protein2 = protein_set.get_protein_by_id(location2.index());
    let length2 = protein2.len();
    
    for offset in 0..(kmer_size - 12 + 1) { // DO NOT

        let start1: i16 = location1.position() as i16 - offset as i16;
        let end1: usize = location1.position() + kmer_size - offset - 1;
        let start2: i16 = location2.position() as i16  - offset as i16;
        let end2: usize = location2.position() + kmer_size - offset - 1;

        if end1 >= length1 || end2 >= length2 {
            continue;
        }

        if start1 < 0 || start2 < 0 {
            break;
        }

        // If the subsequences do not meet the minimum score for this offset,
        // don't bother trying to extend this region
        let score = score_sequences(protein1, protein2, start1 as usize, start2 as usize, kmer_size);
        if score < t_hit {
            continue
        }

        if location1.index() < location2.index() {
            hits.push((location1.index(), location1.position() - offset, location2.index(), location2.position() - offset, score));
        } else {
            hits.push((location2.index(), location2.position() - offset, location1.index(), location1.position() - offset, score));
        }

        break;
    }

    hits
}

pub fn extend_hit(hit: &(usize, usize, usize, usize, i16), protein_set: &ProteinSet, kmer_size: usize, t_hit: i16) -> HSP {

    let protein1 = protein_set.get_protein_by_id(hit.0);
    let protein2 = protein_set.get_protein_by_id(hit.2);
    let mut new_sta1 = hit.1;
    let mut new_sta2 = hit.3;
    let mut new_end1 = hit.1 + kmer_size - 1;
    let mut pos_1l = hit.1;
    let mut pos_2l = hit.3;
    let mut pos_1r = new_end1;
    let mut pos_2r = hit.3 + kmer_size - 1;
    let to_left: usize;
    let to_right: usize;
    let mut current_score = hit.4;

    if protein1.len() - hit.1 > protein2.len() - hit.3 {
        to_right = protein2.len() - (hit.3 + kmer_size);
    } else {
        to_right = protein1.len() - (hit.1 + kmer_size); 
    } 

    for i in 0..to_right {
        current_score = current_score
            - score_position(protein1, protein2, pos_1l, pos_2l)
            + score_position(protein1, protein2, pos_1r + 1, pos_2r + 1);

        if current_score >= t_hit {
            pos_1l += 1;
            pos_2l += 1;
            pos_1r += 1;
            pos_2r += 1;
            if i == (to_right - 1) {
                new_end1 = pos_1r;
            }
        } else {
            new_end1 = pos_1r;
            break;
        }
    }

    pos_1l = hit.1;
    pos_2l = hit.3;
    pos_1r = hit.1 + kmer_size - 1;
    pos_2r = hit.3 + kmer_size - 1;
    current_score = hit.4;

    if hit.1 > hit.3 {
        to_left = hit.3;
    } else {
        to_left = hit.1;
    }

    for i in 0..to_left {
        current_score = current_score 
            - score_position(protein1, protein2, pos_1r, pos_2r)
            + score_position(protein1, protein2, pos_1l - 1, pos_2l - 1);

        if current_score >= t_hit {
            pos_1l -= 1;
            pos_2l -= 1;
            pos_1r -= 1;
            pos_2r -= 1;
            if i == (to_left - 1) {
                new_sta1 = pos_1l;
                new_sta2 = pos_2l;
            }
        } else {
            new_sta1 = pos_1l;
            new_sta2 = pos_2l;
            break;
        }
    }

    let length = new_end1 - new_sta1 + 1;
    HSP::new(Location::new(hit.0, new_sta1), Location::new(hit.2, new_sta2), length)
}
