use std::collections::HashSet;
use crate::rsprint::proteinset::ProteinSet;
use crate::rsprint::hsp::HSP;

use super::location::Location;

pub fn process_hsps(
    protein_set: &ProteinSet,
    hsps: HashSet<HSP>,
    kmer_size: usize,
    t_count: u16,
    verbose: bool
) -> HashSet<HSP> {

    // Retrieve the HSPs
    if verbose {
        println!("Collecting the HSPs...");
    }

    // Count the occurence of amino acids within HSPs
    if verbose {
        println!("Counting residue occurences within HSPs...");
    }
    let counts = count_residue_occurences(&hsps, &protein_set, kmer_size);

    // Process the HSPs
    if verbose {
        println!("Processing HSPs...");
    }
    let mut processed_hsps: HashSet<HSP> = HashSet::new();
    hsps.into_iter().for_each(|hsp| unsafe {
        process_hsp(hsp, &protein_set, &counts, kmer_size, t_count, &mut processed_hsps);
    });

    processed_hsps.into_iter().map(|hsp| {
        hsp.clone()
    }).collect::<HashSet<HSP>>()

}

pub fn count_residue_occurences(hsps: &HashSet<HSP>, protein_set: &ProteinSet, kmer_size: usize) -> Vec<Vec<u16>>{
    let protein_set_size= protein_set.len();
    let mut counts: Vec<Vec<u16>> = Vec::new();

    for _ in 0..protein_set_size {
        counts.push(Vec::new());
    }

    for i in 0..protein_set_size {
        for _ in 0..protein_set.get_protein_by_id(i).len() + 1 {
            counts[i].push(0);
        }
    }

    for hsp in hsps {
        if hsp.len() < kmer_size {
            continue;
        }
        for i in 0..hsp.len() - (kmer_size - 1) {
            //println!("{}", hsp.as_string(&sequence_set).unwrap());
            counts[hsp.location(0).index()][i + hsp.location(0).position()] += 1;
            counts[hsp.location(1).index()][i + hsp.location(1).position()] += 1;
        }
    }
    counts
}

pub unsafe fn process_hsp(
    hsp: HSP,
    protein_set: &ProteinSet,
    counts: &Vec<Vec<u16>>,
    kmer_size: usize,
    t_count: u16,
    output: &mut HashSet<HSP>
) {

    if hsp.len() < kmer_size {
        return;
    }

    for i in 0..hsp.len() - (kmer_size - 1) {
        let last_position = i == hsp.len() - kmer_size;
        let count1 = counts[hsp.location(0).index()][hsp.location(0).position() + i];
        let count2 = counts[hsp.location(1).index()][hsp.location(1).position() + i];

        if last_position {
            if count1 <= t_count && count2 <= t_count {
                output.insert(hsp);
                break;
            }
        }

        if count1 > t_count || count2 > t_count {
            if i != 0 {
                let new = HSP::new(Location::new(hsp.location(0).index(), hsp.location(0).position()),
                    Location::new(hsp.location(1).index(), hsp.location(1).position()), 
                    kmer_size - 1 + i);
                output.insert(new);
            }

            if (hsp.len() - i) >= (kmer_size + 1) {
                let shorter_hsp = HSP::new(
                    Location::new(hsp.location(0).index(), hsp.location(0).position() + i + 1),
                    Location::new(hsp.location(1).index(), hsp.location(1).position() + i + 1),
                    hsp.len() - i - 1);

                process_hsp(shorter_hsp, &protein_set, &counts, kmer_size, t_count, output);
                break;
            }
        }
    }
}