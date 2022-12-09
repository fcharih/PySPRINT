use std::collections::HashSet;
use crate::sprint::seed::Seed;
use crate::sprint::smer::Smer;
use crate::sprint::constants::{PAM120, PAM120_SPRINT, PAM120_ORDERED_SPRINT};
use crate::sprint::utils::{convert_smer_to_sprint, convert_smer_to_pysprint};

pub fn compare_smers(smer1: &Smer, smer2: &Smer, seed: &Seed) -> i16 {
    let mut score = 0;
    for i in 0..seed.len() {
        let offset = i * 5;
        // The AND operation with the seed will ensure that the resulting
        // amino acid index will be 0 if the seed for the corresponding
        // position is null. This will yield a score of 0.
        let amino_acid1 = ((smer1.value() & seed.value()) >> offset) & 31;
        let amino_acid2 = ((smer2.value() & seed.value()) >> offset) & 31;
        score += PAM120[amino_acid1 as usize][amino_acid2 as usize];
    }
    score
}

pub fn compute_similar_smers(smer: u64, seed: &Seed, t_hit: i16) -> HashSet<u64> {
    let sprint_smer = convert_smer_to_sprint(smer);
    let mut similar_smers: HashSet<u64> = HashSet::new();
    let mut smer_score: i16 = 0;
    let mut seed_current_digit: u64;
    let mut smer_current_digit: u64;

    for i in 1..seed.len() {
        seed_current_digit = seed.value() & (31 << (i * 5));
        if seed_current_digit != 0 {
            smer_current_digit = sprint_smer & (31 << (i * 5));
            smer_current_digit = smer_current_digit >> (i * 5);
            smer_score += PAM120_SPRINT[smer_current_digit as usize][smer_current_digit as usize];
        }
    }

    change_smer_digit(sprint_smer, 0, seed, t_hit - smer_score, &mut similar_smers);
    similar_smers
}

pub fn change_smer_digit(smer: u64, digit_pos: usize, seed: &Seed, score_needed: i16, similar_smers: &mut HashSet<u64>) {
    if digit_pos == seed.len() {
        return;
    }
    let seed_current_digit: u64 = seed.value() & (31 << (digit_pos * 5));
    if seed_current_digit == 0 {
        change_smer_digit(smer, digit_pos + 1, seed, score_needed, similar_smers);
    } else if digit_pos == seed.len() - 1 {
        let smer_cunt_dig = (smer & (31 << (digit_pos * 5))) >> (digit_pos * 5);
        let mut y: u64;

        for i in 0..20 {
            if PAM120_SPRINT[smer_cunt_dig as usize][PAM120_ORDERED_SPRINT[smer_cunt_dig as usize][i]] >= score_needed {
                y = smer & (!(31 << (digit_pos * 5)));
                y = y | (PAM120_ORDERED_SPRINT[smer_cunt_dig as usize][i] << (5 * digit_pos)) as u64;
                similar_smers.insert(convert_smer_to_pysprint(y, seed));
            } else {
                break;
            }
        }
    } else {
        let smer_cunt_dig: u64 = (smer & (31 << (digit_pos * 5))) >> (digit_pos * 5);
        let smer_cunt_dig_plus1: u64 = get_smer_dig_plus1(smer, digit_pos, &seed);
        let mut y: u64;
        let mut new_score_needed: i16;
        for i in 0..20 {
            if PAM120_SPRINT[smer_cunt_dig as usize][PAM120_ORDERED_SPRINT[smer_cunt_dig as usize][i]] >= score_needed {
                y = smer & (!(31 << (digit_pos * 5)));
                y = y | (PAM120_ORDERED_SPRINT[smer_cunt_dig as usize][i] << (5 * digit_pos)) as u64;
                similar_smers.insert(convert_smer_to_pysprint(y, seed));
                new_score_needed = score_needed + PAM120_SPRINT[smer_cunt_dig_plus1 as usize][smer_cunt_dig_plus1 as usize] - PAM120_SPRINT[smer_cunt_dig as usize][PAM120_ORDERED_SPRINT[smer_cunt_dig as usize][i]];
                change_smer_digit(y, digit_pos + 1, seed, new_score_needed, similar_smers);
            } else {
                break;
            }
        }
    }
}

pub fn get_smer_dig_plus1(smer: u64, digit_pos: usize, seed: &Seed) -> u64 {
    let seed_current_digit = seed.value() & (31 << ((digit_pos + 1) * 5));
    if seed_current_digit != 0 {
        return (smer & (31 << ((digit_pos + 1) * 5))) >> ((digit_pos + 1) * 5);
    } else {
        return get_smer_dig_plus1(smer, digit_pos + 1, seed);
    }
}
