use std::collections::HashSet;

use crate::sprint::constants::CODE_RESIDUE_MAP;
use crate::sprint::location::Location;
use crate::sprint::seed::Seed;

#[derive(Hash, Clone)]
pub struct Smer {
    value: u64,
    location: Location,
}

impl Smer {
    pub fn new(value: u64, location: Location) -> Self {
        Smer { value, location }
    }

    pub fn value(&self) -> u64 {
        self.value
    }

    pub fn location(&self) -> &Location {
        &self.location
    }

    pub fn get_aa_index_at(value: u64, position: usize, seed: &Seed) -> usize {
        ((value >> ((seed.len() - position - 1) * 5)) & 31) as usize
    }

    pub fn as_string(value: u64, seed: &Seed) -> String {
        let mut sequence = "".to_string();
        for i in 0..seed.len() {
            let value = value >> (5u64 * (seed.len() as u64 - i as u64 - 1)) & 31u64;
            let amino_acid = *CODE_RESIDUE_MAP.get(&(value as u16)).unwrap();
            sequence.push(amino_acid);
        }
        sequence
    }

    pub fn to_sequence(&self, seed: &Seed) -> String {
        let mut sequence = "".to_string();
        for i in 0..seed.len() {
            let value = self.value >> (5u64 * (seed.len() as u64 - i as u64 - 1)) & 31u64;
            let amino_acid = *CODE_RESIDUE_MAP.get(&(value as u16)).unwrap();
            sequence.push(amino_acid);
        }
        sequence
    }

    pub fn mutate(value: u64, position: usize, aa_index: u64, seed: &Seed) -> u64 {
        let original = value & (31 << ((seed.len() - position - 1) * 5));
        let new_value = value - original + (aa_index << ((seed.len() - position - 1) * 5));
        new_value
    }
}

pub struct SmerCollection {
    value: u64,
    locations: HashSet<Location>,
}

impl SmerCollection {
    pub fn new(value: u64, locations: HashSet<Location>) -> Self {
        SmerCollection { value, locations }
    }

    pub fn value(&self) -> u64 {
        self.value
    }

    pub fn locations(&self) -> &HashSet<Location> {
        &self.locations
    }
}

impl std::cmp::Ord for SmerCollection {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.value.cmp(&other.value)
    }
}

impl PartialOrd for SmerCollection {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl std::cmp::PartialEq for SmerCollection {
    fn eq(&self, other: &Self) -> bool {
        self.value == other.value
    }
}

impl std::cmp::Eq for SmerCollection {}

