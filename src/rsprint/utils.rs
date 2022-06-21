use crate::rsprint::seed::Seed;

pub fn matrix_max(matrix: &Vec<Vec<f32>>) -> f32 {
    let mut max = 0f32;
    for row in matrix.iter() {
        for entry in row {
            if *entry > max {
                max = *entry;
            }
        }
    }
    max
}

/// Converts an s-mer from rSPRINT format to SPRINT format, i.e. with
/// -1 for every "matter" position
pub fn convert_smer_to_sprint(smer: u64) -> u64 {
    let mut sprint_smer: u64 = 0;
    for i in 0..12 {
        let val = (smer >> (i * 5)) & 31;
        let amino_acid_index = if val == 0 { 0 } else { val - 1 };
        sprint_smer += amino_acid_index << (i * 5);
    }
    sprint_smer
}

/// Converts an s-mer in SPRINT format to rSPRINT format, i.e. with
/// +1 for every "matter" position
pub fn convert_smer_to_pysprint(smer: u64, seed: &Seed) -> u64 {
    let mut pysprint_smer: u64 = 0;
    for i in 0..12 {
        if ((seed.value() >> (i * 5)) & 31) == 0 {
            continue;
        }
        let amino_acid_index = ((smer >> (i * 5)) & 31) + 1;
        pysprint_smer += amino_acid_index << (i * 5);
    }
    pysprint_smer
}
