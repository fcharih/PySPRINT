use std::collections::HashMap;

// The seeds that were determined to be optimal by
// the authors of SPRINT
pub static SEEDS: [&str; 4] = [
    "11****11***1",
    "1**1*1***1*1",
    "11**1***1**1",
    "1*1******111"
];

lazy_static! {

    /// Maps an amino acid to its index in the PAM120 matrix
    pub static ref RESIDUE_CODES: HashMap<char, usize> = [
        ('A', 1),  ('a', 1),
        ('R', 2),  ('r', 2),
        ('N', 3),  ('n', 3),
        ('D', 4),  ('d', 4),
        ('C', 5),  ('c', 5),
        ('Q', 6),  ('q', 6),
        ('E', 7),  ('e', 7),
        ('G', 8),  ('g', 8),
        ('H', 9),  ('h', 9),
        ('I', 10),  ('i', 10),
        ('L', 11), ('l', 11),
        ('K', 12), ('k', 12),
        ('M', 13), ('m', 13),
        ('F', 14), ('f', 14),
        ('P', 15), ('p', 15),
        ('S', 16), ('s', 16),
        ('T', 17), ('t', 17),
        ('W', 18), ('w', 18),
        ('Y', 19), ('y', 19),
        ('V', 20), ('v', 20),
        ('B', 20), ('b', 20),
        ('Z', 20), ('z', 20),
        ('X', 20), ('x', 20),
        ('U', 20), ('u', 20),
        ('O', 20), ('o', 20),
    ].iter().cloned().collect();

    
    /// Maps an index to the corresponding amino acid
    pub static ref CODE_RESIDUE_MAP: HashMap<u16, char> = [
        (0, '-'), 
        (1, 'A'), 
        (2, 'R'), 
        (3, 'N'), 
        (4, 'D'), 
        (5, 'C'), 
        (6, 'Q'), 
        (7, 'E'), 
        (8, 'G'), 
        (9, 'H'), 
        (10, 'I'), 
        (11, 'L'), 
        (12, 'K'), 
        (13, 'M'), 
        (14, 'F'), 
        (15, 'P'), 
        (16, 'S'), 
        (17, 'T'), 
        (18, 'W'), 
        (19, 'Y'), 
        (20, 'V'),  // assume a V if random letter (same within the algo)
        (21, 'X')
    ].iter().cloned().collect();
}

/// PAM120 matrix
/// Note: The first row/column correspond to non-important positions
/// (indicated by a "-" in the seed)
pub static PAM120: [[i16; 24]; 24] = [
    [ 0, 0, 0, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 ],
    [ 0,  3, -3, -1,  0, -3, -1,  0,  1, -3, -1, -3, -2, -2, -4,  1,  1,  1, -7, -4,  0,  0, -1, -1 ],
    [ 0, -3,  6, -1, -3, -4,  1, -3, -4,  1, -2, -4,  2, -1, -5, -1, -1, -2,  1, -5, -3, -2, -1, -2 ],
    [ 0, -1, -1,  4,  2, -5,  0,  1,  0,  2, -2, -4,  1, -3, -4, -2,  1,  0, -4, -2, -3,  3,  0, -1 ],
    [ 0,  0, -3,  2,  5, -7,  1,  3,  0,  0, -3, -5, -1, -4, -7, -3,  0, -1, -8, -5, -3,  4,  3, -2 ],
    [ 0, -3, -4, -5, -7,  9, -7, -7, -4, -4, -3, -7, -7, -6, -6, -4,  0, -3, -8, -1, -3, -6, -7, -4 ],
    [ 0, -1,  1,  0,  1, -7,  6,  2, -3,  3, -3, -2,  0, -1, -6,  0, -2, -2, -6, -5, -3,  0,  4, -1 ],
    [ 0,  0, -3,  1,  3, -7,  2,  5, -1, -1, -3, -4, -1, -3, -7, -2, -1, -2, -8, -5, -3,  3,  4, -1 ],
    [ 0,  1, -4,  0,  0, -4, -3, -1,  5, -4, -4, -5, -3, -4, -5, -2,  1, -1, -8, -6, -2,  0, -2, -2 ],
    [ 0, -3,  1,  2,  0, -4,  3, -1, -4,  7, -4, -3, -2, -4, -3, -1, -2, -3, -3, -1, -3,  1,  1, -2 ],
    [ 0, -1, -2, -2, -3, -3, -3, -3, -4, -4,  6,  1, -3,  1,  0, -3, -2,  0, -6, -2,  3, -3, -3, -1 ],
    [ 0, -3, -4, -4, -5, -7, -2, -4, -5, -3,  1,  5, -4,  3,  0, -3, -4, -3, -3, -2,  1, -4, -3, -2 ],
    [ 0, -2,  2,  1, -1, -7,  0, -1, -3, -2, -3, -4,  5,  0, -7, -2, -1, -1, -5, -5, -4,  0, -1, -2 ],
    [ 0, -2, -1, -3, -4, -6, -1, -3, -4, -4,  1,  3,  0,  8, -1, -3, -2, -1, -6, -4,  1, -4, -2, -2 ],
    [ 0, -4, -5, -4, -7, -6, -6, -7, -5, -3,  0,  0, -7, -1,  8, -5, -3, -4, -1,  4, -3, -5, -6, -3 ],
    [ 0,  1, -1, -2, -3, -4,  0, -2, -2, -1, -3, -3, -2, -3, -5,  6,  1, -1, -7, -6, -2, -2, -1, -2 ],
    [ 0,  1, -1,  1,  0,  0, -2, -1,  1, -2, -2, -4, -1, -2, -3,  1,  3,  2, -2, -3, -2,  0, -1, -1 ],
    [ 0,  1, -2,  0, -1, -3, -2, -2, -1, -3,  0, -3, -1, -1, -4, -1,  2,  4, -6, -3,  0,  0, -2, -1 ],
    [ 0, -7,  1, -4, -8, -8, -6, -8, -8, -3, -6, -3, -5, -6, -1, -7, -2, -6, 12, -2, -8, -6, -7, -5 ],
    [ 0, -4, -5, -2, -5, -1, -5, -5, -6, -1, -2, -2, -5, -4,  4, -6, -3, -3, -2,  8, -3, -3, -5, -3 ],
    [ 0,  0, -3, -3, -3, -3, -3, -3, -2, -3,  3,  1, -4,  1, -3, -2, -2,  0, -8, -3,  5, -3, -3, -1 ],
    [ 0,  0, -2,  3,  4, -6,  0,  3,  0,  1, -3, -4,  0, -4, -5, -2,  0,  0, -6, -3, -3,  4,  2, -1 ],
    [ 0, -1, -1,  0,  3, -7,  4,  4, -2,  1, -3, -3, -1, -2, -6, -1, -1, -2, -7, -5, -3,  2,  4, -1 ],
    [ 0, -1, -2, -1, -2, -4, -1, -1, -2, -2, -1, -2, -2, -2, -3, -2, -1, -1, -5, -3, -1, -1, -1, -2 ],
];

/// Standard PAM120 matrix (without the ``don't care`` row and column)
pub static PAM120_SPRINT: [[i16; 23]; 23] = [
    [  3, -3, -1,  0, -3, -1,  0,  1, -3, -1, -3, -2, -2, -4,  1,  1,  1, -7, -4,  0,  0, -1, -1 ],
    [ -3,  6, -1, -3, -4,  1, -3, -4,  1, -2, -4,  2, -1, -5, -1, -1, -2,  1, -5, -3, -2, -1, -2 ],
    [ -1, -1,  4,  2, -5,  0,  1,  0,  2, -2, -4,  1, -3, -4, -2,  1,  0, -4, -2, -3,  3,  0, -1 ],
    [  0, -3,  2,  5, -7,  1,  3,  0,  0, -3, -5, -1, -4, -7, -3,  0, -1, -8, -5, -3,  4,  3, -2 ],
    [ -3, -4, -5, -7,  9, -7, -7, -4, -4, -3, -7, -7, -6, -6, -4,  0, -3, -8, -1, -3, -6, -7, -4 ],
    [ -1,  1,  0,  1, -7,  6,  2, -3,  3, -3, -2,  0, -1, -6,  0, -2, -2, -6, -5, -3,  0,  4, -1 ],
    [  0, -3,  1,  3, -7,  2,  5, -1, -1, -3, -4, -1, -3, -7, -2, -1, -2, -8, -5, -3,  3,  4, -1 ],
    [  1, -4,  0,  0, -4, -3, -1,  5, -4, -4, -5, -3, -4, -5, -2,  1, -1, -8, -6, -2,  0, -2, -2 ],
    [ -3,  1,  2,  0, -4,  3, -1, -4,  7, -4, -3, -2, -4, -3, -1, -2, -3, -3, -1, -3,  1,  1, -2 ],
    [ -1, -2, -2, -3, -3, -3, -3, -4, -4,  6,  1, -3,  1,  0, -3, -2,  0, -6, -2,  3, -3, -3, -1 ],
    [ -3, -4, -4, -5, -7, -2, -4, -5, -3,  1,  5, -4,  3,  0, -3, -4, -3, -3, -2,  1, -4, -3, -2 ],
    [ -2,  2,  1, -1, -7,  0, -1, -3, -2, -3, -4,  5,  0, -7, -2, -1, -1, -5, -5, -4,  0, -1, -2 ],
    [ -2, -1, -3, -4, -6, -1, -3, -4, -4,  1,  3,  0,  8, -1, -3, -2, -1, -6, -4,  1, -4, -2, -2 ],
    [ -4, -5, -4, -7, -6, -6, -7, -5, -3,  0,  0, -7, -1,  8, -5, -3, -4, -1,  4, -3, -5, -6, -3 ],
    [  1, -1, -2, -3, -4,  0, -2, -2, -1, -3, -3, -2, -3, -5,  6,  1, -1, -7, -6, -2, -2, -1, -2 ],
    [  1, -1,  1,  0,  0, -2, -1,  1, -2, -2, -4, -1, -2, -3,  1,  3,  2, -2, -3, -2,  0, -1, -1 ],
    [  1, -2,  0, -1, -3, -2, -2, -1, -3,  0, -3, -1, -1, -4, -1,  2,  4, -6, -3,  0,  0, -2, -1 ],
    [ -7,  1, -4, -8, -8, -6, -8, -8, -3, -6, -3, -5, -6, -1, -7, -2, -6, 12, -2, -8, -6, -7, -5 ],
    [ -4, -5, -2, -5, -1, -5, -5, -6, -1, -2, -2, -5, -4,  4, -6, -3, -3, -2,  8, -3, -3, -5, -3 ],
    [  0, -3, -3, -3, -3, -3, -3, -2, -3,  3,  1, -4,  1, -3, -2, -2,  0, -8, -3,  5, -3, -3, -1 ],
    [  0, -2,  3,  4, -6,  0,  3,  0,  1, -3, -4,  0, -4, -5, -2,  0,  0, -6, -3, -3,  4,  2, -1 ],
    [ -1, -1,  0,  3, -7,  4,  4, -2,  1, -3, -3, -1, -2, -6, -1, -1, -2, -7, -5, -3,  2,  4, -1 ],
    [ -1, -2, -1, -2, -4, -1, -1, -2, -2, -1, -2, -2, -2, -3, -2, -1, -1, -5, -3, -1, -1, -1, -2 ],
];

/// Indices of PAM120 matrix for similarity computation
/// TODO: Reverse engineer this...
pub static PAM120_ORDERED_SPRINT: [[usize; 21]; 21] = [
    [ 0, 7, 14, 15, 16, 3, 6, 19, 20, 2, 5, 9, 11, 12, 1, 4, 8, 10, 13, 18, 17, ],
    [ 1, 11, 5, 8, 17, 2, 12, 14, 15, 9, 16, 20, 0, 3, 6, 19, 4, 7, 10, 13, 18, ],
    [ 2, 20, 3, 8, 6, 11, 15, 5, 7, 16, 0, 1, 9, 14, 18, 12, 19, 10, 13, 17, 4, ],
    [ 3, 20, 6, 2, 5, 0, 7, 8, 15, 11, 16, 1, 9, 14, 19, 12, 10, 18, 4, 13, 17, ],
    [ 4, 15, 18, 0, 9, 16, 19, 1, 7, 8, 14, 2, 12, 13, 20, 3, 5, 6, 10, 11, 17, ],
    [ 5, 8, 6, 1, 3, 2, 11, 14, 20, 0, 12, 10, 15, 16, 7, 9, 19, 18, 13, 17, 4, ],
    [ 6, 3, 20, 5, 2, 0, 7, 8, 11, 15, 14, 16, 1, 9, 12, 19, 10, 18, 4, 13, 17, ],
    [ 7, 0, 15, 2, 3, 20, 6, 16, 14, 19, 5, 11, 1, 4, 8, 9, 12, 10, 13, 18, 17, ],
    [ 8, 5, 2, 1, 20, 3, 6, 14, 18, 11, 15, 0, 10, 13, 16, 17, 19, 4, 7, 9, 12, ],
    [ 9, 19, 10, 12, 13, 16, 0, 1, 2, 15, 18, 3, 4, 5, 6, 11, 14, 20, 7, 8, 17, ],
    [ 10, 12, 9, 19, 13, 5, 18, 0, 8, 14, 16, 17, 1, 2, 6, 11, 15, 20, 3, 7, 4, ],
    [ 11, 1, 2, 5, 12, 20, 3, 6, 15, 16, 0, 8, 14, 7, 9, 10, 19, 17, 18, 4, 13, ],
    [ 12, 10, 9, 19, 11, 1, 5, 13, 16, 0, 15, 2, 6, 14, 3, 7, 8, 18, 20, 4, 17, ],
    [ 13, 18, 9, 10, 12, 17, 8, 15, 19, 0, 2, 16, 1, 7, 14, 20, 4, 5, 3, 6, 11, ],
    [ 14, 0, 15, 5, 1, 8, 16, 2, 6, 7, 11, 19, 20, 3, 9, 10, 12, 4, 13, 18, 17, ],
    [ 15, 16, 0, 2, 7, 14, 3, 4, 20, 1, 6, 11, 5, 8, 9, 12, 17, 19, 13, 18, 10, ],
    [ 16, 15, 0, 2, 9, 19, 20, 3, 7, 11, 12, 14, 1, 5, 6, 4, 8, 10, 18, 13, 17, ],
    [ 17, 1, 13, 15, 18, 8, 10, 2, 11, 5, 9, 12, 16, 20, 0, 14, 3, 4, 6, 7, 19, ],
    [ 18, 13, 4, 8, 2, 9, 10, 17, 15, 16, 19, 20, 0, 12, 1, 3, 5, 6, 11, 7, 14, ],
    [ 19, 9, 10, 12, 0, 16, 7, 14, 15, 1, 2, 3, 4, 5, 6, 8, 13, 18, 20, 11, 17, ],
    [ 3, 20, 2, 6, 8, 0, 5, 7, 11, 15, 16, 1, 14, 9, 18, 19, 10, 12, 13, 4, 17, ]
];