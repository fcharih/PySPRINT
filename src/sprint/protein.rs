use crate::sprint::constants::RESIDUE_CODES;

#[derive(Clone)]
pub struct Protein {
    index: usize,
    name: String,
    sequence: String,
    residues: Vec<usize>,
    new: bool
}

impl Protein {
    /// Creates a new protein
    pub fn new(index: usize, name: String, sequence: String, new: bool) -> Self {
        Protein { 
            index,
            name: name.clone(),
            sequence: sequence.clone(),
            residues: sequence.chars()
                .map(|amino_acid| *RESIDUE_CODES.get(&amino_acid).unwrap())
                .collect(),
            new
        }
    }

    pub fn sub(&self, start: usize, end: usize) -> &str {
        &self.sequence[start as usize..end + 1 as usize]
    }

    pub fn seq(&self) -> String {
        self.sequence.clone()
    }

    pub fn name(&self) -> String {
        self.name.clone()
    }

    pub fn index(&self) -> usize {
        self.index
    }

    pub fn len(&self) -> usize {
        self.residues.len()
    }

    pub fn residue(&self, index: usize) -> usize {
        self.residues[index]
    }

    pub fn is_new(&self) -> bool {
        self.new
    }
}
