
#[derive(Hash, Clone)]
pub struct Location {
    protein_index: usize,
    position: usize
}

impl Location {
    pub fn new(protein_index: usize, position: usize) -> Self {
        Location {
            protein_index,
            position
        }
    }

    pub fn index(&self) -> usize {
        self.protein_index
    }

    pub fn position(&self) -> usize {
        self.position
    }
}

impl PartialEq for Location {
    fn eq(&self, other: &Self) -> bool {
        self.protein_index == other.protein_index &&
        self.position == other.position
    }
}

impl Eq for Location {}