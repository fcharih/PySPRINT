use std::collections::HashMap;
use std::io::Error;
use crate::sprint::protein::Protein;

use crate::sprint::fileio::load_fasta;

/// Structure that holds the sequences used to extract HSPs
/// and score protein interactions
#[derive(Clone)]
pub struct ProteinSet {
    // Set of proteins sequences
    proteins: Vec<Protein>,
    // Indices of the proteins
    indices: HashMap<String, usize>
}

/// Builds a name-index map from a vector of Protein
fn get_indices(sequences: &Vec<Protein>) -> HashMap<String, usize> {
    let indices: HashMap<String, usize> = sequences
        .iter()
        .enumerate()
        .map(|(i, protein)| (protein.name(), i))
        .collect();
    indices
}

impl ProteinSet {

    /// Create a protein set from protein sequences
    pub fn new(proteins: Vec<Protein>) -> Self {
        let indices = get_indices(&proteins);

        ProteinSet {
            proteins,
            indices
        }
    }

    /// Create a protein set from a FASTA file
    pub fn from_file(filepath: &str) -> Result<Self, Error> {
        let proteins = load_fasta(filepath, false)?;
        let indices = get_indices(&proteins);

        Ok(ProteinSet {
            proteins,
            indices
        })
    }

    /// Add new proteins
    pub fn add_new(&mut self, proteins: Vec<Protein>, new: bool) {
        let mut new_set = vec![];

        // Add the old proteins
        for p in self.proteins.iter() {
            new_set.push(p.clone());
        }

        // Add the new proteins
        for (i, p) in proteins.iter().enumerate() {
            let new_prot = Protein::new(i + self.proteins.len(), p.name(), p.seq(), new);
            new_set.push(new_prot);
        }

        // Reset the indices table
        let indices = get_indices(&new_set);

        // Update
        self.proteins = new_set;
        self.indices = indices;
    }

    // Add new proteins from a file
    pub fn add_from_file(&mut self, filepath: &String, new: bool) {
        let proteins = load_fasta(filepath, new)
            .expect("Error while opening sequences file.");
        self.add_new(proteins, new);
    }

    pub fn get_protein_by_id(&self, index: usize) -> &Protein {
        &self.proteins[index]
    }

    pub fn get_protein_by_name(&self, name: &String) -> &Protein {
        &self.proteins[*self.indices.get(name).unwrap()]
    }

    pub fn iter(&self) -> impl Iterator<Item = &Protein> {
        self.proteins.iter()
    }

    pub fn is_new(&self, index: usize) -> bool {
        self.get_protein_by_id(index).is_new()
    }

    pub fn len(&self) -> usize {
        self.proteins.len()
    }

    pub fn contains(&self, name: &String) -> bool {
        self.indices.contains_key(name)
    }
}


