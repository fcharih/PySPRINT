<<<<<<< HEAD
use std::collections::HashSet;

use pyo3::prelude::*;

use super::extraction::extract_hsps;
use super::fileio::{load_hsps, load_pairs};
use super::hsp::HSP;
use super::prediction::score_interactions;
use super::protein::Protein;
use super::proteinset::ProteinSet;

#[pyclass]
pub struct SPRINT {
    rank: usize,
    world_size: usize,
    protein_set: ProteinSet,
    hsps: HashSet<HSP>,
    training_pairs: Vec<(String, String)>,
}

#[pymethods]
impl SPRINT {
    #[new]
    pub fn new(
        rank: usize,
        world_size: usize,
        proteins_file: String,
        hsps_file: String,
        training_pairs_file: String,
    ) -> Self {
        let protein_set = ProteinSet::from_file(&proteins_file).unwrap();
        let hsps = load_hsps(&hsps_file, &protein_set);
        let training_pairs = load_pairs(&training_pairs_file);

        SPRINT {
            rank,
            world_size,
            protein_set,
            hsps,
            training_pairs,
        }
    }

    #[args(kmer_size = "20", t_smer = "15", t_hsp = "35")]
    pub fn extract_peptide_hsps(
        &self,
        peptides: Vec<(String, String)>,
        kmer_size: usize,
        t_smer: i16,
        t_hsp: i16,
    ) -> PyResult<HashSet<HSP>> {
        let new_proteins: Vec<Protein> = peptides
            .iter()
            .enumerate()
            .map(|(i, pep)| {
                Protein::new(
                    i + self.protein_set.len(),
                    pep.0.clone(),
                    pep.1.clone(),
                    true,
                )
            })
            .collect();
        let mut protein_set_copy = self.protein_set.clone();
        protein_set_copy.add_new(new_proteins);
        Ok(extract_hsps(
            &protein_set_copy,
            kmer_size,
            t_smer,
            t_hsp,
            self.rank,
            self.world_size,
            true,
            true,
            false,
        ))
    }

    #[args(kmer_size = "20")]
    pub fn score_interactions(
        &self,
        peptides: Vec<(String, String)>,
        peptide_hsps: HashSet<HSP>,
        kmer_size: usize,
    ) -> PyResult<Vec<f32>> {
        let new_proteins: Vec<Protein> = peptides
            .iter()
            .enumerate()
            .map(|(i, pep)| {
                Protein::new(
                    i + self.protein_set.len(),
                    pep.0.clone(),
                    pep.1.clone(),
                    true,
                )
            })
            .collect();
        let mut protein_set_copy = self.protein_set.clone();
        protein_set_copy.add_new(new_proteins);

        let mut hsps = HashSet::new();
        self.hsps.iter().for_each(|hsp| {
            hsps.insert(hsp.clone());
        });
        peptide_hsps.iter().for_each(|hsp| {
            hsps.insert(hsp.clone());
        });

        Ok(score_interactions(
            &protein_set_copy,
            &hsps,
            &self.training_pairs,
            kmer_size,
            self.rank,
            self.world_size,
            false,
        ))
    }
}

=======
//use std::collections::HashSet;
//
//use pyo3::prelude::*;
//
//use super::extraction::extract_hsps;
//use super::fileio::{load_hsps, load_pairs};
//use super::prediction::score_interactions;
//use super::proteinset::ProteinSet;
//use super::protein::Protein;
//use super::hsp::HSP;
//
//#[pyclass]
//pub struct SPRINT {
//    rank: usize,
//    world_size: usize,
//    protein_set: ProteinSet,
//    hsps: HashSet<HSP>,
//    training_pairs: Vec<(String, String)>
//}
//
//#[pymethods]
//impl SPRINT {
//    
//    #[new]
//    pub fn new(rank: usize, world_size: usize, proteins_file: String, hsps_file: String, training_pairs_file: String) -> Self {
//        let protein_set = ProteinSet::from_file(&proteins_file).unwrap();
//        let hsps = load_hsps(&hsps_file, &protein_set);
//        let training_pairs = load_pairs(&training_pairs_file);
//
//        SPRINT {
//            rank,
//            world_size,
//            protein_set,
//            hsps,
//            training_pairs
//        }
//    }
//
//    #[args(
//        kmer_size = "20",
//        t_smer = "15",
//        t_hsp = "35"
//    )]
//    pub fn extract_peptide_hsps(&self, peptides: Vec<(String, String)>, kmer_size: usize, t_smer: i16, t_hsp: i16) -> PyResult<HashSet<HSP>> {
//        let new_proteins: Vec<Protein> = peptides.iter()
//            .enumerate()
//            .map(|(i, pep)| {
//                Protein::new(i + self.protein_set.len(), pep.0.clone(), pep.1.clone(), true)
//            }).collect();
//        let mut protein_set_copy = self.protein_set.clone();
//        protein_set_copy.add_new(new_proteins);
//        Ok(extract_hsps(&protein_set_copy, kmer_size, t_smer, t_hsp, self.rank, self.world_size, true, true, false))
//    }
//
//    #[args(
//        kmer_size = "20",
//    )]
//    pub fn score_interactions(&self, peptides: Vec<(String, String)>, peptide_hsps: HashSet<HSP>, kmer_size: usize) -> PyResult<Vec<f32>> {
//        let new_proteins: Vec<Protein> = peptides.iter()
//            .enumerate()
//            .map(|(i, pep)| {
//                Protein::new(i + self.protein_set.len(), pep.0.clone(), pep.1.clone(), true)
//            }).collect();
//        let mut protein_set_copy = self.protein_set.clone();
//        protein_set_copy.add_new(new_proteins);
//
//        let mut hsps = HashSet::new();
//        self.hsps.iter().for_each(|hsp| { hsps.insert(hsp.clone()); });
//        peptide_hsps.iter().for_each(|hsp| { hsps.insert(hsp.clone()); });
//
//        Ok(score_interactions(&protein_set_copy, hsps, self.training_pairs.clone(), kmer_size, self.rank, self.world_size, false))
//    }
//}
>>>>>>> 8230d32737e34e598414dac9fca0330f9fba71da
