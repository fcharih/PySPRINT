use std::collections::HashSet;

use numpy::ToPyArray;
use numpy::{PyArray2};
use pyo3::prelude::*;
use pyo3::pymodule;

use crate::rsprint::hsp::HSP;
use crate::rsprint::prediction::score_interactions;

use super::{proteinset::ProteinSet, protein::Protein, extraction::extract_hsps};

#[pymodule]
fn rsprint(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    #[pyfn(m)]
    #[pyo3(name = "extract_peptide_hsps")]
    pub fn extract_peptide_hsps_py(
        proteins: Vec<(String, String)>,
        peptides: Vec<(String, String)>,
        process_rank:usize,
        world_size: usize,
        t_smer: i16,
        t_hsp: i16,
        kmer_size: usize
    ) -> PyResult<HashSet<(usize, usize, usize, usize, usize)>> {
        
        let mut protein_set = ProteinSet::new(
            convert_tuples_to_proteins(proteins, false)
        );

        protein_set.add_new(
            convert_tuples_to_proteins(peptides, true),
            true
        );

        let hsps = extract_hsps(
            &protein_set,
            kmer_size,
            t_smer,
            t_hsp,
            process_rank,
            world_size,
            true,
            true,
            false
        );

        Ok(hsps.iter().map(|hsp| hsp.to_tuple().unwrap()).collect())
    }

    #[pyfn(m)]
    #[pyo3(name = "score_peptides")]
    pub fn score_interactions_py<'py>(
        py: Python<'py>,
        proteins: Vec<(String, String)>,
        peptides: Vec<(String, String)>,
        hsps: HashSet<(usize, usize, usize, usize, usize)>,
        training_pairs: Vec<(String, String)>,
        kmer_size: usize,
        process_rank:usize,
        world_size: usize
    ) -> &'py PyArray2<f32> {
        
        let mut protein_set = ProteinSet::new(
            convert_tuples_to_proteins(proteins, false)
        );

        protein_set.add_new(
            convert_tuples_to_proteins(peptides, true),
            true
        );

        let parsed_hsps: HashSet<HSP> = hsps
            .iter()
            .map(|h| HSP::from_tuple(*h).unwrap())
            .collect();

        let matrix = score_interactions(
            &protein_set,
            &parsed_hsps,
            &training_pairs,
            kmer_size,
            process_rank,
            world_size,
            false
        );

        matrix.to_pyarray(py)

    }

    Ok(())
}


fn convert_tuples_to_proteins(tuples: Vec<(String, String)>, new: bool) -> Vec<Protein> {
    tuples
        .into_iter() // TODO is it necessary to yield ownership?
        .enumerate()
        .map(|(i, p)| 
            Protein::new(i, p.0, p.1, new)
        )
        .collect::<Vec<Protein>>()
}