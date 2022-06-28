use pyo3::prelude::*;

pub mod protein;
pub mod proteinset;
pub mod constants;
pub mod fileio;
pub mod extraction;
pub mod seed;
pub mod scoring;
pub mod smer;
pub mod hsp;
pub mod location;
pub mod similarity;
pub mod utils;
pub mod processing;
pub mod prediction;

pub mod pymodules;

//#[pymodule]
//fn darwin(_py: Python, m: &PyModule) -> PyResult<()> {
//    m.add_class::<pymodules::SPRINT>()?;
//    m.add_class::<hsp::HSP>()?;
//    Ok(())
//}
