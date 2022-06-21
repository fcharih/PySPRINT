use pyo3::prelude::*;

use crate::rsprint::proteinset::ProteinSet;
use crate::rsprint::location::Location;

#[pyclass]
#[derive(Hash)]
pub struct HSP {
    location1: Location,
    location2: Location,
    length: usize
}

#[pymethods]
impl HSP {
    #[staticmethod]
    pub fn from_tuple(tuple: (usize, usize, usize, usize, usize)) -> PyResult<HSP> {
        let location1 = Location::new(tuple.0, tuple.1);
        let location2 = Location::new(tuple.2, tuple.3);
        Ok(HSP::new(location1, location2, tuple.4))
    }

    pub fn to_tuple(&self) -> PyResult<(usize, usize, usize, usize, usize)> {
        Ok((self.location1.index(), self.location1.position(), self.location2.index(), self.location2.position(), self.length))
    }
}

impl HSP {
    pub fn new(location1: Location, location2: Location, length: usize) -> Self {

        if location1.index() == location2.index() {
            let min = std::cmp::min(location1.position(), location2.position());
            let max = std::cmp::max(location1.position(), location2.position());
            return HSP {
                location1: Location::new(location1.index(), min),
                location2: Location::new(location2.index(), max),
                length
            }
        }

        match location1.index() < location2.index() {
            true => HSP {
                location1,
                location2,
                length
            },
            false => HSP {
                location2: location1,
                location1: location2,
                length
            }
        }
        
    }
    
    pub fn len(&self) -> usize {
        self.length
    }

    pub fn location(&self, index: usize) -> &Location {
        if index == 0 {
            return &self.location1;
        } else if index == 1 {
            return &self.location2;
        } else {
            panic!("Attempted to get bad location for an HSP.");
        }
    }

    pub fn as_string(&self, protein_set: &ProteinSet) -> String {
        let name1 = protein_set.get_protein_by_id(self.location1.index()).name();
        let name2 = protein_set.get_protein_by_id(self.location2.index()).name();
        if name1 == name2 {
            let min = std::cmp::min(self.location1.position(), self.location2.position());
            let max = std::cmp::max(self.location1.position(), self.location2.position());
            return format!("{} {} {} {} {}", name1, name2, min, max, self.length);
        } else if name1 < name2 {
            return format!("{} {} {} {} {}", name1, name2, self.location1.position(), self.location2.position(), self.length);
        } else {
            return format!("{} {} {} {} {}", name2, name1, self.location2.position(), self.location1.position(), self.length);
        }

    }
}


impl HSP {
    pub fn from(string: String, protein_set: &ProteinSet) -> Self {
        let tokens: Vec<&str> = string.split_whitespace().collect();
        let protein1 = tokens[0].to_owned();
        let protein2 = tokens[1].to_owned();
        let position1 = tokens[2].parse::<usize>().unwrap();
        let position2 = tokens[3].parse::<usize>().unwrap();
        let length = tokens[4].parse::<usize>().unwrap();
        match protein1 < protein2 {
            true => HSP {
                location1: Location::new(protein_set.get_protein_by_name(&protein1).index(), position1),
                location2: Location::new(protein_set.get_protein_by_name(&protein2).index(), position2),
                length
            },
            false => HSP {
                location2: Location::new(protein_set.get_protein_by_name(&protein1).index(), position1),
                location1: Location::new(protein_set.get_protein_by_name(&protein2).index(), position2),
                length
            }
        }
        
    }
}

impl Clone for HSP {
    fn clone(&self) -> HSP {
        HSP {
            location1: self.location1.clone(),
            location2: self.location2.clone(),
            length: self.length
        }
    }
}

impl PartialEq for HSP {
    fn eq(&self, other: &Self) -> bool {
        self.location1 == other.location1 &&
        self.location2 == other.location2 &&
        self.length == other.length
    }
}

impl Eq for HSP {}

/* impl Ord for HSP {
    fn cmp(&self, other: &HSP) -> std::cmp::Ordering {
        (&self.location1.0, &self.location2.0,
         self.location1.1, self.location2.1)
            .cmp(&(&other.location1.0, &other.location2.0,
                 other.location1.1, other.location2.1))
    }
}

impl PartialOrd for HSP {
    fn partial_cmp(&self, other: &HSP) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

 */
