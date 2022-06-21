pub struct Seed {
    seed: String,
    value: u64,
    length: usize,
    non_zero_positions: Vec<usize>
}

impl Seed {
    pub fn new(seed: &str) -> Seed {
        let mut int_seed: u64 = 0;
        let mut non_zero_positions = vec![];
        for (i, character) in seed.chars().enumerate() {
            if character == '*' {
                int_seed = int_seed << 5; // just shift
            } else if character == '1' {
                int_seed = int_seed << 5;
                int_seed = int_seed + 31; // 31 is decimal for 11111
                non_zero_positions.push(i);
            } else {
                panic!("Encountered invalid character {} while converting the seed.", character);
            }
        }

        Seed {
            seed: seed.to_string(),
            value: int_seed,
            length: seed.len(),
            non_zero_positions
        }
    }

    pub fn len(&self) -> usize {
        self.seed.len()
    }

    pub fn value(&self) -> u64 {
        self.value
    }
    
    pub fn clone(&self) -> Seed {
        Seed {
            seed: self.seed.clone(),
            value: self.value,
            length: self.length,
            non_zero_positions: self.non_zero_positions.clone()
        }
    }

    pub fn as_string(&self) -> String {
        self.seed.clone()
    }
}
