/// An Iterator that returns all triples (i, j, k) where 0 <= i < j < k < n for fixed n.
pub struct UnorderedTriples {
    n: usize,
    i: usize,
    j: usize,
    k: usize
}

impl UnorderedTriples {
    pub fn new(n: usize) -> Self {
        Self { n, i: 0, j: 1, k: 2 }
    }
}

impl Iterator for UnorderedTriples {
    type Item = (usize, usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        if self.n <= 2 {
            return None;
        }
        let result = (self.i, self.j, self.k);
        if self.i >= self.n - 2 {
            return None;
        }
        self.k += 1;
        if self.k >= self.n {
            self.j += 1;
            // Reset k.
            self.k = self.j + 1;
            if self.j >= self.n - 1 {
                self.i += 1;
                // Reset j and k.
                self.j = self.i + 1;
                self.k = self.j + 1;
            }
        }
        Some(result)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_iter_unordered_triples_1() {
        let result = UnorderedTriples::new(5).collect::<Vec<_>>();
        assert_eq!(
            result,
            vec![
                (0, 1, 2),
                (0, 1, 3),
                (0, 1, 4),
                (0, 2, 3),
                (0, 2, 4),
                (0, 3, 4),
                (1, 2, 3),
                (1, 2, 4),
                (1, 3, 4),
                (2, 3, 4),
            ]
        );
    }

    #[test]
    fn test_iter_unordered_triples_2() {
        let n = 10;
        let n_choose_3 = 120;
        let result = UnorderedTriples::new(n).collect::<Vec<_>>();
        assert_eq!(result.len(), n_choose_3);
        assert!(result.into_iter().all(|(i, j, k)| {
            i < j && j < k && k < n
        }));
    }
}