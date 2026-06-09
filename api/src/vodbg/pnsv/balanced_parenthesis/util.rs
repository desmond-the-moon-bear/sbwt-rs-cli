pub fn bit_width(value: usize) -> usize {
    64 - u64::leading_zeros(value as u64) as usize
}

pub fn one_indices_iterator<B: Iterator<Item = bool> + Clone>(iterator: B) -> impl Iterator<Item = usize> + Clone {
    iterator.enumerate().filter(|(_, value)| *value).map(|(index, _)| index)
}
