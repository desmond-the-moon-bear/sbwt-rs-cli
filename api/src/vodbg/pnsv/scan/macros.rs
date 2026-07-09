
macro_rules! define_lcs {
    ($name:ident, $word:ty, $element:ty) => {
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct $name {
            pub words: Vec<$word>,
            pub n: usize,
        }
    };
}

pub(super) use define_lcs;

