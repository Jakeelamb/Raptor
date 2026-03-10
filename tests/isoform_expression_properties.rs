use proptest::prelude::*;
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use raptor::io::isoform_runner::load_expression_data;
use std::collections::BTreeMap;
use std::io::Write;
use tempfile::NamedTempFile;

fn write_expression_file(rows: &[(usize, f64)]) -> NamedTempFile {
    let mut file = NamedTempFile::new().expect("create temp expression file");
    writeln!(file, "# synthetic expression table").expect("write comment row");
    writeln!(file).expect("write blank row");
    for (id, coverage) in rows {
        writeln!(file, "{id}\t{coverage:.6}").expect("write expression row");
    }
    file.flush().expect("flush expression file");
    file
}

proptest! {
    #[test]
    fn load_expression_data_is_invariant_to_unique_row_permutation(
        rows in prop::collection::vec((0usize..500usize, 0u32..1_000_000u32), 0..64)
    ) {
        let mut dedup = BTreeMap::new();
        for (id, coverage_scaled) in rows {
            dedup.entry(id).or_insert(coverage_scaled);
        }

        let canonical: Vec<(usize, f64)> = dedup
            .iter()
            .map(|(&id, &coverage_scaled)| (id, coverage_scaled as f64 / 1000.0))
            .collect();
        let mut permuted = canonical.clone();
        let mut rng = StdRng::seed_from_u64(0x15_0F_0F_12_34_56_78_AB);
        permuted.shuffle(&mut rng);

        let file_a = write_expression_file(&canonical);
        let file_b = write_expression_file(&permuted);

        let parsed_a = load_expression_data(file_a.path().to_str().expect("utf8 path"))
            .expect("canonical expression data should parse");
        let parsed_b = load_expression_data(file_b.path().to_str().expect("utf8 path"))
            .expect("permuted expression data should parse");

        prop_assert_eq!(parsed_a, parsed_b);
    }
}
