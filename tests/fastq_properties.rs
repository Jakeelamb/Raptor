use proptest::prelude::*;
use raptor::io::fastq::{stream_fastq_records_checked, stream_paired_fastq_records_checked};
use std::io::Cursor;

fn render_record(id: &str, seq: &str, qual_len: usize) -> String {
    format!("@{}\n{}\n+\n{}\n", id, seq, "I".repeat(qual_len))
}

proptest! {
    #[test]
    fn checked_fastq_parser_accepts_well_formed_records(
        seq1 in "[ACGTN]{1,128}",
        seq2 in "[ACGTN]{1,128}",
    ) {
        let content = format!(
            "{}{}",
            render_record("read1", &seq1, seq1.len()),
            render_record("read2", &seq2, seq2.len())
        );

        let mut iter = stream_fastq_records_checked(Cursor::new(content.into_bytes()));

        let rec1 = iter.next().expect("missing first record").expect("first record should parse");
        let rec2 = iter.next().expect("missing second record").expect("second record should parse");

        prop_assert_eq!(rec1.sequence, seq1);
        prop_assert_eq!(rec2.sequence, seq2);
        prop_assert!(iter.next().is_none());
    }

    #[test]
    fn checked_fastq_parser_rejects_quality_length_mismatch(
        seq in "[ACGTN]{2,128}",
    ) {
        let shorter = render_record("bad", &seq, seq.len() - 1);
        let mut iter = stream_fastq_records_checked(Cursor::new(shorter.into_bytes()));
        let err = iter
            .next()
            .expect("expected parse result")
            .expect_err("expected invalid-data mismatch");

        prop_assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    }

    #[test]
    fn checked_paired_parser_rejects_mismatched_record_counts(
        seq in "[ACGTN]{8,64}",
    ) {
        let r1 = format!(
            "{}{}",
            render_record("pair/1_a", &seq, seq.len()),
            render_record("pair/1_b", &seq, seq.len())
        );
        let r2 = render_record("pair/2_a", &seq, seq.len());

        let mut iter = stream_paired_fastq_records_checked(
            Cursor::new(r1.into_bytes()),
            Cursor::new(r2.into_bytes()),
        );

        let first = iter.next().expect("first pair missing").expect("first pair should parse");
        prop_assert_eq!(first.0.sequence.as_str(), seq.as_str());
        prop_assert_eq!(first.1.sequence.as_str(), seq.as_str());

        let err = iter
            .next()
            .expect("expected mismatch error")
            .expect_err("expected paired-count error");
        prop_assert_eq!(err.kind(), std::io::ErrorKind::InvalidData);
    }
}
