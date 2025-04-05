use raptor::graph::polish::polish_contig;
use raptor::io::fastq::FastqRecord;

#[test]
fn test_polish_simple() {
    let reads = vec![
        FastqRecord {
            header: "@read1".into(),
            sequence: "ACGTACGT".into(),
            plus: "+".into(),
            quality: "FFFFFFFF".into(),
        },
        FastqRecord {
            header: "@read2".into(),
            sequence: "ACGTACGT".into(),
            plus: "+".into(),
            quality: "FFFFFFFF".into(),
        },
    ];

    // The T at the end should be corrected to T since both reads have T there
    let result = polish_contig("ACGTACGA", &reads, 4);
    assert_eq!(result, "ACGTACGT");
}

#[test]
fn test_polish_consensus() {
    let reads = vec![
        FastqRecord {
            header: "@read1".into(),
            sequence: "AAAAAAA".into(),
            plus: "+".into(),
            quality: "FFFFFFF".into(),
        },
        FastqRecord {
            header: "@read2".into(),
            sequence: "AAAAAAA".into(),
            plus: "+".into(),
            quality: "FFFFFFF".into(),
        },
        FastqRecord {
            header: "@read3".into(),
            sequence: "CCCCCCC".into(),
            plus: "+".into(),
            quality: "FFFFFFF".into(),
        },
    ];

    // The consensus should be A (majority) even though there's a C in one read
    let result = polish_contig("AAAAAAA", &reads, 3);
    assert_eq!(result, "AAAAAAA");
}

#[test]
fn test_polish_multiple_errors() {
    let reads = vec![
        FastqRecord {
            header: "@read1".into(),
            sequence: "ACGTACGTACGT".into(),
            plus: "+".into(),
            quality: "FFFFFFFFFFFF".into(),
        },
        FastqRecord {
            header: "@read2".into(),
            sequence: "ACGTACGTACGT".into(),
            plus: "+".into(),
            quality: "FFFFFFFFFFFF".into(),
        },
    ];

    // Multiple errors should be corrected
    let result = polish_contig("ACTTACTTACTT", &reads, 5);
    assert_eq!(result, "ACGTACGTACGT");
} 