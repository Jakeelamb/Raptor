use assembler::graph::simplify::collapse_repeats;
use assembler::graph::assembler::Contig;

#[test]
fn test_repeat_collapse() {
    let contigs = vec![
        Contig { sequence: "ATATATAT".into(), kmer_path: vec![] },
        Contig { sequence: "ATATATAT".into(), kmer_path: vec![] },
        Contig { sequence: "CGTACGTA".into(), kmer_path: vec![] },
    ];

    // Test with a minimum repeat length of 3
    let collapsed = collapse_repeats(contigs, 3);
    assert_eq!(collapsed.len(), 2);
    
    // Verify the first two contigs were collapsed to one
    let sequences: Vec<_> = collapsed.iter().map(|c| &c.sequence).collect();
    assert!(sequences.contains(&&"ATATATAT".to_string()));
    assert!(sequences.contains(&&"CGTACGTA".to_string()));
}

#[test]
fn test_repeat_collapse_with_different_thresholds() {
    let contigs = vec![
        Contig { sequence: "AAAACCCGGTT".into(), kmer_path: vec![] },
        Contig { sequence: "AAACCCCTTT".into(), kmer_path: vec![] },
        Contig { sequence: "GGGGGCCCCC".into(), kmer_path: vec![] },
    ];
    
    // With high threshold, no collapsing should occur
    let collapsed_high = collapse_repeats(contigs.clone(), 6);
    assert_eq!(collapsed_high.len(), 3);
    
    // With lower threshold, some collapsing may occur based on RLE patterns
    let collapsed_low = collapse_repeats(contigs, 2);
    assert!(collapsed_low.len() <= 3);
} 