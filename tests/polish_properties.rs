use proptest::prelude::*;
use raptor::graph::polish::{polish_contig, polish_contig_parallel, polish_contig_string};
use raptor::io::fastq::FastqRecord;
use rayon::ThreadPoolBuilder;

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

fn mutate_read(template: &str, mutation: Option<(usize, u8)>) -> String {
    let mut bytes = template.as_bytes().to_vec();
    if let Some((position, base_selector)) = mutation {
        if !bytes.is_empty() {
            let idx = position % bytes.len();
            let mut replacement = BASES[(base_selector as usize) % BASES.len()];
            if replacement == bytes[idx] {
                replacement = BASES[((base_selector as usize) + 1) % BASES.len()];
            }
            bytes[idx] = replacement;
        }
    }
    String::from_utf8(bytes).expect("synthetic DNA read should stay valid UTF-8")
}

fn synth_fastq_reads(
    template: &str,
    num_reads: usize,
    mutations: &[(usize, u8)],
) -> Vec<FastqRecord> {
    (0..num_reads)
        .map(|idx| {
            let mutation = if mutations.is_empty() {
                None
            } else {
                Some(mutations[idx % mutations.len()])
            };
            let sequence = mutate_read(template, mutation);
            FastqRecord {
                header: format!("@read{}", idx),
                sequence: sequence.clone(),
                plus: "+".to_string(),
                quality: "I".repeat(sequence.len()),
            }
        })
        .collect()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(32))]

    #[test]
    fn parallel_polish_matches_sequential_and_is_thread_count_deterministic(
        draft in "[ACGT]{24,64}",
        window in 3usize..9,
        num_reads in 3usize..8,
        chunk_divisor in 2usize..6,
        mutations in prop::collection::vec((0usize..96, 0u8..4u8), 0..8),
    ) {
        let reads = synth_fastq_reads(&draft, num_reads, &mutations);

        let mut chunk_size = (draft.len() / chunk_divisor.max(1)).max(1);
        if chunk_size * 2 > draft.len() {
            chunk_size = (draft.len() / 2).max(1);
        }
        let window = window.min(draft.len());

        let sequential = polish_contig(&draft, &reads, window);

        let one_thread = ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .expect("build single-thread rayon pool")
            .install(|| polish_contig_parallel(&draft, &reads, window, chunk_size));
        let four_threads = ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("build multi-thread rayon pool")
            .install(|| polish_contig_parallel(&draft, &reads, window, chunk_size));

        prop_assert_eq!(one_thread.as_str(), sequential.as_str());
        prop_assert_eq!(four_threads.as_str(), sequential.as_str());
    }

    #[test]
    fn polish_contig_string_is_invariant_to_read_order(
        draft in "[ACGT]{16,48}",
        num_reads in 3usize..9,
        mutations in prop::collection::vec((0usize..96, 0u8..4u8), 1..8),
    ) {
        let reads = synth_fastq_reads(&draft, num_reads, &mutations);
        let forward_reads: Vec<String> = reads.iter().map(|r| r.sequence.clone()).collect();
        let mut reverse_reads = forward_reads.clone();
        reverse_reads.reverse();

        let forward = polish_contig_string(&draft, &forward_reads, 0.4);
        let reverse = polish_contig_string(&draft, &reverse_reads, 0.4);
        prop_assert_eq!(forward, reverse);
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(16))]

    #[test]
    fn parallel_polish_handles_arbitrary_utf8_inputs_deterministically(
        draft_bytes in prop::collection::vec(any::<u8>(), 1..64),
        read_bytes in prop::collection::vec(any::<u8>(), 1..64),
        num_reads in 3usize..8,
        window in 1usize..8,
        chunk_divisor in 2usize..6,
    ) {
        let draft = String::from_utf8_lossy(&draft_bytes).into_owned();
        let read_template = String::from_utf8_lossy(&read_bytes).into_owned();

        let reads: Vec<FastqRecord> = (0..num_reads)
            .map(|idx| FastqRecord {
                header: format!("@utf8_read{}", idx),
                sequence: read_template.clone(),
                plus: "+".to_string(),
                quality: "I".repeat(read_template.len()),
            })
            .collect();

        let mut chunk_size = (draft.len() / chunk_divisor.max(1)).max(1);
        if chunk_size * 2 > draft.len() {
            chunk_size = (draft.len() / 2).max(1);
        }
        let window = window.min(draft.len());

        let sequential = polish_contig(&draft, &reads, window);

        let one_thread = ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .expect("build single-thread rayon pool")
            .install(|| polish_contig_parallel(&draft, &reads, window, chunk_size));
        let four_threads = ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("build multi-thread rayon pool")
            .install(|| polish_contig_parallel(&draft, &reads, window, chunk_size));

        prop_assert_eq!(one_thread.as_str(), sequential.as_str());
        prop_assert_eq!(four_threads.as_str(), sequential.as_str());
        prop_assert_eq!(sequential.len(), draft.len());
    }
}
