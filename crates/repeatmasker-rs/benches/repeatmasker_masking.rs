use criterion::{black_box, criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use repeatmasker_rs::{mask_fasta_records, FastaRecord, MaskInterval, MaskMode};
use std::collections::BTreeMap;

fn synthetic_records(record_count: usize, sequence_len: usize) -> Vec<FastaRecord> {
    let template = b"ACGT".repeat(sequence_len / 4 + 1);
    (0..record_count)
        .map(|idx| FastaRecord {
            id: format!("seq{idx}"),
            description: String::new(),
            sequence: template[..sequence_len].to_vec(),
        })
        .collect()
}

fn synthetic_intervals(
    record_count: usize,
    sequence_len: usize,
    intervals_per_record: usize,
) -> BTreeMap<String, Vec<MaskInterval>> {
    let mut intervals = BTreeMap::new();
    let span = 48usize;
    let stride = (sequence_len / intervals_per_record.max(1)).max(span + 1);
    for idx in 0..record_count {
        let seq_id = format!("seq{idx}");
        let mut seq_intervals = Vec::with_capacity(intervals_per_record);
        let mut start = 1usize;
        for _ in 0..intervals_per_record {
            let end = (start + span - 1).min(sequence_len);
            seq_intervals.push(MaskInterval { begin: start, end });
            start = (start + stride).min(sequence_len);
        }
        intervals.insert(seq_id, seq_intervals);
    }
    intervals
}

fn bench_repeatmasker_masking(c: &mut Criterion) {
    let mut group = c.benchmark_group("repeatmasker_masking");
    for (record_count, sequence_len, intervals_per_record) in
        [(8usize, 250_000usize, 256usize), (16, 250_000, 512)]
    {
        let records = synthetic_records(record_count, sequence_len);
        let intervals = synthetic_intervals(record_count, sequence_len, intervals_per_record);
        group.bench_with_input(
            BenchmarkId::new(
                format!("records_{record_count}"),
                format!("len_{sequence_len}_ints_{intervals_per_record}"),
            ),
            &(records, intervals),
            |b, (records, intervals)| {
                b.iter_batched(
                    || records.clone(),
                    |mut cloned_records| {
                        let stats = mask_fasta_records(&mut cloned_records, intervals, MaskMode::N);
                        black_box((cloned_records[0].sequence[0], stats.masked_bases))
                    },
                    BatchSize::LargeInput,
                );
            },
        );
    }
    group.finish();
}

criterion_group!(benches, bench_repeatmasker_masking);
criterion_main!(benches);
