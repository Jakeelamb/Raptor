use crate::io::fastq::FastqRecord;

pub fn polish_sequence(seq: &str, reads: &[FastqRecord], k: usize) -> String {
    let mut counts = vec![[0u32; 4]; seq.len()];
    for read in reads {
        if let Some(pos) = seq.find(&read.sequence) {
            for (i, c) in read.sequence.bytes().enumerate() {
                let idx = match c {
                    b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => continue,
                };
                counts[pos + i][idx] += 1;
            }
        }
    }

    counts.iter().map(|arr| {
        match arr.iter().enumerate().max_by_key(|(_, &v)| v) {
            Some((0, _)) => 'A',
            Some((1, _)) => 'C',
            Some((2, _)) => 'G',
            Some((3, _)) => 'T',
            _ => 'N',
        }
    }).collect()
} 