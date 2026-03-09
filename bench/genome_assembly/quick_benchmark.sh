#!/bin/bash
# Quick benchmark using simulated data - tests Raptor only if SPAdes unavailable
# Usage: ./quick_benchmark.sh [threads]

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data/quick_test"
RESULTS_DIR="${SCRIPT_DIR}/results/quick_test_$(date +%Y%m%d_%H%M%S)"
RAPTOR_BIN="${SCRIPT_DIR}/../../target/release/raptor"
THREADS="${1:-8}"

echo "=== Quick Genome Assembly Benchmark ==="
echo "Threads: ${THREADS}"
echo ""

# Build Raptor if needed
if [ ! -f "${RAPTOR_BIN}" ]; then
    echo "Building Raptor..."
    (cd "${SCRIPT_DIR}/../.." && cargo build --release)
fi

# Generate simulated data
mkdir -p "${DATA_DIR}"
mkdir -p "${RESULTS_DIR}"

echo "Generating simulated genome and reads..."
python3 << 'EOF'
import random
import gzip
import os

random.seed(42)  # Reproducible

data_dir = os.environ.get('DATA_DIR', 'data/quick_test')
os.makedirs(data_dir, exist_ok=True)

# Generate a 2 Mb reference genome
print("Generating 2 Mb reference genome...")
genome_size = 2_000_000
bases = ['A', 'C', 'G', 'T']

# Add some repeats to make it more realistic
repeat_unit = ''.join(random.choices(bases, k=500))
genome_parts = []

for i in range(0, genome_size, 10000):
    chunk_size = min(10000, genome_size - i)
    if random.random() < 0.1:  # 10% chance of repeat
        genome_parts.append(repeat_unit * (chunk_size // 500 + 1))
    else:
        genome_parts.append(''.join(random.choices(bases, k=chunk_size)))

genome = ''.join(genome_parts)[:genome_size]

with open(f'{data_dir}/reference.fa', 'w') as f:
    f.write('>chr1\n')
    for i in range(0, len(genome), 80):
        f.write(genome[i:i+80] + '\n')

# Generate paired-end reads (150bp, 50x coverage)
print("Generating reads (50x coverage)...")
read_len = 150
insert_size = 400
insert_std = 50
coverage = 50
num_pairs = (genome_size * coverage) // (2 * read_len)

with gzip.open(f'{data_dir}/reads_1.fastq.gz', 'wt') as r1, \
     gzip.open(f'{data_dir}/reads_2.fastq.gz', 'wt') as r2:

    for i in range(num_pairs):
        insert = int(random.gauss(insert_size, insert_std))
        insert = max(2 * read_len, min(insert, 600))
        pos = random.randint(0, genome_size - insert)

        seq1 = genome[pos:pos + read_len]
        qual1 = 'I' * read_len

        seq2_raw = genome[pos + insert - read_len:pos + insert]
        seq2 = seq2_raw.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        qual2 = 'I' * read_len

        # Add 1% sequencing errors
        def add_errors(seq, rate=0.01):
            seq = list(seq)
            for j in range(len(seq)):
                if random.random() < rate:
                    seq[j] = random.choice([b for b in bases if b != seq[j]])
            return ''.join(seq)

        seq1 = add_errors(seq1)
        seq2 = add_errors(seq2)

        r1.write(f'@read_{i}/1\n{seq1}\n+\n{qual1}\n')
        r2.write(f'@read_{i}/2\n{seq2}\n+\n{qual2}\n')

print(f"Generated {num_pairs} read pairs")
print(f"Data saved to {data_dir}/")
EOF

export DATA_DIR="${DATA_DIR}"

# Calculate input stats
echo ""
echo "=== Input Statistics ==="
READS1="${DATA_DIR}/reads_1.fastq.gz"
READS2="${DATA_DIR}/reads_2.fastq.gz"
NUM_READS=$(zcat "${READS1}" | wc -l)
NUM_READS=$((NUM_READS / 4))
echo "Read pairs: ${NUM_READS}"
echo "Reference size: $(wc -c < "${DATA_DIR}/reference.fa") bytes"

# Run Raptor
echo ""
echo "=== Running Raptor ==="
RAPTOR_OUT="${RESULTS_DIR}/raptor"
mkdir -p "${RAPTOR_OUT}"

START_TIME=$(date +%s)

"${RAPTOR_BIN}" assemble-large \
    -i "${READS1}" \
    --input2 "${READS2}" \
    -o "${RAPTOR_OUT}/contigs.fa" \
    -t "${THREADS}" \
    -k 31 \
    --min-count 3 \
    --scaffold \
    --polish \
    2>&1 | tee "${RAPTOR_OUT}/log.txt"

END_TIME=$(date +%s)
RAPTOR_TIME=$((END_TIME - START_TIME))

# Extract memory usage (approximate from log if available)
RAPTOR_MEM=$(grep -oP 'Disk used: \K[0-9.]+' "${RAPTOR_OUT}/log.txt" 2>/dev/null || echo "0")
RAPTOR_MEM_MB="N/A"

# Calculate Raptor assembly stats
echo ""
echo "=== Raptor Assembly Statistics ==="
python3 << EOF
def calc_n50(lengths):
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    cumsum = 0
    for l in lengths:
        cumsum += l
        if cumsum >= total / 2:
            return l
    return 0

lengths = []
current = 0
with open("${RAPTOR_OUT}/contigs.fa") as f:
    for line in f:
        if line.startswith('>'):
            if current > 0:
                lengths.append(current)
            current = 0
        else:
            current += len(line.strip())
if current > 0:
    lengths.append(current)

lengths.sort(reverse=True)
n50 = calc_n50(lengths)

print(f"Contigs:      {len(lengths)}")
print(f"Total length: {sum(lengths):,} bp")
print(f"N50:          {n50:,} bp")
print(f"Largest:      {lengths[0]:,} bp")
print(f"Contigs ≥1kb: {sum(1 for l in lengths if l >= 1000)}")
EOF

echo ""
echo "Time:   ${RAPTOR_TIME}s"
echo "Memory: ${RAPTOR_MEM_MB} MB"

# Run SPAdes if available
if command -v spades.py &> /dev/null; then
    echo ""
    echo "=== Running SPAdes ==="
    SPADES_OUT="${RESULTS_DIR}/spades"
    mkdir -p "${SPADES_OUT}"

    START_TIME=$(date +%s)

    spades.py \
        -1 "${READS1}" \
        -2 "${READS2}" \
        -o "${SPADES_OUT}" \
        -t "${THREADS}" \
        --careful \
        2>&1 | tee "${SPADES_OUT}/log.txt"

    END_TIME=$(date +%s)
    SPADES_TIME=$((END_TIME - START_TIME))

    echo ""
    echo "=== SPAdes Assembly Statistics ==="
    python3 << EOF
def calc_n50(lengths):
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    cumsum = 0
    for l in lengths:
        cumsum += l
        if cumsum >= total / 2:
            return l
    return 0

lengths = []
current = 0
with open("${SPADES_OUT}/contigs.fasta") as f:
    for line in f:
        if line.startswith('>'):
            if current > 0:
                lengths.append(current)
            current = 0
        else:
            current += len(line.strip())
if current > 0:
    lengths.append(current)

lengths.sort(reverse=True)
n50 = calc_n50(lengths)

print(f"Contigs:      {len(lengths)}")
print(f"Total length: {sum(lengths):,} bp")
print(f"N50:          {n50:,} bp")
print(f"Largest:      {lengths[0]:,} bp")
print(f"Contigs ≥1kb: {sum(1 for l in lengths if l >= 1000)}")
EOF

    echo ""
    echo "Time: ${SPADES_TIME}s"

    # Summary comparison
    echo ""
    echo "=============================================="
    echo "           BENCHMARK SUMMARY"
    echo "=============================================="
    echo ""
    echo "                Raptor          SPAdes"
    echo "Time:           ${RAPTOR_TIME}s          ${SPADES_TIME}s"
    echo "Memory:         ${RAPTOR_MEM_MB} MB"
    echo ""
else
    echo ""
    echo "SPAdes not available - skipping comparison"
    echo "To install: conda install -c bioconda spades"
fi

echo ""
echo "Results saved to: ${RESULTS_DIR}"
