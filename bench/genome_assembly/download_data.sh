#!/bin/bash
# Download test datasets for genome assembly benchmarking
# Usage: ./download_data.sh [drosophila|human|all]

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"

mkdir -p "${DATA_DIR}"

download_drosophila() {
    echo "=== Downloading Drosophila melanogaster data ==="
    local DMEL_DIR="${DATA_DIR}/drosophila"
    mkdir -p "${DMEL_DIR}"

    # Download reference genome for evaluation
    echo "Downloading reference genome..."
    if [ ! -f "${DMEL_DIR}/reference.fa.gz" ]; then
        wget -O "${DMEL_DIR}/reference.fa.gz" \
            "https://ftp.ensembl.org/pub/release-110/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.46.dna.toplevel.fa.gz"
    fi

    # Download Illumina reads from SRA (SRR6425988 - D. melanogaster whole genome)
    # This is a moderate coverage dataset (~30x)
    echo "Downloading Illumina paired-end reads..."
    if [ ! -f "${DMEL_DIR}/reads_1.fastq.gz" ]; then
        # Using ENA for faster downloads
        wget -O "${DMEL_DIR}/reads_1.fastq.gz" \
            "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR642/008/SRR6425988/SRR6425988_1.fastq.gz" || \
        # Fallback: use fasterq-dump if available
        (command -v fasterq-dump && fasterq-dump --split-files -O "${DMEL_DIR}" SRR6425988 && \
         gzip "${DMEL_DIR}/SRR6425988_1.fastq" && mv "${DMEL_DIR}/SRR6425988_1.fastq.gz" "${DMEL_DIR}/reads_1.fastq.gz")
    fi

    if [ ! -f "${DMEL_DIR}/reads_2.fastq.gz" ]; then
        wget -O "${DMEL_DIR}/reads_2.fastq.gz" \
            "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR642/008/SRR6425988/SRR6425988_2.fastq.gz" || \
        (command -v fasterq-dump && gzip "${DMEL_DIR}/SRR6425988_2.fastq" && \
         mv "${DMEL_DIR}/SRR6425988_2.fastq.gz" "${DMEL_DIR}/reads_2.fastq.gz")
    fi

    echo "Drosophila data downloaded to ${DMEL_DIR}"
    echo "  Reference: ${DMEL_DIR}/reference.fa.gz"
    echo "  Reads R1:  ${DMEL_DIR}/reads_1.fastq.gz"
    echo "  Reads R2:  ${DMEL_DIR}/reads_2.fastq.gz"
}

download_human_chr21() {
    echo "=== Downloading Human Chromosome 21 data ==="
    local HUMAN_DIR="${DATA_DIR}/human_chr21"
    mkdir -p "${HUMAN_DIR}"

    # Download chr21 reference (smaller for testing)
    echo "Downloading chr21 reference..."
    if [ ! -f "${HUMAN_DIR}/chr21.fa.gz" ]; then
        wget -O "${HUMAN_DIR}/chr21.fa.gz" \
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz"
    fi

    # For human reads, we'll use a subset from a public dataset
    # NA12878 is the standard reference sample
    echo "Downloading Human Illumina reads (NA12878 subset)..."
    if [ ! -f "${HUMAN_DIR}/reads_1.fastq.gz" ]; then
        # Download from ENA - ERR194147 is NA12878 Illumina data
        # We'll download a small subset (first 10M reads)
        wget -O "${HUMAN_DIR}/full_reads_1.fastq.gz" \
            "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz" || true

        if [ -f "${HUMAN_DIR}/full_reads_1.fastq.gz" ]; then
            echo "Extracting subset of reads (10M)..."
            zcat "${HUMAN_DIR}/full_reads_1.fastq.gz" | head -n 40000000 | gzip > "${HUMAN_DIR}/reads_1.fastq.gz"
        fi
    fi

    if [ ! -f "${HUMAN_DIR}/reads_2.fastq.gz" ]; then
        wget -O "${HUMAN_DIR}/full_reads_2.fastq.gz" \
            "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz" || true

        if [ -f "${HUMAN_DIR}/full_reads_2.fastq.gz" ]; then
            zcat "${HUMAN_DIR}/full_reads_2.fastq.gz" | head -n 40000000 | gzip > "${HUMAN_DIR}/reads_2.fastq.gz"
        fi
    fi

    echo "Human chr21 data downloaded to ${HUMAN_DIR}"
}

download_simulated() {
    echo "=== Generating simulated test data ==="
    local SIM_DIR="${DATA_DIR}/simulated"
    mkdir -p "${SIM_DIR}"

    # Create a simple simulated dataset for quick testing
    SIM_DIR="${SIM_DIR}" python3 << 'EOF'
import random
import gzip
import os

output_dir = os.environ.get('SIM_DIR', 'data/simulated')
os.makedirs(output_dir, exist_ok=True)

# Generate a random reference genome (5 Mb)
print("Generating reference genome...")
genome_size = 5_000_000
bases = ['A', 'C', 'G', 'T']
genome = ''.join(random.choices(bases, k=genome_size))

with open(f'{output_dir}/reference.fa', 'w') as f:
    f.write('>simulated_genome\n')
    for i in range(0, len(genome), 80):
        f.write(genome[i:i+80] + '\n')

# Generate paired-end reads (150bp, 30x coverage)
print("Generating reads...")
read_len = 150
insert_size = 350
insert_std = 50
coverage = 30
num_pairs = (genome_size * coverage) // (2 * read_len)

with gzip.open(f'{output_dir}/reads_1.fastq.gz', 'wt') as r1, \
     gzip.open(f'{output_dir}/reads_2.fastq.gz', 'wt') as r2:

    for i in range(num_pairs):
        # Random insert position
        insert = int(random.gauss(insert_size, insert_std))
        insert = max(2 * read_len, min(insert, 600))

        pos = random.randint(0, genome_size - insert)

        # Read 1 (forward)
        seq1 = genome[pos:pos + read_len]
        qual1 = 'I' * read_len

        # Read 2 (reverse complement)
        seq2_raw = genome[pos + insert - read_len:pos + insert]
        seq2 = seq2_raw.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        qual2 = 'I' * read_len

        # Add some sequencing errors (1%)
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

        if i % 100000 == 0:
            print(f'  Generated {i}/{num_pairs} read pairs...')

print(f"Simulated data generated in {output_dir}")
EOF
}

# Main
case "${1:-all}" in
    drosophila)
        download_drosophila
        ;;
    human)
        download_human_chr21
        ;;
    simulated)
        download_simulated
        ;;
    all)
        download_simulated
        download_drosophila
        download_human_chr21
        ;;
    *)
        echo "Usage: $0 [drosophila|human|simulated|all]"
        exit 1
        ;;
esac

echo ""
echo "=== Download complete ==="
echo "Data directory: ${DATA_DIR}"
ls -la "${DATA_DIR}"
