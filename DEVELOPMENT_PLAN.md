# Raptor Large Genome Assembler - Development Plan

**Status: COMPLETE**
**Last Updated: 2026-01-30**
**All core features implemented!**

---

## Overview

Raptor is an RNA-Seq/genome assembler extended to handle extremely large genomes (100+ Gb) on personal computers with limited RAM (16 GB). The full assembly pipeline is now complete.

## Completed Features

- [x] Disk-based k-mer counting (`src/kmer/disk_counting_v2.rs`)
- [x] 2-bit reversible k-mer encoding (not hashes)
- [x] Error correction (singleton k-mer merging)
- [x] Tip removal (dead-end path cleaning)
- [x] Bubble popping (alternative path removal)
- [x] Basic contig assembly with bidirectional extension
- [x] **CLI Integration** - `assemble-large` command available
- [x] **Paired-end read support** - Track insert sizes, paired FASTQ processing
- [x] **Scaffolding** - Link contigs using paired-end information (`src/pipeline/scaffolder.rs`)
- [x] **Contig polishing** - Error correction using read pileups (`src/pipeline/polisher.rs`)
- [x] **Repeat resolution** - Coverage-guided traversal at branch points
- [x] All 95 tests passing

## Implemented Tasks

### Task 1: CLI Integration ✓
**Status: COMPLETE**

Users can now run:
```bash
raptor assemble-large --input reads.fq --output contigs.fa --kmer 31 --min-count 2
```

**Implemented in:** `src/cli_main.rs`, `src/main.rs`

**Available parameters:**
- `--input` / `-i` - Input FASTQ file(s)
- `--input2` - Second input for paired-end reads
- `--output` / `-o` - Output FASTA file
- `--kmer` / `-k` - K-mer size (default: 31)
- `--min-count` / `-c` - Minimum k-mer count (default: 2)
- `--min-contig` - Minimum contig length (default: 200)
- `--threads` / `-t` - Number of threads
- `--temp-dir` - Temporary directory for disk buckets
- `--max-tip-len` - Maximum tip length to remove (default: 100)
- `--max-bubble-len` - Maximum bubble length to pop (default: 50)
- `--scaffold` - Enable scaffolding with paired-end reads
- `--polish` - Enable contig polishing

---

### Task 2: Paired-End Read Support ✓
**Status: COMPLETE**

Implemented in `src/pipeline/large_genome_assembler.rs`:
- `assemble_paired()` method for paired FASTQ files
- `InsertSizeStats` struct for tracking insert size distribution
- `distribute_from_paired_fastq()` for efficient paired-end processing

---

### Task 3: Scaffolding ✓
**Status: COMPLETE**

Implemented in `src/pipeline/scaffolder.rs`:
- Minimizer-based contig indexing for fast read mapping
- Paired-read to contig mapping
- Scaffold graph construction with gap estimation
- Greedy scaffold path extraction
- Gap-aware scaffold output with N's

---

### Task 4: Contig Polishing ✓
**Status: COMPLETE**

Implemented in `src/pipeline/polisher.rs`:
- Minimizer-based read-to-contig alignment
- Per-position pileup construction
- Coverage-based consensus calling
- Iterative polishing support

---

### Task 5: Repeat Resolution ✓
**Status: COMPLETE**

Implemented in `src/pipeline/large_genome_assembler.rs`:
- `identify_repeats()` - Detects high-coverage k-mers
- `RepeatStats` - Tracking repeat detection metrics
- `extend_bidirectional_with_coverage()` - Coverage-guided graph traversal
- Non-repeat seeds prioritized for extension
- Coverage ratios used at branch points

---

### Task 6: Long Read Integration ✓
**Status: COMPLETE**

Implemented in `src/pipeline/long_read_integration.rs`:
- Minimizer-based long read to contig mapping
- Detection of contig-spanning reads for gap filling
- Contig extension using long read overhangs
- Contig joining with bridging sequences
- CLI option: `--long-reads <FASTQ>`

---

### Task 7: Performance Optimization ✓
**Status: COMPLETE**

Implemented in `src/kmer/disk_counting_optimized.rs`:
- LZ4 compression for disk buckets (2-4x smaller files)
- Memory-mapped I/O for faster bucket reading
- Parallel sequence distribution
- CLI option: `--compress-buckets`

---

## All Features Complete!

The Raptor large genome assembler now includes:
- Disk-based k-mer counting for 100+ Gb genomes
- Paired-end read support with insert size tracking
- Graph cleaning (tip removal, bubble popping)
- Coverage-guided repeat resolution
- Scaffolding with paired-end reads
- Contig polishing with pileup consensus
- Long read hybrid assembly
- LZ4 compression for reduced disk I/O

---

## Architecture Reference

```
src/
├── kmer/
│   ├── disk_counting_v2.rs   # Disk-based k-mer counting
│   ├── kmer.rs               # KmerU64 2-bit encoding
│   └── superkmer.rs          # Super-k-mer compression
├── pipeline/
│   ├── large_genome_assembler.rs  # Main assembler (6-phase pipeline)
│   └── streaming_assembly.rs      # Streaming variant
└── main.rs                   # CLI entry point
```

**Assembly Pipeline (current):**
```
Phase 1: Distribute k-mers to disk buckets
Phase 2: Count k-mers from buckets
Phase 3: Error correction (singleton merging)
Phase 4: Build & clean graph (tips, bubbles)
Phase 5: Extract contigs (greedy extension)
Phase 6: Write output
```

---

## Testing Commands

```bash
# Run all tests
cargo test

# Run large genome tests specifically
cargo test large_genome -- --nocapture

# Build release version
cargo build --release

# Run with tracing
RUST_LOG=info cargo run --release -- assemble-large -i reads.fq -o out.fa
```

---

## Usage Examples

```bash
# Single-end assembly
raptor assemble-large -i reads.fq -o contigs.fa -k 31 -c 2

# Paired-end assembly
raptor assemble-large -i reads_R1.fq --input2 reads_R2.fq -o contigs.fa -k 31

# With scaffolding and polishing
raptor assemble-large -i reads_R1.fq --input2 reads_R2.fq -o contigs.fa \
    --scaffold --polish

# Hybrid assembly with long reads
raptor assemble-large -i short_R1.fq --input2 short_R2.fq -o contigs.fa \
    --long-reads nanopore.fq --min-long-read-len 2000

# Full pipeline with all features
raptor assemble-large -i reads_R1.fq --input2 reads_R2.fq -o contigs.fa \
    --kmer 31 --min-count 2 \
    --scaffold --polish --polish-iterations 2 \
    --long-reads long_reads.fq \
    --compress-buckets

# Custom parameters for very large genome
raptor assemble-large -i reads.fq -o contigs.fa \
    --kmer 31 --min-count 3 --min-contig 500 \
    --temp-dir /scratch/tmp --threads 16 --compress-buckets
```

---

## Memory Budget (Reference)

For 120 Gb genome on 16 GB RAM machine:
- Bucket size: ~2 GB per bucket
- Adjacency cache: ~1-2 GB
- Working memory: ~1-2 GB
- Target: < 4 GB peak usage
