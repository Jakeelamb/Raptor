# Raptor Large Genome Assembler - Development Plan

**Status: IN PROGRESS**
**Last Updated: 2026-01-29**
**Next Session: Start with Task 1 - CLI Integration**

---

## Overview

Raptor is an RNA-Seq/genome assembler being extended to handle extremely large genomes (100+ Gb) on personal computers with limited RAM (16 GB). The disk-based k-mer counting and graph cleaning infrastructure is complete. This plan outlines the remaining work.

## Completed Features

- [x] Disk-based k-mer counting (`src/kmer/disk_counting_v2.rs`)
- [x] 2-bit reversible k-mer encoding (not hashes)
- [x] Error correction (singleton k-mer merging)
- [x] Tip removal (dead-end path cleaning)
- [x] Bubble popping (alternative path removal)
- [x] Basic contig assembly with bidirectional extension
- [x] All 85 tests passing

## Pending Tasks

### Task 1: CLI Integration (START HERE)
**Priority: HIGH** | **Complexity: Low**

Add the large genome assembler to the main CLI so users can run:
```bash
raptor assemble-large --input reads.fq --output contigs.fa --kmer 31 --min-count 2
```

**Files to modify:**
- `src/main.rs` - Add new subcommand
- `src/cli/` - Add argument parsing for large genome options
- Consider adding `src/bin/assemble_large.rs` as standalone binary

**Key parameters to expose:**
- `--input` / `-i` - Input FASTQ file(s)
- `--output` / `-o` - Output FASTA file
- `--kmer` / `-k` - K-mer size (default: 31)
- `--min-count` / `-c` - Minimum k-mer count (default: 2)
- `--min-contig` - Minimum contig length (default: 200)
- `--threads` / `-t` - Number of threads
- `--temp-dir` - Temporary directory for disk buckets
- `--max-tip-len` - Maximum tip length to remove (default: 100)
- `--max-bubble-len` - Maximum bubble length to pop (default: 50)

---

### Task 2: Paired-End Read Support
**Priority: HIGH** | **Complexity: Medium**

Currently only handles single-end reads. Need to:
1. Parse paired FASTQ files (R1/R2)
2. Track insert size distribution
3. Use pair information during graph traversal

**Files to create/modify:**
- `src/io/fastq.rs` - Add paired-end streaming
- `src/pipeline/large_genome_assembler.rs` - Use pair constraints

---

### Task 3: Scaffolding
**Priority: HIGH** | **Complexity: Medium-High**

Link contigs into scaffolds using paired-end information:
1. Align reads to assembled contigs
2. Identify contig pairs connected by read pairs
3. Estimate gap sizes from insert size distribution
4. Build scaffold graph and extract paths

**Files to create:**
- `src/pipeline/scaffolder.rs`
- `src/graph/scaffold_graph.rs`

**Algorithm outline:**
```
1. Index contigs (minimizer-based)
2. Map paired reads to contigs
3. For each read pair spanning two contigs:
   - Record (contig_A, contig_B, orientation, gap_estimate)
4. Build scaffold graph with edge weights = number of supporting pairs
5. Filter edges with < N supporting pairs
6. Extract scaffold paths (longest paths through graph)
7. Output scaffolds with N's for gaps
```

---

### Task 4: Contig Polishing
**Priority: MEDIUM** | **Complexity: Medium**

Correct errors in assembled contigs using raw reads:
1. Align reads back to contigs
2. Build pileup at each position
3. Call consensus based on coverage

**Files to create:**
- `src/pipeline/polisher.rs`

**Existing code to leverage:**
- `src/graph/polish.rs` - Has some polishing logic already

---

### Task 5: Repeat Resolution
**Priority: MEDIUM** | **Complexity: High**

Better handling of repetitive regions:
1. Identify high-multiplicity k-mers (repeats)
2. Use coverage depth to estimate copy number
3. Resolve repeat boundaries using unique flanking k-mers

**Approach:**
- Track k-mer counts during assembly
- At branch points, use coverage ratios to guide traversal
- Consider "repeat graph" representation

---

### Task 6: Long Read Integration (Future)
**Priority: LOW** | **Complexity: High**

Hybrid assembly with PacBio/Nanopore:
1. Use short reads for accurate k-mer graph
2. Use long reads to span repeats
3. Anchor long reads to short-read contigs

---

### Task 7: Performance Optimization
**Priority: LOW** | **Complexity: Medium**

- Parallel bucket processing (currently sequential)
- Memory-mapped I/O for bucket files
- SIMD k-mer encoding/decoding
- Compressed bucket storage (LZ4)

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

## Notes for Next Session

1. **Start with Task 1 (CLI Integration)** - This is the quickest win and makes all the new features accessible to users.

2. The `LargeGenomeConfig` struct in `src/pipeline/large_genome_assembler.rs` already has all the parameters - just need to wire them to CLI arguments.

3. Look at existing CLI structure in `src/main.rs` for patterns to follow.

4. After CLI integration, move to Task 2 (paired-end support) since most real data is paired.

---

## Memory Budget (Reference)

For 120 Gb genome on 16 GB RAM machine:
- Bucket size: ~2 GB per bucket
- Adjacency cache: ~1-2 GB
- Working memory: ~1-2 GB
- Target: < 4 GB peak usage
