#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
out_dir="${1:-$root/target/rescue-baseline}"
reads="${out_dir}/synthetic_${READS:-20000}x${READ_LEN:-150}.fastq"
summary="${out_dir}/kmer_cpu_opencl.tsv"
k="${K:-25}"
num_reads="${READS:-20000}"
read_len="${READ_LEN:-150}"
genome_len="${GENOME_LEN:-50000}"

mkdir -p "$out_dir"

python3 - "$reads" "$num_reads" "$read_len" "$genome_len" <<'PY'
import random
import sys

path = sys.argv[1]
num_reads = int(sys.argv[2])
read_len = int(sys.argv[3])
genome_len = int(sys.argv[4])

rng = random.Random(7)
bases = "ACGT"
genome = "".join(rng.choice(bases) for _ in range(genome_len))

with open(path, "w", encoding="ascii") as handle:
    for idx in range(num_reads):
        start = rng.randrange(0, genome_len - read_len + 1)
        read = genome[start : start + read_len]
        handle.write(f"@read_{idx}\n{read}\n+\n{'I' * read_len}\n")
PY

cargo build --manifest-path "$root/Cargo.toml" --bin count_kmers --bin count_gpu --features gpu

echo -e "backend\tstatus\telapsed_seconds\treads\tread_len\tgenome_len\tk\tinput" > "$summary"

run_and_record() {
  local backend="$1"
  shift
  local log="${out_dir}/${backend}.log"
  local status="ok"
  local start end elapsed
  start="$(date +%s.%N)"
  if "$@" > "$log" 2>&1; then
    status="ok"
  else
    status="failed"
  fi
  end="$(date +%s.%N)"
  elapsed="$(awk -v s="$start" -v e="$end" 'BEGIN { printf "%.6f", e - s }')"
  echo -e "${backend}\t${status}\t${elapsed}\t${num_reads}\t${read_len}\t${genome_len}\t${k}\t${reads}" >> "$summary"
  cat "$log"
  test "$status" = "ok"
}

run_and_record cpu "$root/target/debug/count_kmers" "$reads" "$k"
run_and_record opencl "$root/target/debug/count_gpu" "$reads" "$k"

echo "summary=$summary"
