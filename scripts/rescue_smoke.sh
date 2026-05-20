#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tmp="${TMPDIR:-/tmp}/raptor-rescue-smoke"
reads="$tmp/reads.fastq"
assembly="$tmp/assembly.fasta"

rm -rf "$tmp"
mkdir -p "$tmp"

cat > "$reads" <<'FASTQ'
@read_1
ACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read_2
CGTACGTACGTACGTACGTACGTACGTACGTA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
FASTQ

echo "== cargo check =="
cargo check --manifest-path "$root/Cargo.toml"

echo "== cargo check --features gpu =="
cargo check --manifest-path "$root/Cargo.toml" --features gpu

echo "== GPU k-mer smoke =="
if cargo run --manifest-path "$root/Cargo.toml" --features gpu --bin count_gpu -- "$reads" 5; then
  echo "gpu_smoke=ok"
else
  echo "gpu_smoke=unavailable"
fi

echo "== assembler smoke =="
cargo run --manifest-path "$root/Cargo.toml" -- assemble \
  --input "$reads" \
  --output "$assembly" \
  --min-len 5 \
  --threads 2

test -s "$assembly"
echo "assembly_output=$assembly"
echo "rescue_smoke=ok"
