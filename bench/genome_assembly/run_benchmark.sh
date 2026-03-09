#!/bin/bash
# Run genome assembly benchmark comparing Raptor vs SPAdes
# Usage: ./run_benchmark.sh [dataset] [threads]
#   dataset: simulated, drosophila, human_chr21
#   threads: number of threads (default: 8)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
RESULTS_DIR="${SCRIPT_DIR}/results"
RAPTOR_BIN="${SCRIPT_DIR}/../../target/release/raptor"
SUMMARY_SCRIPT="${SCRIPT_DIR}/summarize_results.py"
TIME_BIN=""

DATASET="${1:-simulated}"
THREADS="${2:-8}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Check dependencies
check_dependencies() {
    echo "Checking dependencies..."

    if [ ! -f "${RAPTOR_BIN}" ]; then
        echo "Building Raptor in release mode..."
        (cd "${SCRIPT_DIR}/../.." && cargo build --release)
    fi

    if ! command -v spades.py &> /dev/null; then
        echo "ERROR: SPAdes not found. Please install SPAdes:"
        echo "  conda install -c bioconda spades"
        echo "  or"
        echo "  wget https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Linux.tar.gz"
        exit 1
    fi

    if ! command -v quast.py &> /dev/null; then
        echo "WARNING: QUAST not found. Install for assembly evaluation:"
        echo "  conda install -c bioconda quast"
        QUAST_AVAILABLE=false
    else
        QUAST_AVAILABLE=true
    fi

    if command -v /usr/bin/time &> /dev/null; then
        TIME_BIN="/usr/bin/time"
    elif command -v gtime &> /dev/null; then
        TIME_BIN="$(command -v gtime)"
    else
        TIME_BIN=""
    fi

    echo "Dependencies OK"
}

# Get memory usage (Linux)
get_peak_memory() {
    local pid=$1
    local peak=0
    while kill -0 "$pid" 2>/dev/null; do
        if [ -f "/proc/${pid}/status" ]; then
            local mem=$(grep VmPeak /proc/${pid}/status 2>/dev/null | awk '{print $2}')
            if [ -n "$mem" ] && [ "$mem" -gt "$peak" ]; then
                peak=$mem
            fi
        fi
        sleep 1
    done
    echo $peak
}

# Run SPAdes assembly
run_spades() {
    local reads1=$1
    local reads2=$2
    local output_dir=$3
    local threads=$4

    echo "=== Running SPAdes ==="
    mkdir -p "${output_dir}"

    local start_time=$(date +%s.%N)

    # Run SPAdes
    spades.py \
        -1 "${reads1}" \
        -2 "${reads2}" \
        -o "${output_dir}" \
        -t "${threads}" \
        --careful \
        2>&1 | tee "${output_dir}/spades.log" &

    local pid=$!
    local peak_mem=$(get_peak_memory $pid) &
    local mem_pid=$!

    wait $pid
    local exit_code=$?
    wait $mem_pid 2>/dev/null || true

    local end_time=$(date +%s.%N)
    local elapsed=$(echo "$end_time - $start_time" | bc)

    # Get peak memory from log if available
    if [ -f "${output_dir}/spades.log" ]; then
        peak_mem=$(grep -oP 'Maximum memory used: \K[0-9.]+' "${output_dir}/spades.log" 2>/dev/null || echo "$peak_mem")
    fi

    echo "SPAdes completed in ${elapsed}s"
    echo "${elapsed}" > "${output_dir}/time.txt"
    echo "${peak_mem}" > "${output_dir}/memory.txt"

    return $exit_code
}

# Run Raptor assembly
run_raptor() {
    local reads1=$1
    local reads2=$2
    local output_dir=$3
    local threads=$4

    echo "=== Running Raptor ==="
    mkdir -p "${output_dir}"

    local start_time=$(date +%s.%N)

    local exit_code=0
    local peak_mem="N/A"

    if [ -n "${TIME_BIN}" ]; then
        "${TIME_BIN}" -v "${RAPTOR_BIN}" assemble-large \
            -i "${reads1}" \
            --input2 "${reads2}" \
            -o "${output_dir}/contigs.fa" \
            -t "${threads}" \
            --min-count 2 \
            --scaffold \
            --polish \
            --compress-buckets \
            2>&1 | tee "${output_dir}/raptor.log"
        exit_code=${PIPESTATUS[0]}
        peak_mem=$(grep "Maximum resident set size" "${output_dir}/raptor.log" | awk '{print $6}' | tail -1)
        peak_mem=${peak_mem:-N/A}
    else
        (
            SECONDS=0
            "${RAPTOR_BIN}" assemble-large \
                -i "${reads1}" \
                --input2 "${reads2}" \
                -o "${output_dir}/contigs.fa" \
                -t "${threads}" \
                --min-count 2 \
                --scaffold \
                --polish \
                --compress-buckets \
                2>&1 | tee "${output_dir}/raptor.log"
            echo "ELAPSED_SECONDS=${SECONDS}" >> "${output_dir}/raptor.log"
        )
        exit_code=$?
    fi

    local end_time=$(date +%s.%N)
    local elapsed=$(echo "$end_time - $start_time" | bc)

    echo "Raptor completed in ${elapsed}s"
    echo "${elapsed}" > "${output_dir}/time.txt"
    echo "${peak_mem}" > "${output_dir}/memory.txt"

    return $exit_code
}

parse_quast() {
    local quast_dir=$1
    local output_prefix=$2

    if [ ! -f "${quast_dir}/report.tsv" ]; then
        return 0
    fi

    python3 << EOF
from pathlib import Path

report = Path("${quast_dir}/report.tsv")
rows = []
for line in report.read_text().splitlines():
    if line.strip():
        rows.append(line.split("\t"))

if len(rows) < 2:
    raise SystemExit(0)

header = rows[0]
metrics = rows[1:]
assemblies = header[1:]

for idx, assembly in enumerate(assemblies, start=1):
    name = assembly.strip().lower()
    out = Path(f"${output_prefix}/{name}_quast.tsv")
    with out.open("w") as handle:
        for metric in metrics:
            if idx < len(metric):
                handle.write(f"{metric[0]}\t{metric[idx]}\n")
EOF
}

# Calculate assembly statistics
calc_stats() {
    local assembly=$1
    local output=$2

    python3 << EOF
import sys

def calc_n50(lengths):
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    cumsum = 0
    for l in lengths:
        cumsum += l
        if cumsum >= total / 2:
            return l
    return 0

def calc_stats(fasta_path):
    lengths = []
    current_len = 0

    try:
        with open(fasta_path) as f:
            for line in f:
                if line.startswith('>'):
                    if current_len > 0:
                        lengths.append(current_len)
                    current_len = 0
                else:
                    current_len += len(line.strip())
            if current_len > 0:
                lengths.append(current_len)
    except FileNotFoundError:
        return None

    if not lengths:
        return None

    lengths.sort(reverse=True)

    return {
        'num_contigs': len(lengths),
        'total_length': sum(lengths),
        'largest': lengths[0],
        'n50': calc_n50(lengths),
        'n90': calc_n50(lengths[:int(len(lengths)*0.9)]) if len(lengths) > 10 else lengths[-1],
        'avg_length': sum(lengths) / len(lengths),
        'contigs_1k': sum(1 for l in lengths if l >= 1000),
        'contigs_10k': sum(1 for l in lengths if l >= 10000),
        'contigs_100k': sum(1 for l in lengths if l >= 100000),
    }

stats = calc_stats("${assembly}")
if stats:
    with open("${output}", 'w') as f:
        for k, v in stats.items():
            f.write(f"{k}\t{v}\n")
    print(f"Contigs: {stats['num_contigs']}")
    print(f"Total length: {stats['total_length']:,} bp")
    print(f"Largest: {stats['largest']:,} bp")
    print(f"N50: {stats['n50']:,} bp")
else:
    print("ERROR: Could not read assembly file")
    sys.exit(1)
EOF
}

# Run QUAST comparison
run_quast() {
    local reference=$1
    local spades_assembly=$2
    local raptor_assembly=$3
    local output_dir=$4

    if [ "$QUAST_AVAILABLE" = true ]; then
        echo "=== Running QUAST comparison ==="
        quast.py \
            -r "${reference}" \
            -o "${output_dir}" \
            -t "${THREADS}" \
            --labels "SPAdes,Raptor" \
            "${spades_assembly}" \
            "${raptor_assembly}" \
            2>&1 | tee "${output_dir}/quast.log"
    else
        echo "Skipping QUAST (not installed)"
    fi
}

# Generate comparison report
generate_report() {
    local results_dir=$1
    local dataset=$2

    echo ""
    echo "=============================================="
    echo "       BENCHMARK RESULTS: ${dataset}"
    echo "=============================================="
    echo ""

    # Read stats
    if [ -f "${results_dir}/spades/stats.txt" ]; then
        echo "=== SPAdes ==="
        cat "${results_dir}/spades/stats.txt"
        echo "Time: $(cat ${results_dir}/spades/time.txt 2>/dev/null || echo 'N/A')s"
        echo "Memory: $(cat ${results_dir}/spades/memory.txt 2>/dev/null || echo 'N/A') KB"
        echo ""
    fi

    if [ -f "${results_dir}/raptor/stats.txt" ]; then
        echo "=== Raptor ==="
        cat "${results_dir}/raptor/stats.txt"
        echo "Time: $(cat ${results_dir}/raptor/time.txt 2>/dev/null || echo 'N/A')s"
        echo "Memory: $(cat ${results_dir}/raptor/memory.txt 2>/dev/null || echo 'N/A') KB"
        echo ""
    fi

    # Generate markdown report
    python3 << EOF
import os

results_dir = "${results_dir}"
dataset = "${dataset}"

def read_stats(path):
    stats = {}
    try:
        with open(path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    stats[parts[0]] = parts[1]
    except:
        pass
    return stats

def read_single(path):
    try:
        with open(path) as f:
            return f.read().strip()
    except:
        return 'N/A'

spades = read_stats(f"{results_dir}/spades/stats.txt")
raptor = read_stats(f"{results_dir}/raptor/stats.txt")

spades_time = read_single(f"{results_dir}/spades/time.txt")
raptor_time = read_single(f"{results_dir}/raptor/time.txt")
spades_mem = read_single(f"{results_dir}/spades/memory.txt")
raptor_mem = read_single(f"{results_dir}/raptor/memory.txt")

report = f"""# Genome Assembly Benchmark: {dataset}

## Summary

| Metric | SPAdes | Raptor | Winner |
|--------|--------|--------|--------|
| Contigs | {spades.get('num_contigs', 'N/A')} | {raptor.get('num_contigs', 'N/A')} | - |
| Total Length | {int(spades.get('total_length', 0)):,} bp | {int(raptor.get('total_length', 0)):,} bp | - |
| N50 | {int(spades.get('n50', 0)):,} bp | {int(raptor.get('n50', 0)):,} bp | {'SPAdes' if int(spades.get('n50', 0)) > int(raptor.get('n50', 0)) else 'Raptor'} |
| Largest Contig | {int(spades.get('largest', 0)):,} bp | {int(raptor.get('largest', 0)):,} bp | - |
| Time | {spades_time}s | {raptor_time}s | {'SPAdes' if float(spades_time.replace('N/A','999999')) < float(raptor_time.replace('N/A','999999')) else 'Raptor'} |
| Peak Memory | {spades_mem} KB | {raptor_mem} KB | - |

## Detailed Statistics

### SPAdes
- Contigs >= 1kb: {spades.get('contigs_1k', 'N/A')}
- Contigs >= 10kb: {spades.get('contigs_10k', 'N/A')}
- Contigs >= 100kb: {spades.get('contigs_100k', 'N/A')}

### Raptor
- Contigs >= 1kb: {raptor.get('contigs_1k', 'N/A')}
- Contigs >= 10kb: {raptor.get('contigs_10k', 'N/A')}
- Contigs >= 100kb: {raptor.get('contigs_100k', 'N/A')}
"""

with open(f"{results_dir}/report.md", 'w') as f:
    f.write(report)

print(report)
EOF

    echo ""
    echo "Report saved to: ${results_dir}/report.md"
}

write_run_metadata() {
    local run_dir=$1
    local dataset=$2
    local threads=$3

    cat > "${run_dir}/metadata.txt" << EOF
dataset=${dataset}
threads=${threads}
timestamp=${TIMESTAMP}
raptor_bin=${RAPTOR_BIN}
hostname=$(hostname)
kernel=$(uname -sr)
EOF

    cat > "${run_dir}/commands.sh" << EOF
#!/bin/bash
set -e

spades.py -1 "${DATA_DIR}/${dataset}/reads_1.fastq.gz" -2 "${DATA_DIR}/${dataset}/reads_2.fastq.gz" -o "${run_dir}/spades" -t "${threads}" --careful
"${RAPTOR_BIN}" assemble-large -i "${DATA_DIR}/${dataset}/reads_1.fastq.gz" --input2 "${DATA_DIR}/${dataset}/reads_2.fastq.gz" -o "${run_dir}/raptor/contigs.fa" -t "${threads}" --min-count 2 --scaffold --polish --compress-buckets
EOF
    chmod +x "${run_dir}/commands.sh"
}

# Main benchmark function
run_benchmark() {
    local dataset=$1
    local threads=$2

    local dataset_dir="${DATA_DIR}/${dataset}"
    local run_dir="${RESULTS_DIR}/${dataset}_${TIMESTAMP}"

    echo "=============================================="
    echo "  Genome Assembly Benchmark"
    echo "  Dataset: ${dataset}"
    echo "  Threads: ${threads}"
    echo "  Output:  ${run_dir}"
    echo "=============================================="

    # Check data exists
    if [ ! -d "${dataset_dir}" ]; then
        echo "ERROR: Dataset not found at ${dataset_dir}"
        echo "Run ./download_data.sh ${dataset} first"
        exit 1
    fi

    local reads1="${dataset_dir}/reads_1.fastq.gz"
    local reads2="${dataset_dir}/reads_2.fastq.gz"
    local reference="${dataset_dir}/reference.fa"

    if [ ! -f "${reads1}" ] || [ ! -f "${reads2}" ]; then
        echo "ERROR: Read files not found"
        exit 1
    fi

    mkdir -p "${run_dir}"
    write_run_metadata "${run_dir}" "${dataset}" "${threads}"

    # Run SPAdes
    run_spades "${reads1}" "${reads2}" "${run_dir}/spades" "${threads}"
    calc_stats "${run_dir}/spades/contigs.fasta" "${run_dir}/spades/stats.txt"

    # Run Raptor
    run_raptor "${reads1}" "${reads2}" "${run_dir}/raptor" "${threads}"
    calc_stats "${run_dir}/raptor/contigs.fa" "${run_dir}/raptor/stats.txt"

    # Run QUAST if reference available
    if [ -f "${reference}" ] || [ -f "${reference}.gz" ]; then
        if [ -f "${reference}.gz" ]; then
            gunzip -k "${reference}.gz" 2>/dev/null || true
        fi
        run_quast "${reference}" \
            "${run_dir}/spades/contigs.fasta" \
            "${run_dir}/raptor/contigs.fa" \
            "${run_dir}/quast"
        parse_quast "${run_dir}/quast" "${run_dir}"
    fi

    # Generate report
    generate_report "${run_dir}" "${dataset}"

    if [ -f "${SUMMARY_SCRIPT}" ]; then
        python3 "${SUMMARY_SCRIPT}"
    fi
}

# Main
check_dependencies
run_benchmark "${DATASET}" "${THREADS}"
