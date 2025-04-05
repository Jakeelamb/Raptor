#!/bin/bash

# Simple memory monitoring script for Raptor pipeline
# Usage: ./monitor_memory.sh <process_name> <output_file>

if [ $# -lt 2 ]; then
    echo "Usage: $0 <process_name> <output_file>"
    echo "Example: $0 normalize_paired_reads memory_usage.csv"
    exit 1
fi

PROCESS_NAME=$1
OUTPUT_FILE=$2

# Create CSV header
echo "timestamp,process,rss_mb,vsz_mb,cpu_percent" > "$OUTPUT_FILE"

echo "Starting memory monitoring for processes containing '$PROCESS_NAME'"
echo "Press Ctrl+C to stop monitoring"

# Monitor memory usage every 2 seconds
while true; do
    TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
    # Use ps to find process and get memory info
    ps aux | grep "$PROCESS_NAME" | grep -v "grep" | grep -v "monitor_memory" | 
    while read -r line; do
        # Extract process info
        PROC_USER=$(echo "$line" | awk '{print $1}')
        PROC_PID=$(echo "$line" | awk '{print $2}')
        PROC_CPU=$(echo "$line" | awk '{print $3}')
        PROC_MEM=$(echo "$line" | awk '{print $4}')
        PROC_VSZ=$(echo "$line" | awk '{print $5}')
        PROC_RSS=$(echo "$line" | awk '{print $6}')
        PROC_CMD=$(echo "$line" | awk '{$1=$2=$3=$4=$5=$6=$7=$8=$9=$10=""; print $0}' | sed 's/^ *//')
        
        # Convert KB to MB
        RSS_MB=$(echo "scale=2; $PROC_RSS / 1024" | bc)
        VSZ_MB=$(echo "scale=2; $PROC_VSZ / 1024" | bc)
        
        # Append to CSV
        echo "$TIMESTAMP,$PROC_CMD,$RSS_MB,$VSZ_MB,$PROC_CPU" >> "$OUTPUT_FILE"
    done
    
    sleep 2
done 