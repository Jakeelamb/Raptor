#!/bin/bash
# Usage: ./monitor_hpc.sh <job-id>

if [ $# -ne 1 ]; then
    echo "Usage: ./monitor_hpc.sh <job-id>"
    exit 1
fi

JOB_ID=$1
INTERVAL=10  # seconds between checks
LOG_FILE="raptor_monitor_${JOB_ID}.log"

echo "Starting monitoring for Raptor job $JOB_ID"
echo "Results will be logged to $LOG_FILE"
echo "Press Ctrl+C to stop monitoring"

# Header for log file
echo "Timestamp,CPU%,Memory(GB),GPU%,GPUMem(GB),DiskIO(MB/s)" > $LOG_FILE

# Monitor continuously
while true; do
    # Check if job is still running
    squeue -j $JOB_ID &>/dev/null
    if [ $? -ne 0 ]; then
        echo "Job $JOB_ID is no longer running. Monitoring stopped."
        break
    fi
    
    # Get timestamp
    TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
    
    # Get CPU and memory usage using sstat
    STATS=$(sstat --format=AveCPU,AveRSS -j $JOB_ID.batch -n)
    CPU=$(echo $STATS | awk '{print $1}')
    MEM=$(echo $STATS | awk '{print $2}' | sed 's/K//' | awk '{printf "%.2f", $1/1024/1024}')
    
    # Get GPU usage if available
    if command -v nvidia-smi &>/dev/null; then
        NODE=$(squeue -j $JOB_ID -h -o %N)
        if [ ! -z "$NODE" ]; then
            GPU_STATS=$(ssh $NODE "nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader,nounits")
            GPU_UTIL=$(echo $GPU_STATS | awk -F, '{print $1}')
            GPU_MEM=$(echo $GPU_STATS | awk -F, '{print $2}' | awk '{printf "%.2f", $1/1024}')
        else
            GPU_UTIL="N/A"
            GPU_MEM="N/A"
        fi
    else
        GPU_UTIL="N/A"
        GPU_MEM="N/A"
    fi
    
    # Get disk I/O
    DISK_IO=$(ssh $NODE "iostat -m | grep -A 1 avg-cpu | tail -1" | awk '{print $6}')
    
    # Log to file
    echo "$TIMESTAMP,$CPU,$MEM,$GPU_UTIL,$GPU_MEM,$DISK_IO" >> $LOG_FILE
    
    # Display current stats
    echo -e "Time: $TIMESTAMP | CPU: $CPU% | Memory: ${MEM}GB | GPU: $GPU_UTIL% | GPU Mem: ${GPU_MEM}GB | Disk: ${DISK_IO}MB/s"
    
    sleep $INTERVAL
done

echo "Monitoring complete. Results saved to $LOG_FILE" 