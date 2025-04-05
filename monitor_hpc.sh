#!/bin/bash

# monitor_hpc.sh - Tool for monitoring Raptor jobs on HPC systems
# Usage: ./monitor_hpc.sh [job_id]

BLUE='\033[0;34m'
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[0;33m'
NC='\033[0m' # No Color

print_help() {
    echo -e "${BLUE}Usage:${NC} ./monitor_hpc.sh [job_id]"
    echo ""
    echo "Monitor Raptor jobs running on HPC systems."
    echo ""
    echo -e "${BLUE}Options:${NC}"
    echo "  -h, --help    Show this help message and exit"
    echo "  -v, --verbose Show more detailed information"
    echo "  -f, --follow  Continuously monitor job (like 'tail -f')"
    echo ""
    echo -e "${BLUE}Examples:${NC}"
    echo "  ./monitor_hpc.sh 1234          # Monitor SLURM job 1234"
    echo "  ./monitor_hpc.sh 1234.compute  # Monitor PBS job 1234.compute"
    echo "  ./monitor_hpc.sh -f 1234       # Continuously monitor job 1234"
}

detect_scheduler() {
    if command -v squeue &> /dev/null; then
        echo "slurm"
    elif command -v qstat &> /dev/null; then
        if qstat -help 2>&1 | grep -q PBS; then
            echo "pbs"
        else
            echo "sge"
        fi
    else
        echo "unknown"
    fi
}

monitor_slurm_job() {
    local job_id=$1
    local follow=$2
    
    echo -e "${BLUE}SLURM Job Information:${NC}"
    scontrol show job "$job_id"
    
    echo -e "\n${BLUE}Resource Usage:${NC}"
    sstat --format=AveCPU,AveRSS,AveVMSize,MaxRSS,MaxVMSize -j "$job_id"
    
    echo -e "\n${BLUE}Job Output:${NC}"
    if [ "$follow" = true ]; then
        tail -f "slurm-${job_id}.out"
    else
        tail -n 20 "slurm-${job_id}.out"
    fi
}

monitor_pbs_job() {
    local job_id=$1
    local follow=$2
    
    echo -e "${BLUE}PBS Job Information:${NC}"
    qstat -f "$job_id"
    
    # Try to find output file
    local output_file=$(qstat -f "$job_id" | grep -oP "Output_Path = \K.*")
    
    if [ -n "$output_file" ] && [ -f "$output_file" ]; then
        echo -e "\n${BLUE}Job Output:${NC}"
        if [ "$follow" = true ]; then
            tail -f "$output_file"
        else
            tail -n 20 "$output_file"
        fi
    else
        echo -e "\n${YELLOW}Cannot locate output file${NC}"
    fi
}

monitor_sge_job() {
    local job_id=$1
    local follow=$2
    
    echo -e "${BLUE}SGE Job Information:${NC}"
    qstat -j "$job_id"
    
    # Try to find output file based on job ID
    local output_file=$(find . -name "*o${job_id}" -type f | head -n 1)
    
    if [ -n "$output_file" ] && [ -f "$output_file" ]; then
        echo -e "\n${BLUE}Job Output:${NC}"
        if [ "$follow" = true ]; then
            tail -f "$output_file"
        else
            tail -n 20 "$output_file"
        fi
    else
        echo -e "\n${YELLOW}Cannot locate output file${NC}"
    fi
}

# Parse arguments
verbose=false
follow=false
job_id=""

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            print_help
            exit 0
            ;;
        -v|--verbose)
            verbose=true
            shift
            ;;
        -f|--follow)
            follow=true
            shift
            ;;
        *)
            job_id=$1
            shift
            ;;
    esac
done

# Check if job ID was provided
if [ -z "$job_id" ]; then
    echo -e "${RED}Error: No job ID provided${NC}"
    print_help
    exit 1
fi

# Detect which scheduler is being used
scheduler=$(detect_scheduler)

echo -e "${GREEN}Monitoring job ${job_id} on ${scheduler} system${NC}"

# Monitor job based on scheduler type
case $scheduler in
    slurm)
        monitor_slurm_job "$job_id" "$follow"
        ;;
    pbs)
        monitor_pbs_job "$job_id" "$follow"
        ;;
    sge)
        monitor_sge_job "$job_id" "$follow"
        ;;
    *)
        echo -e "${RED}Error: Could not detect HPC scheduler. Please make sure you're on a SLURM, PBS, or SGE system.${NC}"
        exit 1
        ;;
esac

# If not following, check job status
if [ "$follow" = false ]; then
    case $scheduler in
        slurm)
            status=$(squeue -h -j "$job_id" -o "%T" 2>/dev/null)
            if [ -z "$status" ]; then
                echo -e "\n${YELLOW}Job $job_id is not in the queue. It may have completed or failed.${NC}"
                echo "Check sacct -j $job_id for more information."
            else
                echo -e "\n${GREEN}Job $job_id status: $status${NC}"
            fi
            ;;
        pbs|sge)
            status=$(qstat -j "$job_id" 2>&1)
            if [[ "$status" == *"not found"* ]]; then
                echo -e "\n${YELLOW}Job $job_id is not in the queue. It may have completed or failed.${NC}"
            else
                echo -e "\n${GREEN}Job $job_id is still running${NC}"
            fi
            ;;
    esac
fi

exit 0 