#!/bin/bash

# Test script for running the low-memory pipeline with real data
# and monitoring memory usage

set -e

# Start the memory monitor in the background
./monitor_memory.sh normalize_paired_reads memory_normalize.csv &
MONITOR_PID=$!

# Ensure we kill the memory monitor on exit
trap "kill $MONITOR_PID 2>/dev/null || true" EXIT

# Edit the low_memory_run.sh script to use real data
sed -i 's/TEST_MODE=true/TEST_MODE=false/' low_memory_run.sh

echo "Starting low-memory pipeline with real data..."
echo "Memory monitoring active, results will be saved to memory_normalize.csv"

# Run the low-memory pipeline
./low_memory_run.sh

# Wait a moment to ensure final memory readings are captured
sleep 5

# Restore the test mode setting
sed -i 's/TEST_MODE=false/TEST_MODE=true/' low_memory_run.sh

echo "Test complete!"
echo "Memory usage data saved to memory_normalize.csv"
echo "You can analyze this file to see the memory usage pattern." 