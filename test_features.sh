#!/bin/bash

# Test script for verifying optional feature builds

echo "======================================================================="
echo "Testing Raptor builds with different feature combinations"
echo "======================================================================="

# Function to run a test and report result
run_test() {
    local test_name=$1
    local command=$2
    
    echo "-----------------------------------------------------------------------"
    echo "Test: $test_name"
    echo "Command: $command"
    echo "-----------------------------------------------------------------------"
    
    if eval "$command"; then
        echo "RESULT: ✅ PASSED"
    else
        echo "RESULT: ❌ FAILED"
    fi
    echo ""
}

# Clean any previous builds
cargo clean

# Test 1: Default build (with MPI enabled)
run_test "Default Build (MPI enabled)" "cargo build --release --lib"

# Test 2: No MPI build
run_test "No MPI Build" "cargo build --release --no-default-features --lib"

# Test 3: GPU-enabled build
run_test "GPU-enabled Build" "cargo build --release --features gpu --lib"

# Test 4: MPI and GPU enabled
run_test "MPI and GPU enabled" "cargo build --release --features \"mpi-support gpu\" --lib"

echo "======================================================================="
echo "Test Summary"
echo "======================================================================="
echo "All tests completed. Check above for results."
echo "If any test failed, check the error messages for details."
echo ""
echo "If you're on an HPC system, you might need to load appropriate modules:"
echo "  module load openmpi  # For MPI support"
echo "  module load cuda     # For GPU support"
echo "=======================================================================" 