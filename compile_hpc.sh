#!/bin/bash

# Script to compile Raptor with optimizations for HPC systems
# Usage: ./compile_hpc.sh [--static] [--gpu]

# Parse options
STATIC=false
GPU=false

for arg in "$@"; do
  case $arg in
    --static)
      STATIC=true
      shift
      ;;
    --gpu)
      GPU=true
      shift
      ;;
  esac
done

# Set optimization flags
export RUSTFLAGS="-C target-cpu=native -C opt-level=3 -C codegen-units=1 -C lto=fat"

echo "=========================================================="
echo "Compiling Raptor with HPC optimizations"
echo "=========================================================="
echo "Optimization flags: $RUSTFLAGS"
echo "Static linking: $STATIC"
echo "GPU support: $GPU"
echo "=========================================================="

# Create build command based on options
BUILD_CMD="cargo build --release"

# Add GPU feature if requested
if [ "$GPU" = true ]; then
  BUILD_CMD="$BUILD_CMD --features=gpu"
  echo "Including GPU acceleration support"
fi

# Add static target if requested
if [ "$STATIC" = true ]; then
  BUILD_CMD="$BUILD_CMD --target x86_64-unknown-linux-musl"
  echo "Creating statically linked binary (better portability across nodes)"
fi

echo "Running: $BUILD_CMD"
echo "=========================================================="

# Execute build
eval $BUILD_CMD

# Check build status
if [ $? -eq 0 ]; then
  echo "=========================================================="
  echo "Compilation successful!"
  
  if [ "$STATIC" = true ]; then
    echo "Binary location: target/x86_64-unknown-linux-musl/release/raptor"
    echo "Copy this binary to your HPC system for better portability"
  else
    echo "Binary location: target/release/raptor"
  fi
  
  echo "=========================================================="
else
  echo "=========================================================="
  echo "Compilation failed!"
  echo "=========================================================="
  exit 1
fi 