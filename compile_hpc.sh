#!/bin/bash

# Parse command line arguments
GPU_SUPPORT=false
STATIC_BUILD=false
NO_MPI=false
VERBOSE=false

for arg in "$@"
do
    case $arg in
        --gpu)
        GPU_SUPPORT=true
        shift
        ;;
        --static)
        STATIC_BUILD=true
        shift
        ;;
        --no-mpi)
        NO_MPI=true
        shift
        ;;
        --verbose)
        VERBOSE=true
        shift
        ;;
        --help)
        echo "Usage: ./compile_hpc.sh [options]"
        echo "Options:"
        echo "  --gpu      Enable GPU support"
        echo "  --static   Create a static binary for better portability"
        echo "  --no-mpi   Disable MPI support (use if MPI is not available)"
        echo "  --verbose  Show more detailed output"
        echo "  --help     Show this help message"
        exit 0
        ;;
    esac
done

# Module loading examples for different HPC environments
print_module_examples() {
    echo "======================================================================"
    echo "Example module loads for common HPC environments:"
    echo "----------------------------------------------------------------------"
    echo "SLURM / PBS / SGE systems:"
    echo "  module purge"
    echo "  module load gcc/latest     # or gcc/11.2.0"
    echo ""
    echo "For MPI support:"
    echo "  module load openmpi/latest # or mpich/latest"
    echo ""
    echo "For GPU support:"
    echo "  module load cuda/latest    # or cuda/12.2"
    echo "  module load opencl         # if available"
    echo ""
    echo "If you're encountering MPI-related build errors, you can try:"
    echo "1. Check available MPI modules with: module avail mpi"
    echo "2. Load the appropriate MPI module"
    echo "3. Or use the --no-mpi flag to disable MPI support"
    echo "======================================================================"
}

# Check for MPI module
check_mpi() {
    if command -v mpicc >/dev/null 2>&1; then
        echo "MPI compiler found: $(mpicc --version | head -n 1)"
        return 0
    else
        echo "Warning: MPI compiler not found."
        echo "If you encounter MPI-related build errors, try the following:"
        echo "1. Load an MPI module (e.g., 'module load openmpi' or 'module load mpich')"
        echo "2. Or use the --no-mpi flag to disable MPI support"
        if [ "$NO_MPI" = false ]; then
            echo "Proceeding with MPI support enabled..."
        else
            echo "MPI support has been disabled via --no-mpi flag."
        fi
        return 1
    fi
}

# Function to set up build environment
setup_environment() {
    echo "Setting up build environment..."
    
    # Check for required build tools
    if ! command -v cargo >/dev/null 2>&1; then
        echo "Error: Rust/Cargo not found. Please install Rust before proceeding."
        echo "Visit https://rustup.rs/ for installation instructions."
        exit 1
    fi
    
    # Print module examples
    print_module_examples
    
    # Check for GPU toolkit if GPU support is requested
    if [ "$GPU_SUPPORT" = true ]; then
        if ! command -v nvcc >/dev/null 2>&1; then
            echo "Warning: CUDA toolkit not found but GPU support requested."
            echo "If you encounter CUDA-related build errors, load the CUDA module"
            echo "or make sure CUDA toolkit is installed and in your PATH."
        else
            echo "CUDA toolkit found: $(nvcc --version | head -n 1)"
        fi
    fi
    
    # Check for MPI if not disabled
    if [ "$NO_MPI" = false ]; then
        check_mpi
    fi
    
    # Print environment information
    if [ "$VERBOSE" = true ]; then
        echo "Environment details:"
        echo "- Rust version: $(rustc --version)"
        echo "- Cargo version: $(cargo --version)"
        echo "- System: $(uname -a)"
        if [ "$GPU_SUPPORT" = true ] && command -v nvidia-smi >/dev/null 2>&1; then
            echo "- GPU info: $(nvidia-smi --query-gpu=name --format=csv,noheader)"
        fi
        echo "- Build directory: $(pwd)"
    fi
}

# Main build function
build_raptor() {
    echo "Building Raptor for HPC environment..."
    
    # Set up environment
    setup_environment
    
    # Base build command
    BUILD_CMD="cargo build --release"
    
    # Add features for GPU if requested
    if [ "$GPU_SUPPORT" = true ]; then
        BUILD_CMD="$BUILD_CMD --features gpu"
    fi
    
    # Disable MPI if requested
    if [ "$NO_MPI" = false ]; then
        BUILD_CMD="$BUILD_CMD --features mpi-support"
    else
        BUILD_CMD="$BUILD_CMD --no-default-features"
    fi
    
    # Add rustflags for static build if requested
    if [ "$STATIC_BUILD" = true ]; then
        # Note: This works for Linux systems
        export RUSTFLAGS="-C target-feature=+crt-static"
        echo "Static build enabled. Adding RUSTFLAGS: $RUSTFLAGS"
    fi
    
    # Set environment variables for build optimization
    export RUSTFLAGS="${RUSTFLAGS} -C opt-level=3 -C target-cpu=native"
    export CARGO_PROFILE_RELEASE_LTO=true
    export CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1
    
    # Print the final build command
    echo "Executing: $BUILD_CMD"
    
    # Execute the build
    if [ "$VERBOSE" = true ]; then
        $BUILD_CMD -v
    else
        $BUILD_CMD
    fi
    
    # Check if build was successful
    if [ $? -eq 0 ]; then
        echo "Build successful!"
        echo "Binary location: $(pwd)/target/release/raptor"
        
        # Create a directory to copy files for HPC deployment
        mkdir -p raptor_hpc_deploy
        cp target/release/raptor raptor_hpc_deploy/
        cp raptor_hpc.sh raptor_hpc_deploy/
        cp monitor_hpc.sh raptor_hpc_deploy/
        chmod +x raptor_hpc_deploy/monitor_hpc.sh
        
        echo "HPC deployment files prepared in: raptor_hpc_deploy/"
        echo "Transfer this directory to your HPC system and run with:"
        echo "  sbatch raptor_hpc.sh <read1.fastq.gz> <read2.fastq.gz>"
    else
        echo "Build failed!"
        return 1
    fi
}

# Execute the build
build_raptor 