# Development Notes: HPC Support

## Overview of Changes

### 1. Optional Dependencies
- Modified `Cargo.toml` to make MPI and GPU support optional features
- Added feature flags to control which dependencies are included in the build
- Default build includes MPI support, which can be disabled with `--no-default-features`

### 2. HPC Build Scripts
- Created `compile_hpc.sh` with environment module detection for different HPC systems
- Added command-line flags (`--no-mpi`, `--gpu`, `--static`, `--verbose`) for flexible compilation
- Implemented module loading suggestions for various HPC environments

### 3. Job Monitoring
- Developed `monitor_hpc.sh` to track and display job statistics across different schedulers (SLURM, PBS, SGE)
- Added support for continuous monitoring mode with `-f/--follow` flag
- Provided output file inspection for easier debugging

### 4. Feature Testing
- Created `test_features.sh` to validate builds with different feature combinations
- Tests all combinations: default, no-MPI, GPU-only, and MPI+GPU
- Ensures that all feature combinations compile correctly

### 5. Documentation
- Added `HPC_INSTRUCTIONS.md` with detailed setup instructions for HPC environments
- Updated main README.md to include information about HPC support
- Added troubleshooting guidance for common MPI and GPU issues

## Current Limitations

1. **Executable Errors**: The main executable still has compilation errors that need to be fixed
2. **CI/CD Integration**: No automated testing for MPI functionality in CI pipeline
3. **Performance Benchmarks**: Need to establish baseline performance metrics for MPI vs non-MPI modes

## Future Development Tasks

### Immediate Priorities
- [ ] Fix compilation errors in the main executable
- [ ] Add error handling for MPI initialization failures
- [ ] Create more sophisticated MPI partitioning strategies for large graphs

### Medium-term Tasks
- [ ] Add automated tests for MPI support in CI/CD pipeline
- [ ] Implement dynamic load balancing for MPI processes
- [ ] Optimize GPU kernel for different CUDA architectures
- [ ] Add support for OpenCL acceleration as an alternative to CUDA

### Long-term Goals
- [ ] Develop hybrid MPI+GPU computation model for maximum performance
- [ ] Implement checkpoint/restart functionality for long-running jobs
- [ ] Create benchmark suite specific to HPC environments
- [ ] Support for container-based deployment (Singularity/Apptainer)

## Contributing to HPC Support

When working on HPC support features:

1. Always test with and without MPI/GPU features enabled
2. Use the `test_features.sh` script to validate your changes
3. Update documentation when adding new flags or features
4. Consider compatibility with different module systems and schedulers
5. Follow existing error handling patterns for consistency

## Architecture Considerations

### Conditional Compilation
The codebase uses Rust's feature flags for conditional compilation:
```rust
#[cfg(feature = "mpi-support")]
use mpi::traits::*;

#[cfg(feature = "mpi-support")]
fn distribute_work() {
    // MPI implementation
}

#[cfg(not(feature = "mpi-support"))]
fn distribute_work() {
    // Non-MPI fallback implementation
}
```

### Function Arguments
When MPI functionality requires additional arguments, consider using Option types:
```rust
fn process_data(
    data: &[u8],
    #[cfg(feature = "mpi-support")] mpi_comm: Option<&mpi::Communicator>,
) {
    // Function implementation
}
```

This approach maintains compatibility while allowing for MPI-specific enhancements. 