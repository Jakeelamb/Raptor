# Raptor HPC Setup Instructions

This guide provides detailed instructions for compiling and running Raptor on HPC systems.

## Compilation Options

Raptor supports various compilation options to adapt to different HPC environments:

| Option | Description |
|--------|-------------|
| `--no-mpi` | Disables MPI support (use when MPI libraries are unavailable) |
| `--gpu` | Enables GPU acceleration using OpenCL |
| `--static` | Creates a statically linked binary for better portability across nodes |
| `--verbose` | Shows more detailed output during compilation |

## Common HPC Environment Setup

Different HPC systems have different module systems. Here are common setups:

### Loading Required Modules

Before compiling or running Raptor, load the appropriate modules:

```bash
# Clean the environment
module purge

# Load essential modules
module load gcc/latest      # or a specific version like gcc/11.2.0

# For MPI support
module load openmpi/latest  # or mpich/latest

# For GPU support
module load cuda/latest     # or a specific version like cuda/12.2
module load opencl          # if available
```

### Checking Available Modules

If you're unsure which modules are available:

```bash
# List all available modules
module avail

# Search for specific modules
module avail mpi
module avail cuda
```

## Handling MPI Issues

If you encounter MPI-related compilation errors:

1. **Check if MPI is available**:
   ```bash
   which mpicc
   mpicc --version
   ```

2. **Try loading different MPI implementations**:
   ```bash
   module load openmpi  # Try OpenMPI
   # or
   module load mpich    # Try MPICH
   ```

3. **Compile without MPI**:
   ```bash
   ./compile_hpc.sh --no-mpi
   ```

4. **Set MPI library path manually**:
   ```bash
   # For OpenMPI
   export MPI_ROOT=/path/to/openmpi
   export PATH=$MPI_ROOT/bin:$PATH
   export LD_LIBRARY_PATH=$MPI_ROOT/lib:$LD_LIBRARY_PATH
   ```

## Handling GPU Issues

If you encounter GPU-related compilation errors:

1. **Check if CUDA is available**:
   ```bash
   which nvcc
   nvcc --version
   ```

2. **Load required modules**:
   ```bash
   module load cuda
   module load opencl  # if available
   ```

3. **Compile without GPU support** if not needed:
   ```bash
   ./compile_hpc.sh --no-mpi  # Just omit the --gpu flag
   ```

## Running on HPC

### SLURM Submission

```bash
# Submit the job
sbatch raptor_hpc.sh /path/to/read1.fastq.gz /path/to/read2.fastq.gz

# Monitor job
./monitor_hpc.sh <job-id>
```

### PBS/Torque Submission

```bash
# Submit the job 
qsub raptor_hpc.sh /path/to/read1.fastq.gz /path/to/read2.fastq.gz

# Monitor job
qstat -f <job-id>
```

## Troubleshooting

### MPI Error: Could not find library

```
error: could not find native static library `mpicxx`, perhaps an -L flag is missing?
error: could not find system library 'mpicxx' required by the 'mpi-sys' crate
```

**Solution**: Try using the `--no-mpi` flag:
```bash
./compile_hpc.sh --no-mpi
```

### Static Linking Errors

If static linking fails, compile without the `--static` flag:
```bash
./compile_hpc.sh --no-mpi --gpu  # Omit --static
```

### Memory Optimization

For extremely large datasets:
1. Enable streaming mode in the job script (enabled by default)
2. Consider distributed mode for assembly with `--distributed --buckets 64` 
3. Monitor memory usage with `./monitor_hpc.sh <job-id>`

## Testing Feature Combinations

You can test different feature combinations with:

```bash
./test_features.sh
```

This will verify that Raptor compiles with various combinations of features. 