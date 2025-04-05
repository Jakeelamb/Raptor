use std::process::{Command, Child};
use std::path::Path;
use rayon::prelude::*;
use std::collections::HashMap;
use std::time::Duration;
use tracing::info;

/// Configuration for distributed job scheduling
pub struct DistributedConfig {
    /// Number of buckets for partitioning
    pub buckets: usize,
    /// Number of threads per job
    pub threads_per_job: usize,
    /// Use slurm (true) or local parallel execution (false)
    pub use_slurm: bool,
    /// Maximum runtime per job in minutes
    pub max_runtime: u32,
    /// Memory per job in GB
    pub memory_per_job: u32,
}

impl Default for DistributedConfig {
    fn default() -> Self {
        Self {
            buckets: 16,
            threads_per_job: 1,
            use_slurm: false,
            max_runtime: 60,
            memory_per_job: 4,
        }
    }
}

/// Run assembly jobs on partitioned data in parallel using rayon
pub fn run_parallel_assembly(
    partition_files: &[String],
    output_dir: &str,
    min_len: usize,
    threads_per_job: usize,
) -> Vec<String> {
    info!("Running parallel assembly on {} partitions", partition_files.len());
    
    let outputs: Vec<String> = partition_files
        .par_iter()
        .map(|input_file| {
            let file_stem = Path::new(input_file)
                .file_stem()
                .unwrap_or_default()
                .to_string_lossy();
                
            let output_file = format!("{}/{}.fasta", output_dir, file_stem);
            let output_gfa = format!("{}/{}.gfa", output_dir, file_stem);
            
            // Run the assembly for this partition
            let status = Command::new("cargo")
                .args(&[
                    "run", "--release", "--",
                    "assemble",
                    "-i", input_file,
                    "-o", &output_file,
                    "--min-len", &min_len.to_string(),
                    "--threads", &threads_per_job.to_string(),
                    "--gfa",
                ])
                .status()
                .expect("Failed to execute assembly command");
                
            if !status.success() {
                eprintln!("Warning: Assembly of {} failed with status: {}", input_file, status);
            }
            
            output_file
        })
        .collect();
        
    outputs
}

/// Generate a Slurm job script
fn generate_slurm_script(
    job_name: &str,
    input_file: &str,
    output_file: &str,
    min_len: usize,
    threads: usize,
    runtime_mins: u32,
    memory_gb: u32,
) -> String {
    format!(
        "#!/bin/bash\n\
        #SBATCH --job-name={}\n\
        #SBATCH --output=%x.%j.out\n\
        #SBATCH --error=%x.%j.err\n\
        #SBATCH --time={}:00\n\
        #SBATCH --cpus-per-task={}\n\
        #SBATCH --mem={}G\n\
        \n\
        cargo run --release -- assemble \\\n\
            -i {} \\\n\
            -o {} \\\n\
            --min-len {} \\\n\
            --threads {} \\\n\
            --gfa\n",
        job_name, runtime_mins, threads, memory_gb, input_file, output_file, min_len, threads
    )
}

/// Submit jobs to Slurm for distributed execution
pub fn submit_slurm_jobs(
    partition_files: &[String],
    output_dir: &str,
    min_len: usize,
    config: &DistributedConfig,
) -> Vec<String> {
    info!("Submitting {} Slurm jobs", partition_files.len());
    
    let mut job_ids = Vec::new();
    let mut output_files = Vec::new();
    
    for (i, input_file) in partition_files.iter().enumerate() {
        let job_name = format!("assemble_{}", i);
        let file_stem = Path::new(input_file)
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy();
            
        let output_file = format!("{}/{}.fasta", output_dir, file_stem);
        output_files.push(output_file.clone());
        
        let script = generate_slurm_script(
            &job_name,
            input_file,
            &output_file,
            min_len,
            config.threads_per_job,
            config.max_runtime,
            config.memory_per_job,
        );
        
        // Write the job script to a file
        let script_file = format!("{}/job_{}.sh", output_dir, i);
        std::fs::write(&script_file, script).expect("Failed to write job script");
        
        // Submit the job
        let output = Command::new("sbatch")
            .arg(&script_file)
            .output()
            .expect("Failed to submit Slurm job");
            
        if output.status.success() {
            let stdout = String::from_utf8_lossy(&output.stdout);
            if let Some(job_id) = stdout.trim().split_whitespace().last() {
                job_ids.push(job_id.to_string());
            }
        } else {
            eprintln!("Failed to submit job for {}: {:?}", input_file, output);
        }
    }
    
    info!("Submitted {} Slurm jobs with IDs: {:?}", job_ids.len(), job_ids);
    output_files
}

/// Run the distributed assembly process
pub fn run_distributed_assembly(
    input_sequences: &[String],
    output_dir: &str,
    min_len: usize,
    config: &DistributedConfig,
) -> Vec<String> {
    // Create the output directory
    std::fs::create_dir_all(output_dir).expect("Failed to create output directory");
    
    // Partition the sequences
    let partitions = crate::dist::partition::partition_by_minimizer(
        input_sequences, 
        21, // k-mer size for partitioning
        config.buckets
    );
    
    // Save partitions to files
    let partition_dir = format!("{}/partitions", output_dir);
    let partition_files = crate::dist::partition::save_partitions(
        &partitions,
        &partition_dir,
        "part"
    );
    
    // Run the assemblies
    let output_files = if config.use_slurm {
        submit_slurm_jobs(&partition_files, output_dir, min_len, config)
    } else {
        run_parallel_assembly(&partition_files, output_dir, min_len, config.threads_per_job)
    };
    
    output_files
}
