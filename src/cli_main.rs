use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(name = "Raptor", version, about = "High-performance Rust-based assembler", long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Normalize input reads with optional GPU acceleration
    Normalize {
        /// Input read 1 FASTQ(.gz)
        #[arg(short, long)]
        input1: String,

        /// Optional input read 2 for paired-end
        #[arg(short = 'I', long)]
        input2: Option<String>,

        /// Output prefix for normalized files
        #[arg(short, long)]
        output: String,

        /// Enable GPU-based k-mer counting
        #[arg(long)]
        gpu: bool,

        /// Number of threads to use
        #[arg(long, default_value_t = num_cpus::get())]
        threads: usize,
        
        /// Enable streaming mode for large datasets
        #[arg(long)]
        streaming: bool,
    },

    /// Assemble normalized reads into contigs
    Assemble {
        /// Input FASTQ(.gz) file
        #[arg(short, long)]
        input: String,

        /// Output FASTA(.gz) file
        #[arg(short, long)]
        output: String,

        /// Minimum contig length
        #[arg(long, default_value_t = 150)]
        min_len: usize,

        /// Number of threads
        #[arg(long, default_value_t = num_cpus::get())]
        threads: usize,

        /// Output GFA format as well
        #[arg(long)]
        gfa: bool,
        
        /// Enable adaptive k-mer selection
        #[arg(long)]
        adaptive_k: bool,
        
        /// Enable run-length encoding for compression
        #[arg(long)]
        rle: bool,
        
        /// Enable distributed assembly
        #[arg(long)]
        distributed: bool,
        
        /// Number of buckets for distributed assembly
        #[arg(long, default_value_t = 16)]
        buckets: usize,
        
        /// Output GFA2 format as well
        #[arg(long)]
        gfa2: bool,
        
        /// Enable repeat collapsing using RLE
        #[arg(long)]
        collapse_repeats: bool,
        
        /// Minimum repeat length to collapse (in RLE tuples)
        #[arg(long, default_value_t = 20)]
        min_repeat_len: usize,
        
        /// Enable contig polishing
        #[arg(long)]
        polish: bool,
        
        /// Window size for polishing
        #[arg(long, default_value_t = 25)]
        polish_window: usize,
        
        /// Enable streaming mode for large datasets
        #[arg(long)]
        streaming: bool,
        
        /// Export metadata in JSON format
        #[arg(long)]
        export_metadata: bool,

        /// Optional path to write contig metadata as JSON
        #[arg(long)]
        json_metadata: Option<String>,

        /// Optional path to write contig metadata as TSV
        #[arg(long)]
        tsv_metadata: Option<String>,
        
        /// Enable isoform inference and transcript path export
        #[arg(long)]
        isoforms: bool,
        
        /// Optional path to write isoform GTF
        #[arg(long)]
        gtf: Option<String>,
        
        /// Export transcript counts matrix
        #[arg(long)]
        counts_matrix: bool,
        
        /// Optional path to write isoform GFF3
        #[arg(long)]
        gff3: Option<String>,
        
        /// Maximum path depth for isoform traversal
        #[arg(long, default_value_t = 20)]
        max_path_depth: usize,
        
        /// Minimum confidence for keeping isoform paths
        #[arg(long, default_value_t = 0.25)]
        min_confidence: f64,
        
        /// Compute TPM expression values for transcripts
        #[arg(long)]
        compute_tpm: bool,
        
        /// Polish isoform sequences using aligned reads
        #[arg(long)]
        polish_isoforms: bool,
        
        /// CSV file with sample name and SAM alignment path
        #[arg(long, value_name = "CSV")]
        samples: Option<String>,
        
        /// Minimum TPM value for keeping transcripts
        #[arg(long, default_value_t = 0.1)]
        min_tpm: f64,
        
        /// SAM/BAM file with long reads mapped to transcripts for polishing
        #[arg(long)]
        polish_reads: Option<String>,
    },
    
    /// Calculate statistics for assembly output
    Stats {
        /// Input file (FASTA or GFA)
        #[arg(short, long)]
        input: String,

        /// Output format (json or tsv)
        #[arg(long, default_value = "json")]
        format: String,
        
        /// Enable graph stats (branchiness, bubbles)
        #[arg(long)]
        graph: bool,
    },

    /// Benchmark k-mer counting performance
    Benchmark {
        /// Input FASTQ
        #[arg(short, long)]
        input: String,

        /// K-mer size to test
        #[arg(long, default_value_t = 25)]
        k: usize,

        /// Threads
        #[arg(long, default_value_t = num_cpus::get())]
        threads: usize,
    },
    
    /// Reconstruct isoforms from GFA graph and expression data
    Isoform {
        /// Input GFA file containing contigs and overlaps
        #[arg(short, long)]
        input: String,
        
        /// Expression data file (TSV format with contig_id, coverage)
        #[arg(short, long)]
        expression: String,
        
        /// Output prefix for generated files
        #[arg(short, long)]
        output: String,
        
        /// Minimum confidence score for transcript paths (0.0-1.0)
        #[arg(long, default_value_t = 0.25)]
        min_confidence: f64,
        
        /// Maximum depth for graph traversal
        #[arg(long, default_value_t = 50)]
        max_depth: usize,
        
        /// Output format options: fasta,gfa,gtf (comma-separated)
        #[arg(long, default_value = "fasta,gfa")]
        formats: String,
        
        /// Number of threads to use
        #[arg(long, default_value_t = num_cpus::get())]
        threads: usize,
        
        /// Output transcript statistics
        #[arg(long)]
        stats: bool,
        
        /// Enable similarity filtering to remove redundant transcripts
        #[arg(long)]
        filter_similar: bool,
        
        /// Similarity threshold for filtering (0.0-1.0)
        #[arg(long, default_value_t = 0.8)]
        similarity_threshold: f64,
        
        /// Merge similar transcripts instead of filtering them
        #[arg(long)]
        merge_similar: bool,
    },
    
    /// Perform differential expression analysis on transcript counts matrix
    DiffExp {
        /// Input counts matrix file (e.g., output_isoform.counts.matrix)
        #[arg(short, long)]
        matrix: String,
        
        /// Comma-separated list of sample names for group A
        #[arg(long)]
        group_a: String,
        
        /// Comma-separated list of sample names for group B
        #[arg(long)]
        group_b: String,
        
        /// Output file for differential expression results
        #[arg(short, long)]
        output: String,
        
        /// P-value threshold for significance (default: 0.05)
        #[arg(long, default_value_t = 0.05)]
        p_value: f64,
        
        /// Log2 fold-change threshold for significance (default: 1.0)
        #[arg(long, default_value_t = 1.0)]
        fold_change: f64,
    },
    
    /// Compare predicted and truth GTF files to evaluate transcript accuracy
    GtfCompare {
        /// Truth/reference GTF file
        #[arg(short, long)]
        truth: String,
        
        /// Predicted GTF file
        #[arg(short, long)]
        predicted: String,
        
        /// Output file for comparison metrics (optional)
        #[arg(short, long)]
        output: Option<String>,
    },
    
    /// Evaluate assembly results against ground truth
    Eval {
        /// Truth/reference GTF file
        #[arg(long)]
        truth: String,
        
        /// Predicted GTF file
        #[arg(long)]
        pred: String,
        
        /// Output file for evaluation metrics (optional)
        #[arg(long)]
        output: Option<String>,
    },
    
    /// Visualize TPM matrix with PCA plot and heatmap
    Visualize {
        /// Path to TPM matrix
        #[arg(long)]
        matrix: String,
        
        /// Output PCA plot file (SVG or PNG)
        #[arg(long)]
        output: String,
        
        /// Output heatmap file (PNG)
        #[arg(long)]
        heatmap: Option<String>,
        
        /// Output PCA file (PNG)
        #[arg(long)]
        pca: Option<String>,
        
        /// Number of components for PCA (default: 2)
        #[arg(long, default_value_t = 2)]
        components: usize,
    },
    
    /// Traverse paths in a GFA file and export sequences
    Traverse {
        /// Input GFA file with path definitions
        #[arg(short, long)]
        input: String,
        
        /// Segments sequence file (TSV format: segment_id\tsequence)
        #[arg(short, long)]
        segments: String,
        
        /// Output file prefix
        #[arg(short, long)]
        output: String,
        
        /// Export formats (comma-separated, e.g., "fasta,dot,json")
        #[arg(long, default_value = "fasta")]
        formats: String,
        
        /// Include edge information in path
        #[arg(long, default_value_t = false)]
        include_edges: bool,
        
        /// Generate DOT graph visualization
        #[arg(long, default_value_t = false)]
        visualize: bool,
        
        /// Export path metadata
        #[arg(long, default_value_t = false)]
        metadata: bool,
    },
}
