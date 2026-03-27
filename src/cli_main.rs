use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(name = "Raptor", version, about = "High-performance Rust-based assembler", long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum ComponentBenchCommands {
    /// Evaluate the shared read mapper against a task directory or task root
    ReadMapping {
        /// Task directory, or a parent directory containing multiple task directories
        #[arg(long)]
        task: String,

        /// Minimizer k-mer size
        #[arg(long)]
        k: Option<usize>,

        /// Minimizer window size
        #[arg(long)]
        w: Option<usize>,

        /// Minimum supporting minimizer hits for a primary mapping
        #[arg(long, default_value_t = 3)]
        min_primary_matches: usize,

        /// Minimum supporting minimizer hits for a scaffold-level contig hit
        #[arg(long, default_value_t = 2)]
        min_scaffold_matches: usize,

        /// Positional tolerance in base pairs for the "near" score
        #[arg(long, default_value_t = 8)]
        position_tolerance: usize,

        /// Emit JSON to stdout instead of the human-readable summary
        #[arg(long)]
        json: bool,

        /// Optional path to write the JSON report
        #[arg(long)]
        output: Option<String>,
    },

    /// Evaluate branch-choice heuristics against a task directory or task root
    BranchResolution {
        /// Task directory, or a parent directory containing multiple task directories
        #[arg(long)]
        task: String,

        /// Minimum read support required before branch threading overrides coverage
        #[arg(long, default_value_t = 1)]
        branch_support_min_win: u32,

        /// Minimum lead over the runner-up support before branch threading overrides coverage
        #[arg(long, default_value_t = 2)]
        branch_support_min_margin: u32,

        /// Disable the non-repeat preference in coverage-driven branch choice
        #[arg(long)]
        disable_prefer_non_repeat: bool,

        /// Emit JSON to stdout instead of the human-readable summary
        #[arg(long)]
        json: bool,

        /// Optional path to write the JSON report
        #[arg(long)]
        output: Option<String>,
    },

    /// Evaluate the shared scaffold+polish postprocess against a task directory or task root
    ScaffoldPolish {
        /// Task directory, or a parent directory containing multiple task directories
        #[arg(long)]
        task: String,

        /// Minimum number of read pairs required to keep a scaffold link
        #[arg(long, default_value_t = 3)]
        min_scaffold_links: usize,

        /// Minimum supporting minimizer hits for a primary mapping
        #[arg(long, default_value_t = 2)]
        min_primary_matches: usize,

        /// Minimum supporting minimizer hits for a scaffold-level contig hit
        #[arg(long, default_value_t = 4)]
        min_scaffold_matches: usize,

        /// Emit JSON to stdout instead of the human-readable summary
        #[arg(long)]
        json: bool,

        /// Optional path to write the JSON report
        #[arg(long)]
        output: Option<String>,
    },

    /// Evaluate phase-3 error correction against a task directory or task root
    ErrorCorrection {
        /// Task directory, or a parent directory containing multiple task directories
        #[arg(long)]
        task: String,

        /// Minimum count threshold used to derive trusted k-mers
        #[arg(long, default_value_t = 1)]
        min_count: u32,

        /// Minimum count floor required for a k-mer to be considered trusted
        #[arg(long, default_value_t = 4)]
        min_trusted_count: u32,

        /// Emit JSON to stdout instead of the human-readable summary
        #[arg(long)]
        json: bool,

        /// Optional path to write the JSON report
        #[arg(long)]
        output: Option<String>,
    },

    /// Evaluate contig extraction against a task directory or task root
    ContigExtraction {
        /// Task directory, or a parent directory containing multiple task directories
        #[arg(long)]
        task: String,

        /// Disable prioritizing higher-coverage seeds first
        #[arg(long)]
        disable_prefer_high_count_seeds: bool,

        /// Disable prioritizing non-repeat seeds before repeat seeds
        #[arg(long)]
        disable_prefer_non_repeat_seeds: bool,

        /// Disable the repeat-seed completion pass
        #[arg(long)]
        disable_repeat_seed_completion: bool,

        /// Suppress redundant contained or same-flank branch-alternative contigs
        #[arg(long)]
        suppress_redundant_contigs: bool,

        /// Emit JSON to stdout instead of the human-readable summary
        #[arg(long)]
        json: bool,

        /// Optional path to write the JSON report
        #[arg(long)]
        output: Option<String>,
    },

    /// Generate a small synthetic read-mapping task directory
    PrepareReadMapping {
        /// Output directory for the generated task
        #[arg(long)]
        output: String,

        /// Generate paired-end reads instead of single-end reads
        #[arg(long)]
        paired: bool,

        /// Number of synthetic contigs
        #[arg(long, default_value_t = 4)]
        num_contigs: usize,

        /// Length of each synthetic contig
        #[arg(long, default_value_t = 3000)]
        contig_len: usize,

        /// Read length
        #[arg(long, default_value_t = 150)]
        read_len: usize,

        /// Number of reads to generate, or number of pairs when --paired is set
        #[arg(long, default_value_t = 256)]
        reads: usize,

        /// Insert size for paired-end generation
        #[arg(long, default_value_t = 450)]
        insert_size: usize,

        /// Length of the shared repeat planted into each contig
        #[arg(long, default_value_t = 80)]
        repeat_len: usize,

        /// Per-base substitution error rate
        #[arg(long, default_value_t = 0.01)]
        error_rate: f64,

        /// Ratio of extra random decoy templates to add as unmapped reads
        #[arg(long, default_value_t = 0.0)]
        decoy_rate: f64,

        /// Ratio of repeat-only ambiguous decoys that should stay unmapped
        #[arg(long, default_value_t = 0.0)]
        ambiguous_repeat_decoy_rate: f64,

        /// Deterministic RNG seed
        #[arg(long, default_value_t = 7)]
        seed: u64,

        /// Minimizer k-mer size to record in task metadata
        #[arg(long, default_value_t = 15)]
        k: usize,

        /// Minimizer window size to record in task metadata
        #[arg(long, default_value_t = 10)]
        w: usize,
    },

    /// Generate a small synthetic branch-resolution task directory
    PrepareBranchResolution {
        /// Output directory for the generated task
        #[arg(long)]
        output: String,

        /// Number of branch-choice cases to generate
        #[arg(long, default_value_t = 64)]
        cases: usize,

        /// Deterministic RNG seed
        #[arg(long, default_value_t = 7)]
        seed: u64,

        /// Scenario profile to generate: balanced or support_stress
        #[arg(long, default_value = "balanced", value_parser = ["balanced", "support_stress"])]
        scenario_profile: String,
    },

    /// Generate a synthetic shared scaffold+polish task directory
    PrepareScaffoldPolish {
        /// Output directory for the generated task
        #[arg(long)]
        output: String,

        /// Number of independent scaffold groups in the synthetic truth
        #[arg(long, default_value_t = 2)]
        scaffold_groups: usize,

        /// Number of contigs per scaffold group
        #[arg(long, default_value_t = 3)]
        contigs_per_scaffold: usize,

        /// Length of each synthetic contig
        #[arg(long, default_value_t = 400)]
        contig_len: usize,

        /// Read length for generated pairs
        #[arg(long, default_value_t = 80)]
        read_len: usize,

        /// Insert size for same-contig estimation pairs
        #[arg(long, default_value_t = 200)]
        insert_size: usize,

        /// Number of same-contig pairs per contig for insert estimation and polishing
        #[arg(long, default_value_t = 16)]
        internal_pairs_per_contig: usize,

        /// Number of true cross-contig pairs generated per real adjacency
        #[arg(long, default_value_t = 5)]
        true_link_pairs: usize,

        /// Number of false cross-contig decoy pairs generated across scaffold groups
        #[arg(long, default_value_t = 2)]
        decoy_link_pairs: usize,

        /// Number of substitutions introduced into each contig before polishing
        #[arg(long, default_value_t = 1)]
        mutations_per_contig: usize,

        /// Deterministic RNG seed
        #[arg(long, default_value_t = 7)]
        seed: u64,
    },

    /// Generate a synthetic phase-3 error-correction task directory
    PrepareErrorCorrection {
        /// Output directory for the generated task
        #[arg(long)]
        output: String,

        /// K-mer size for the synthetic count table
        #[arg(long, default_value_t = 21)]
        k: usize,

        /// Number of trusted root k-mers with attached singleton errors
        #[arg(long, default_value_t = 64)]
        trusted_roots: usize,

        /// Number of weak root k-mers with protected singleton neighbors
        #[arg(long, default_value_t = 32)]
        weak_roots: usize,

        /// Number of correctable singleton errors generated per trusted root
        #[arg(long, default_value_t = 2)]
        correctable_singletons_per_trusted: usize,

        /// Number of protected singleton neighbors generated per weak root
        #[arg(long, default_value_t = 1)]
        protected_singletons_per_weak: usize,

        /// Deterministic RNG seed
        #[arg(long, default_value_t = 7)]
        seed: u64,
    },

    /// Generate a synthetic contig-extraction task directory
    PrepareContigExtraction {
        /// Output directory for the generated task
        #[arg(long)]
        output: String,

        /// Scenario profile to generate
        #[arg(
            long,
            default_value = "branching",
            value_parser = [
                "branching",
                "repeat_fallback",
                "repeat_completion",
                "repeat_priority",
            ]
        )]
        profile: String,

        /// K-mer size for the synthetic graph
        #[arg(long, default_value_t = 11)]
        k: usize,

        /// Number of branching components for branching-profile tasks
        #[arg(long, default_value_t = 6)]
        component_count: usize,

        /// Number of primary-path reads per branching component
        #[arg(long, default_value_t = 6)]
        primary_reads_per_component: usize,

        /// Number of alternate-path reads per branching component
        #[arg(long, default_value_t = 2)]
        alternate_reads_per_component: usize,

        /// Deterministic RNG seed
        #[arg(long, default_value_t = 7)]
        seed: u64,
    },

    /// Generate a panel of synthetic read-mapping task directories
    PrepareReadMappingPanel {
        /// Output directory for the task panel
        #[arg(long)]
        output: String,

        /// Number of task directories to create
        #[arg(long, default_value_t = 16)]
        tasks: usize,

        /// Generate paired-end reads instead of single-end reads
        #[arg(long)]
        paired: bool,

        /// Number of synthetic contigs per task
        #[arg(long, default_value_t = 4)]
        num_contigs: usize,

        /// Length of each synthetic contig
        #[arg(long, default_value_t = 3000)]
        contig_len: usize,

        /// Read length
        #[arg(long, default_value_t = 150)]
        read_len: usize,

        /// Number of reads per task, or number of pairs when --paired is set
        #[arg(long, default_value_t = 256)]
        reads: usize,

        /// Insert size for paired-end generation
        #[arg(long, default_value_t = 450)]
        insert_size: usize,

        /// Length of the shared repeat planted into each contig
        #[arg(long, default_value_t = 80)]
        repeat_len: usize,

        /// Per-base substitution error rate
        #[arg(long, default_value_t = 0.01)]
        error_rate: f64,

        /// Ratio of extra random decoy templates to add as unmapped reads
        #[arg(long, default_value_t = 0.0)]
        decoy_rate: f64,

        /// Ratio of repeat-only ambiguous decoys that should stay unmapped
        #[arg(long, default_value_t = 0.0)]
        ambiguous_repeat_decoy_rate: f64,

        /// Starting RNG seed for task_000
        #[arg(long, default_value_t = 7)]
        seed: u64,

        /// Increment applied to the seed for each subsequent task
        #[arg(long, default_value_t = 1)]
        seed_step: u64,

        /// Minimizer k-mer size to record in task metadata
        #[arg(long, default_value_t = 15)]
        k: usize,

        /// Minimizer window size to record in task metadata
        #[arg(long, default_value_t = 10)]
        w: usize,
    },

    /// Generate a panel of synthetic branch-resolution task directories
    PrepareBranchResolutionPanel {
        /// Output directory for the task panel
        #[arg(long)]
        output: String,

        /// Number of task directories to create
        #[arg(long, default_value_t = 16)]
        tasks: usize,

        /// Number of branch-choice cases per task
        #[arg(long, default_value_t = 64)]
        cases_per_task: usize,

        /// Starting RNG seed for task_000
        #[arg(long, default_value_t = 7)]
        seed: u64,

        /// Increment applied to the seed for each subsequent task
        #[arg(long, default_value_t = 1)]
        seed_step: u64,

        /// Scenario profile to generate: balanced or support_stress
        #[arg(long, default_value = "balanced", value_parser = ["balanced", "support_stress"])]
        scenario_profile: String,
    },

    /// Generate a panel of synthetic shared scaffold+polish task directories
    PrepareScaffoldPolishPanel {
        /// Output directory for the task panel
        #[arg(long)]
        output: String,

        /// Number of task directories to create
        #[arg(long, default_value_t = 16)]
        tasks: usize,

        /// Number of independent scaffold groups in each task
        #[arg(long, default_value_t = 2)]
        scaffold_groups: usize,

        /// Number of contigs per scaffold group
        #[arg(long, default_value_t = 3)]
        contigs_per_scaffold: usize,

        /// Length of each synthetic contig
        #[arg(long, default_value_t = 400)]
        contig_len: usize,

        /// Read length for generated pairs
        #[arg(long, default_value_t = 80)]
        read_len: usize,

        /// Insert size for same-contig estimation pairs
        #[arg(long, default_value_t = 200)]
        insert_size: usize,

        /// Number of same-contig pairs per contig for insert estimation and polishing
        #[arg(long, default_value_t = 16)]
        internal_pairs_per_contig: usize,

        /// Number of true cross-contig pairs generated per real adjacency
        #[arg(long, default_value_t = 5)]
        true_link_pairs: usize,

        /// Number of false cross-contig decoy pairs generated across scaffold groups
        #[arg(long, default_value_t = 2)]
        decoy_link_pairs: usize,

        /// Number of substitutions introduced into each contig before polishing
        #[arg(long, default_value_t = 1)]
        mutations_per_contig: usize,

        /// Starting RNG seed for task_000
        #[arg(long, default_value_t = 7)]
        seed: u64,

        /// Increment applied to the seed for each subsequent task
        #[arg(long, default_value_t = 1)]
        seed_step: u64,
    },

    /// Generate a panel of synthetic phase-3 error-correction task directories
    PrepareErrorCorrectionPanel {
        /// Output directory for the task panel
        #[arg(long)]
        output: String,

        /// Number of task directories to create
        #[arg(long, default_value_t = 16)]
        tasks: usize,

        /// K-mer size for the synthetic count table
        #[arg(long, default_value_t = 21)]
        k: usize,

        /// Number of trusted root k-mers per task
        #[arg(long, default_value_t = 64)]
        trusted_roots: usize,

        /// Number of weak root k-mers per task
        #[arg(long, default_value_t = 32)]
        weak_roots: usize,

        /// Number of correctable singleton errors generated per trusted root
        #[arg(long, default_value_t = 2)]
        correctable_singletons_per_trusted: usize,

        /// Number of protected singleton neighbors generated per weak root
        #[arg(long, default_value_t = 1)]
        protected_singletons_per_weak: usize,

        /// Starting RNG seed for task_000
        #[arg(long, default_value_t = 7)]
        seed: u64,

        /// Increment applied to the seed for each subsequent task
        #[arg(long, default_value_t = 1)]
        seed_step: u64,
    },

    /// Generate a panel of synthetic contig-extraction task directories
    PrepareContigExtractionPanel {
        /// Output directory for the task panel
        #[arg(long)]
        output: String,

        /// Number of task directories to create
        #[arg(long, default_value_t = 16)]
        tasks: usize,

        /// Scenario profile to generate: mixed, branching, repeat_fallback, repeat_completion, or repeat_priority
        #[arg(
            long,
            default_value = "mixed",
            value_parser = [
                "mixed",
                "branching",
                "repeat_fallback",
                "repeat_completion",
                "repeat_priority",
            ]
        )]
        profile: String,

        /// K-mer size for the synthetic graph
        #[arg(long, default_value_t = 11)]
        k: usize,

        /// Number of branching components for branching-profile tasks
        #[arg(long, default_value_t = 6)]
        component_count: usize,

        /// Number of primary-path reads per branching component
        #[arg(long, default_value_t = 6)]
        primary_reads_per_component: usize,

        /// Number of alternate-path reads per branching component
        #[arg(long, default_value_t = 2)]
        alternate_reads_per_component: usize,

        /// Starting RNG seed for task_000
        #[arg(long, default_value_t = 7)]
        seed: u64,

        /// Increment applied to the seed for each subsequent task
        #[arg(long, default_value_t = 1)]
        seed_step: u64,
    },
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

        /// Target coverage threshold for normalization
        #[arg(long, default_value_t = 500)]
        coverage_target: usize,

        /// Maximum number of reads to process (for subsampling)
        #[arg(long, default_value_t = 5_000_000)]
        max_reads: usize,
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
        #[arg(long, default_value_t = 50)]
        min_len: usize,

        /// Number of threads
        #[arg(long, default_value_t = num_cpus::get())]
        threads: usize,

        /// Output GFA format as well
        #[arg(long)]
        gfa: bool,

        /// Enable GPU acceleration for k-mer counting and overlap detection
        #[arg(long)]
        gpu: bool,

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
        #[arg(long, default_value_t = 0.9)]
        min_confidence: f64,

        /// Minimum path length for isoform traversal (shorter paths can be kept if highly confident)
        #[arg(long, default_value_t = 50)]
        min_path_len: usize,

        /// Enable development mode (bypasses filters and outputs debug information)
        #[arg(long)]
        dev_mode: bool,

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

    /// Evaluate or generate isolated component-level tasks for optimizer loops
    ComponentBench {
        #[command(subcommand)]
        command: ComponentBenchCommands,
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

    /// Assemble large genomes (100+ Gb) using disk-based k-mer counting
    AssembleLarge {
        /// Input FASTQ file(s) - use comma-separated paths for paired-end
        #[arg(short, long)]
        input: String,

        /// Second input FASTQ for paired-end reads (optional)
        #[arg(long)]
        input2: Option<String>,

        /// Output FASTA file
        #[arg(short, long)]
        output: String,

        /// K-mer size (recommend 31 for large genomes)
        #[arg(short, long, default_value_t = 31)]
        kmer: usize,

        /// Minimum k-mer count (filters sequencing errors, 0 = auto)
        #[arg(short = 'c', long, default_value_t = 0)]
        min_count: u32,

        /// Minimum count floor for trusted singleton-rescue neighbors
        #[arg(long, default_value_t = 4)]
        error_correction_min_trusted_count: u32,

        /// Minimum contig length to output
        #[arg(long, default_value_t = 200)]
        min_contig: usize,

        /// Number of threads
        #[arg(short, long, default_value_t = num_cpus::get())]
        threads: usize,

        /// Temporary directory for disk buckets
        #[arg(long)]
        temp_dir: Option<String>,

        /// Number of disk buckets (default: auto based on k-mer size)
        #[arg(long)]
        num_buckets: Option<usize>,

        /// Maximum tip length to remove during graph cleaning
        #[arg(long, default_value_t = 100)]
        max_tip_len: usize,

        /// Maximum bubble length to pop during graph cleaning
        #[arg(long, default_value_t = 50)]
        max_bubble_len: usize,

        /// Minimum read support required before branch threading overrides coverage
        #[arg(long, default_value_t = 1)]
        branch_support_min_win: u32,

        /// Minimum lead over the runner-up support before branch threading overrides coverage
        #[arg(long, default_value_t = 2)]
        branch_support_min_margin: u32,

        /// Disable the non-repeat preference in coverage-driven branch choice
        #[arg(long)]
        disable_prefer_non_repeat: bool,

        /// Disable prioritizing higher-coverage seeds first during contig extraction
        #[arg(long)]
        disable_prefer_high_count_seeds: bool,

        /// Disable prioritizing non-repeat seeds before repeat seeds during contig extraction
        #[arg(long)]
        disable_prefer_non_repeat_seeds: bool,

        /// Disable the repeat-seed completion pass after the primary extraction pass
        #[arg(long)]
        disable_repeat_seed_completion: bool,

        /// Suppress enclosed leftover branch-path contigs after extraction
        #[arg(long)]
        suppress_redundant_contigs: bool,

        /// Enable scaffolding using paired-end information
        #[arg(long)]
        scaffold: bool,

        /// Minimum number of read pairs to support a scaffold link
        #[arg(long, default_value_t = 3)]
        min_scaffold_links: usize,

        /// Enable contig polishing
        #[arg(long)]
        polish: bool,

        /// Number of polishing iterations
        #[arg(long, default_value_t = 1)]
        polish_iterations: usize,

        /// Minimizer k-mer size for post-assembly read mapping
        #[arg(long, default_value_t = 15)]
        read_mapping_k: usize,

        /// Minimizer window size for post-assembly read mapping
        #[arg(long, default_value_t = 10)]
        read_mapping_w: usize,

        /// Minimum supporting minimizer hits for a primary post-assembly mapping
        #[arg(long, default_value_t = 3)]
        read_mapping_min_primary_matches: usize,

        /// Minimum supporting minimizer hits for a scaffold-level contig hit
        #[arg(long, default_value_t = 2)]
        read_mapping_min_scaffold_matches: usize,

        /// Long reads file (FASTQ) for hybrid assembly
        #[arg(long)]
        long_reads: Option<String>,

        /// Minimum long read length to use
        #[arg(long, default_value_t = 1000)]
        min_long_read_len: usize,

        /// Enable LZ4 compression for disk buckets (faster I/O)
        #[arg(long)]
        compress_buckets: bool,
    },
}
