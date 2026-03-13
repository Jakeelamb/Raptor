use crate::io::fasta::FastaWriter;
use crate::io::fastq::{stream_fastq_records_checked, try_open_fastq, FastqRecord, FastqWriter};
use crate::kmer::kmer::KmerU64;
use crate::pipeline::large_genome_assembler::{
    LargeGenomeAssembler, LargeGenomeBenchBranchCandidate, LargeGenomeBenchBranchResolutionCase,
    LargeGenomeConfig, LargeGenomeStageBenchFixture,
};
use crate::pipeline::polisher::{read_contigs, MinimizerIndex, ReadMappingConfig};
use crate::pipeline::scaffolder::scaffold_and_polish_contigs;
use ahash::{AHashMap, AHashSet};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;
use tempfile::TempDir;

const DEFAULT_MINIMIZER_K: usize = 15;
const DEFAULT_MINIMIZER_W: usize = 10;
const DEFAULT_BRANCH_KMER_LEN: usize = 5;

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct ReadMappingTaskMetadata {
    component: String,
    name: String,
    description: String,
    paired: bool,
    minimizer_k: usize,
    minimizer_w: usize,
    num_contigs: usize,
    contig_len: usize,
    read_len: usize,
    reads_generated: usize,
    decoy_rate: f64,
    ambiguous_repeat_decoy_rate: f64,
    repeat_len: usize,
    error_rate: f64,
    seed: u64,
}

#[derive(Debug, Clone)]
struct ExpectedPrimaryHit {
    contig_id: usize,
    start: usize,
    is_reverse: bool,
}

#[derive(Debug, Clone)]
struct ExpectedScaffoldHit {
    contig_id: usize,
    is_reverse: bool,
}

#[derive(Debug, Clone, Default)]
struct ReadExpectation {
    primary: Option<ExpectedPrimaryHit>,
    scaffold: Option<ExpectedScaffoldHit>,
}

#[derive(Debug, Clone, Default)]
struct MappingCounters {
    reads_evaluated: u64,
    primary_expected: u64,
    primary_mapped: u64,
    primary_exact: u64,
    primary_near: u64,
    primary_contig_correct: u64,
    primary_orientation_correct: u64,
    primary_unexpected_total: u64,
    primary_unexpected_mapped: u64,
    scaffold_expected: u64,
    scaffold_mapped: u64,
    scaffold_exact: u64,
    scaffold_unexpected_total: u64,
    scaffold_unexpected_mapped: u64,
    primary_position_abs_error_sum: u64,
    primary_position_abs_error_count: u64,
}

#[derive(Debug, Clone, Default)]
struct PairMappingCounters {
    pair_expected: u64,
    pair_mapped: u64,
    pair_cross_contig_correct: u64,
    pair_unexpected_total: u64,
    pair_unexpected_mapped: u64,
    pair_cross_contig_errors: u64,
}

impl MappingCounters {
    fn observe(
        &mut self,
        expected: &ReadExpectation,
        primary_hit: Option<(usize, usize, bool)>,
        scaffold_hit: Option<(usize, usize, bool)>,
        position_tolerance: usize,
    ) {
        self.reads_evaluated += 1;

        match (&expected.primary, primary_hit) {
            (Some(primary), Some((contig_id, start, is_reverse))) => {
                self.primary_expected += 1;
                self.primary_mapped += 1;

                if contig_id == primary.contig_id {
                    self.primary_contig_correct += 1;
                }
                if is_reverse == primary.is_reverse {
                    self.primary_orientation_correct += 1;
                }
                if contig_id == primary.contig_id && is_reverse == primary.is_reverse {
                    let abs_error = start.abs_diff(primary.start) as u64;
                    self.primary_position_abs_error_sum += abs_error;
                    self.primary_position_abs_error_count += 1;
                    if abs_error == 0 {
                        self.primary_exact += 1;
                    }
                    if abs_error <= position_tolerance as u64 {
                        self.primary_near += 1;
                    }
                }
            }
            (Some(_), None) => {
                self.primary_expected += 1;
            }
            (None, Some(_)) => {
                self.primary_unexpected_total += 1;
                self.primary_unexpected_mapped += 1;
            }
            (None, None) => {
                self.primary_unexpected_total += 1;
            }
        }

        match (&expected.scaffold, scaffold_hit) {
            (Some(scaffold), Some((contig_id, _support, is_reverse))) => {
                self.scaffold_expected += 1;
                self.scaffold_mapped += 1;
                if contig_id == scaffold.contig_id && is_reverse == scaffold.is_reverse {
                    self.scaffold_exact += 1;
                }
            }
            (Some(_), None) => {
                self.scaffold_expected += 1;
            }
            (None, Some(_)) => {
                self.scaffold_unexpected_total += 1;
                self.scaffold_unexpected_mapped += 1;
            }
            (None, None) => {
                self.scaffold_unexpected_total += 1;
            }
        }
    }

    fn merge(&mut self, other: &Self) {
        self.reads_evaluated += other.reads_evaluated;
        self.primary_expected += other.primary_expected;
        self.primary_mapped += other.primary_mapped;
        self.primary_exact += other.primary_exact;
        self.primary_near += other.primary_near;
        self.primary_contig_correct += other.primary_contig_correct;
        self.primary_orientation_correct += other.primary_orientation_correct;
        self.primary_unexpected_total += other.primary_unexpected_total;
        self.primary_unexpected_mapped += other.primary_unexpected_mapped;
        self.scaffold_expected += other.scaffold_expected;
        self.scaffold_mapped += other.scaffold_mapped;
        self.scaffold_exact += other.scaffold_exact;
        self.scaffold_unexpected_total += other.scaffold_unexpected_total;
        self.scaffold_unexpected_mapped += other.scaffold_unexpected_mapped;
        self.primary_position_abs_error_sum += other.primary_position_abs_error_sum;
        self.primary_position_abs_error_count += other.primary_position_abs_error_count;
    }
}

impl PairMappingCounters {
    fn observe(&mut self, observation: &PairObservation) {
        match (
            &observation.mate1_expected,
            &observation.mate2_expected,
            observation.mate1_scaffold_hit,
            observation.mate2_scaffold_hit,
        ) {
            (Some(expected1), Some(expected2), Some((contig1, _)), Some((contig2, _))) => {
                self.pair_expected += 1;
                self.pair_mapped += 1;

                let expected_cross_contig = expected1.contig_id != expected2.contig_id;
                let observed_cross_contig = contig1 != contig2;
                if expected_cross_contig == observed_cross_contig {
                    self.pair_cross_contig_correct += 1;
                }
                if !expected_cross_contig && observed_cross_contig {
                    self.pair_cross_contig_errors += 1;
                }
            }
            (Some(_), Some(_), _, _) => {
                self.pair_expected += 1;
            }
            (None, None, Some(_), Some(_)) => {
                self.pair_unexpected_total += 1;
                self.pair_unexpected_mapped += 1;
            }
            (None, None, _, _) => {
                self.pair_unexpected_total += 1;
            }
            _ => {}
        }
    }

    fn merge(&mut self, other: &Self) {
        self.pair_expected += other.pair_expected;
        self.pair_mapped += other.pair_mapped;
        self.pair_cross_contig_correct += other.pair_cross_contig_correct;
        self.pair_unexpected_total += other.pair_unexpected_total;
        self.pair_unexpected_mapped += other.pair_unexpected_mapped;
        self.pair_cross_contig_errors += other.pair_cross_contig_errors;
    }
}

#[derive(Debug, Clone, Default)]
struct ObservedReadMapping {
    scaffold_hit: Option<(usize, bool)>,
}

#[derive(Debug, Clone, Default)]
struct PairObservation {
    mate1_expected: Option<ExpectedScaffoldHit>,
    mate2_expected: Option<ExpectedScaffoldHit>,
    mate1_scaffold_hit: Option<(usize, bool)>,
    mate2_scaffold_hit: Option<(usize, bool)>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ReadMappingTaskReport {
    pub task_dir: String,
    pub task_name: String,
    pub reads_evaluated: u64,
    pub primary_expected: u64,
    pub primary_mapped: u64,
    pub primary_exact: u64,
    pub primary_near: u64,
    pub primary_contig_correct: u64,
    pub primary_orientation_correct: u64,
    pub primary_unexpected_total: u64,
    pub primary_unexpected_mapped: u64,
    pub scaffold_expected: u64,
    pub scaffold_mapped: u64,
    pub scaffold_exact: u64,
    pub scaffold_unexpected_total: u64,
    pub scaffold_unexpected_mapped: u64,
    pub pair_expected: u64,
    pub pair_mapped: u64,
    pub pair_cross_contig_correct: u64,
    pub pair_unexpected_total: u64,
    pub pair_unexpected_mapped: u64,
    pub pair_cross_contig_errors: u64,
    pub primary_mapped_rate: Option<f64>,
    pub primary_exact_rate: Option<f64>,
    pub primary_near_rate: Option<f64>,
    pub primary_contig_rate: Option<f64>,
    pub primary_orientation_rate: Option<f64>,
    pub primary_unexpected_rate: Option<f64>,
    pub primary_specificity_rate: Option<f64>,
    pub scaffold_mapped_rate: Option<f64>,
    pub scaffold_exact_rate: Option<f64>,
    pub scaffold_unexpected_rate: Option<f64>,
    pub scaffold_specificity_rate: Option<f64>,
    pub pair_mapped_rate: Option<f64>,
    pub pair_cross_contig_rate: Option<f64>,
    pub pair_unexpected_rate: Option<f64>,
    pub pair_specificity_rate: Option<f64>,
    pub primary_position_mae: Option<f64>,
    pub primary_position_samples: u64,
    pub runtime_seconds: f64,
    pub reads_per_second: f64,
    pub score: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ReadMappingAggregateReport {
    pub task_count: usize,
    pub reads_evaluated: u64,
    pub primary_expected: u64,
    pub primary_mapped: u64,
    pub primary_exact: u64,
    pub primary_near: u64,
    pub primary_contig_correct: u64,
    pub primary_orientation_correct: u64,
    pub primary_unexpected_total: u64,
    pub primary_unexpected_mapped: u64,
    pub scaffold_expected: u64,
    pub scaffold_mapped: u64,
    pub scaffold_exact: u64,
    pub scaffold_unexpected_total: u64,
    pub scaffold_unexpected_mapped: u64,
    pub pair_expected: u64,
    pub pair_mapped: u64,
    pub pair_cross_contig_correct: u64,
    pub pair_unexpected_total: u64,
    pub pair_unexpected_mapped: u64,
    pub pair_cross_contig_errors: u64,
    pub primary_mapped_rate: Option<f64>,
    pub primary_exact_rate: Option<f64>,
    pub primary_near_rate: Option<f64>,
    pub primary_contig_rate: Option<f64>,
    pub primary_orientation_rate: Option<f64>,
    pub primary_unexpected_rate: Option<f64>,
    pub primary_specificity_rate: Option<f64>,
    pub scaffold_mapped_rate: Option<f64>,
    pub scaffold_exact_rate: Option<f64>,
    pub scaffold_unexpected_rate: Option<f64>,
    pub scaffold_specificity_rate: Option<f64>,
    pub pair_mapped_rate: Option<f64>,
    pub pair_cross_contig_rate: Option<f64>,
    pub pair_unexpected_rate: Option<f64>,
    pub pair_specificity_rate: Option<f64>,
    pub primary_position_mae: Option<f64>,
    pub primary_position_samples: u64,
    pub runtime_seconds: f64,
    pub reads_per_second: f64,
    pub score: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ReadMappingBenchmarkReport {
    pub component: &'static str,
    pub minimizer_k: usize,
    pub minimizer_w: usize,
    pub min_primary_matches: usize,
    pub min_scaffold_matches: usize,
    pub position_tolerance: usize,
    pub tasks: Vec<ReadMappingTaskReport>,
    pub summary: ReadMappingAggregateReport,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct BranchResolutionTaskMetadata {
    component: String,
    name: String,
    description: String,
    kmer_len: usize,
    cases_generated: usize,
    seed: u64,
    scenario_profile: BranchScenarioProfile,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
enum BranchScenarioProfile {
    #[default]
    Balanced,
    SupportStress,
}

impl BranchScenarioProfile {
    fn parse(raw: &str) -> io::Result<Self> {
        match raw {
            "balanced" => Ok(Self::Balanced),
            "support_stress" | "support-stress" => Ok(Self::SupportStress),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "unknown branch scenario profile '{raw}'; expected 'balanced' or 'support_stress'"
                ),
            )),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct BranchResolutionCandidateSpec {
    base_idx: usize,
    next_kmer: String,
    count: u32,
    is_repeat: bool,
    read_support: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct BranchResolutionCaseSpec {
    scenario: String,
    current_count: u32,
    expected_base_idx: usize,
    expected_next_kmer: String,
    candidates: Vec<BranchResolutionCandidateSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct BranchResolutionTaskFile {
    cases: Vec<BranchResolutionCaseSpec>,
}

#[derive(Debug, Clone, Default)]
struct BranchResolutionCounters {
    cases_evaluated: u64,
    cases_with_choice: u64,
    exact_matches: u64,
    scenario_counts: AHashMap<String, u64>,
    scenario_correct: AHashMap<String, u64>,
}

impl BranchResolutionCounters {
    fn observe(&mut self, scenario: &str, matched: bool, chose_any: bool) {
        self.cases_evaluated += 1;
        *self
            .scenario_counts
            .entry(scenario.to_string())
            .or_default() += 1;
        if chose_any {
            self.cases_with_choice += 1;
        }
        if matched {
            self.exact_matches += 1;
            *self
                .scenario_correct
                .entry(scenario.to_string())
                .or_default() += 1;
        }
    }

    fn merge(&mut self, other: &Self) {
        self.cases_evaluated += other.cases_evaluated;
        self.cases_with_choice += other.cases_with_choice;
        self.exact_matches += other.exact_matches;
        for (scenario, count) in &other.scenario_counts {
            *self.scenario_counts.entry(scenario.clone()).or_default() += count;
        }
        for (scenario, count) in &other.scenario_correct {
            *self.scenario_correct.entry(scenario.clone()).or_default() += count;
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct BranchResolutionScenarioMetric {
    pub cases_evaluated: u64,
    pub exact_matches: u64,
    pub exact_rate: Option<f64>,
}

#[derive(Debug, Clone, Serialize)]
pub struct BranchResolutionTaskReport {
    pub task_dir: String,
    pub task_name: String,
    pub cases_evaluated: u64,
    pub cases_with_choice: u64,
    pub exact_matches: u64,
    pub choice_rate: Option<f64>,
    pub exact_rate: Option<f64>,
    pub runtime_seconds: f64,
    pub cases_per_second: f64,
    pub score: f64,
    pub scenario_metrics: BTreeMap<String, BranchResolutionScenarioMetric>,
}

#[derive(Debug, Clone, Serialize)]
pub struct BranchResolutionAggregateReport {
    pub task_count: usize,
    pub cases_evaluated: u64,
    pub cases_with_choice: u64,
    pub exact_matches: u64,
    pub choice_rate: Option<f64>,
    pub exact_rate: Option<f64>,
    pub runtime_seconds: f64,
    pub cases_per_second: f64,
    pub score: f64,
    pub scenario_metrics: BTreeMap<String, BranchResolutionScenarioMetric>,
}

#[derive(Debug, Clone, Serialize)]
pub struct BranchResolutionBenchmarkReport {
    pub component: &'static str,
    pub branch_support_min_win: u32,
    pub branch_support_min_margin: u32,
    pub prefer_non_repeat_branches: bool,
    pub tasks: Vec<BranchResolutionTaskReport>,
    pub summary: BranchResolutionAggregateReport,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct ScaffoldPolishTaskMetadata {
    component: String,
    name: String,
    description: String,
    scaffold_groups: usize,
    contigs_per_scaffold: usize,
    contigs_generated: usize,
    contig_len: usize,
    read_len: usize,
    insert_size: usize,
    internal_pairs_per_contig: usize,
    true_link_pairs: usize,
    decoy_link_pairs: usize,
    mutations_per_contig: usize,
    pairs_generated: usize,
    seed: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
struct ScaffoldPolishPlacementSpec {
    contig: String,
    is_reverse: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ScaffoldPolishTruthScaffold {
    contigs: Vec<ScaffoldPolishPlacementSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ScaffoldPolishTruthContig {
    header: String,
    sequence: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ScaffoldPolishTruthFile {
    scaffolds: Vec<ScaffoldPolishTruthScaffold>,
    polished_contigs: Vec<ScaffoldPolishTruthContig>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ScaffoldPolishTaskReport {
    pub task_dir: String,
    pub task_name: String,
    pub pairs_evaluated: u64,
    pub expected_links: u64,
    pub observed_links: u64,
    pub correct_links: u64,
    pub expected_scaffolds: u64,
    pub observed_scaffolds: u64,
    pub exact_scaffolds: u64,
    pub polished_contigs_expected: u64,
    pub polished_contigs_observed: u64,
    pub polished_contigs_exact: u64,
    pub polished_bases_expected: u64,
    pub polished_bases_correct: u64,
    pub scaffold_link_precision: Option<f64>,
    pub scaffold_link_recall: Option<f64>,
    pub scaffold_link_f1: Option<f64>,
    pub scaffold_exact_rate: Option<f64>,
    pub polished_exact_rate: Option<f64>,
    pub polished_base_accuracy: Option<f64>,
    pub runtime_seconds: f64,
    pub pairs_per_second: f64,
    pub score: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ScaffoldPolishAggregateReport {
    pub task_count: usize,
    pub pairs_evaluated: u64,
    pub expected_links: u64,
    pub observed_links: u64,
    pub correct_links: u64,
    pub expected_scaffolds: u64,
    pub observed_scaffolds: u64,
    pub exact_scaffolds: u64,
    pub polished_contigs_expected: u64,
    pub polished_contigs_observed: u64,
    pub polished_contigs_exact: u64,
    pub polished_bases_expected: u64,
    pub polished_bases_correct: u64,
    pub scaffold_link_precision: Option<f64>,
    pub scaffold_link_recall: Option<f64>,
    pub scaffold_link_f1: Option<f64>,
    pub scaffold_exact_rate: Option<f64>,
    pub polished_exact_rate: Option<f64>,
    pub polished_base_accuracy: Option<f64>,
    pub runtime_seconds: f64,
    pub pairs_per_second: f64,
    pub score: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ScaffoldPolishBenchmarkReport {
    pub component: &'static str,
    pub min_scaffold_links: usize,
    pub min_primary_matches: usize,
    pub min_scaffold_matches: usize,
    pub tasks: Vec<ScaffoldPolishTaskReport>,
    pub summary: ScaffoldPolishAggregateReport,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct ErrorCorrectionTaskMetadata {
    component: String,
    name: String,
    description: String,
    k: usize,
    trusted_roots: usize,
    weak_roots: usize,
    correctable_singletons_per_trusted: usize,
    protected_singletons_per_weak: usize,
    kmers_generated: usize,
    seed: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ErrorCorrectionCountSpec {
    kmer: String,
    count: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ErrorCorrectionTaskFile {
    counts: Vec<ErrorCorrectionCountSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ErrorCorrectionTruthFile {
    expected_removed_singletons: Vec<String>,
    expected_preserved_singletons: Vec<String>,
    expected_trusted_totals: BTreeMap<String, u32>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ErrorCorrectionTaskReport {
    pub task_dir: String,
    pub task_name: String,
    pub input_kmers: u64,
    pub input_singletons: u64,
    pub corrected_kmers: u64,
    pub remaining_kmers: u64,
    pub expected_removed_singletons: u64,
    pub observed_removed_singletons: u64,
    pub correctly_removed_singletons: u64,
    pub expected_preserved_singletons: u64,
    pub preserved_singletons_retained: u64,
    pub trusted_targets: u64,
    pub trusted_targets_exact: u64,
    pub correction_precision: Option<f64>,
    pub correction_recall: Option<f64>,
    pub correction_f1: Option<f64>,
    pub preserved_retention_rate: Option<f64>,
    pub trusted_target_exact_rate: Option<f64>,
    pub runtime_seconds: f64,
    pub kmers_per_second: f64,
    pub score: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ErrorCorrectionAggregateReport {
    pub task_count: usize,
    pub input_kmers: u64,
    pub input_singletons: u64,
    pub corrected_kmers: u64,
    pub remaining_kmers: u64,
    pub expected_removed_singletons: u64,
    pub observed_removed_singletons: u64,
    pub correctly_removed_singletons: u64,
    pub expected_preserved_singletons: u64,
    pub preserved_singletons_retained: u64,
    pub trusted_targets: u64,
    pub trusted_targets_exact: u64,
    pub correction_precision: Option<f64>,
    pub correction_recall: Option<f64>,
    pub correction_f1: Option<f64>,
    pub preserved_retention_rate: Option<f64>,
    pub trusted_target_exact_rate: Option<f64>,
    pub runtime_seconds: f64,
    pub kmers_per_second: f64,
    pub score: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ErrorCorrectionBenchmarkReport {
    pub component: &'static str,
    pub min_count: u32,
    pub min_trusted_count: u32,
    pub tasks: Vec<ErrorCorrectionTaskReport>,
    pub summary: ErrorCorrectionAggregateReport,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct ContigExtractionTaskMetadata {
    component: String,
    name: String,
    description: String,
    k: usize,
    profile: ContigExtractionScenarioProfile,
    component_count: usize,
    primary_reads_per_component: usize,
    alternate_reads_per_component: usize,
    expected_contigs: usize,
    seed: u64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
enum ContigExtractionScenarioProfile {
    #[default]
    Branching,
    RepeatFallback,
    RepeatCompletion,
    RepeatPriority,
}

impl ContigExtractionScenarioProfile {
    fn parse(raw: &str) -> io::Result<Self> {
        match raw {
            "branching" => Ok(Self::Branching),
            "repeat_fallback" | "repeat-fallback" => Ok(Self::RepeatFallback),
            "repeat_completion" | "repeat-completion" => Ok(Self::RepeatCompletion),
            "repeat_priority" | "repeat-priority" => Ok(Self::RepeatPriority),
            _ => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "unknown contig extraction profile '{raw}'; expected 'branching', 'repeat_fallback', 'repeat_completion', or 'repeat_priority'"
                ),
            )),
        }
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Branching => "branching",
            Self::RepeatFallback => "repeat_fallback",
            Self::RepeatCompletion => "repeat_completion",
            Self::RepeatPriority => "repeat_priority",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ContigExtractionCountSpec {
    kmer: String,
    count: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ContigExtractionBranchSupportSpec {
    from_kmer: String,
    to_kmer: String,
    support: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ContigExtractionGraphFile {
    valid_kmers: Vec<String>,
    branch_support: Vec<ContigExtractionBranchSupportSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ContigExtractionTruthFile {
    expected_contigs: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ContigExtractionTaskReport {
    pub task_dir: String,
    pub task_name: String,
    pub profile: String,
    pub graph_nodes: u64,
    pub expected_contigs: u64,
    pub observed_contigs: u64,
    pub exact_contigs: u64,
    pub truth_kmers_expected: u64,
    pub truth_kmers_observed: u64,
    pub truth_kmers_correct: u64,
    pub exact_contig_rate: Option<f64>,
    pub contig_count_agreement: Option<f64>,
    pub truth_kmer_precision: Option<f64>,
    pub truth_kmer_recall: Option<f64>,
    pub truth_kmer_f1: Option<f64>,
    pub runtime_seconds: f64,
    pub graph_nodes_per_second: f64,
    pub score: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ContigExtractionAggregateReport {
    pub task_count: usize,
    pub graph_nodes: u64,
    pub expected_contigs: u64,
    pub observed_contigs: u64,
    pub exact_contigs: u64,
    pub truth_kmers_expected: u64,
    pub truth_kmers_observed: u64,
    pub truth_kmers_correct: u64,
    pub exact_contig_rate: Option<f64>,
    pub contig_count_agreement: Option<f64>,
    pub truth_kmer_precision: Option<f64>,
    pub truth_kmer_recall: Option<f64>,
    pub truth_kmer_f1: Option<f64>,
    pub runtime_seconds: f64,
    pub graph_nodes_per_second: f64,
    pub score: f64,
}

#[derive(Debug, Clone, Serialize)]
pub struct ContigExtractionBenchmarkReport {
    pub component: &'static str,
    pub prefer_high_count_seeds: bool,
    pub prefer_non_repeat_seeds: bool,
    pub enable_repeat_seed_completion: bool,
    pub suppress_redundant_contigs: bool,
    pub tasks: Vec<ContigExtractionTaskReport>,
    pub summary: ContigExtractionAggregateReport,
}

pub fn run_contig_extraction_bench(
    task_path: &str,
    prefer_high_count_seeds: bool,
    prefer_non_repeat_seeds: bool,
    enable_repeat_seed_completion: bool,
    suppress_redundant_contigs: bool,
    json: bool,
    output: Option<&str>,
) -> io::Result<()> {
    let report = evaluate_contig_extraction_bench(
        Path::new(task_path),
        prefer_high_count_seeds,
        prefer_non_repeat_seeds,
        enable_repeat_seed_completion,
        suppress_redundant_contigs,
    )?;

    if let Some(output_path) = output {
        let serialized = serde_json::to_string_pretty(&report)
            .map_err(|err| io::Error::other(format!("serialize report: {err}")))?;
        if let Some(parent) = Path::new(output_path)
            .parent()
            .filter(|path| !path.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }
        fs::write(output_path, serialized)?;
    }

    if json {
        let serialized = serde_json::to_string_pretty(&report)
            .map_err(|err| io::Error::other(format!("serialize report: {err}")))?;
        println!("{serialized}");
    } else {
        print_contig_extraction_report(&report);
    }

    Ok(())
}

pub fn run_read_mapping_bench(
    task_path: &str,
    minimizer_k: Option<usize>,
    minimizer_w: Option<usize>,
    min_primary_matches: usize,
    min_scaffold_matches: usize,
    position_tolerance: usize,
    json: bool,
    output: Option<&str>,
) -> io::Result<()> {
    let report = evaluate_read_mapping_bench(
        Path::new(task_path),
        minimizer_k.unwrap_or(DEFAULT_MINIMIZER_K),
        minimizer_w.unwrap_or(DEFAULT_MINIMIZER_W),
        min_primary_matches,
        min_scaffold_matches,
        position_tolerance,
    )?;

    if let Some(output_path) = output {
        let serialized = serde_json::to_string_pretty(&report)
            .map_err(|err| io::Error::other(format!("serialize report: {err}")))?;
        if let Some(parent) = Path::new(output_path)
            .parent()
            .filter(|path| !path.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }
        fs::write(output_path, serialized)?;
    }

    if json {
        let serialized = serde_json::to_string_pretty(&report)
            .map_err(|err| io::Error::other(format!("serialize report: {err}")))?;
        println!("{serialized}");
    } else {
        print_read_mapping_report(&report);
    }

    Ok(())
}

pub fn prepare_read_mapping_task(
    output_dir: &str,
    paired: bool,
    num_contigs: usize,
    contig_len: usize,
    read_len: usize,
    reads: usize,
    insert_size: usize,
    repeat_len: usize,
    error_rate: f64,
    decoy_rate: f64,
    ambiguous_repeat_decoy_rate: f64,
    seed: u64,
    minimizer_k: usize,
    minimizer_w: usize,
) -> io::Result<()> {
    validate_prepare_args(
        paired,
        num_contigs,
        contig_len,
        read_len,
        reads,
        insert_size,
        repeat_len,
        error_rate,
        decoy_rate,
        ambiguous_repeat_decoy_rate,
    )?;

    let output_dir = Path::new(output_dir);
    fs::create_dir_all(output_dir)?;

    let mut rng = StdRng::seed_from_u64(seed);
    let repeat = random_dna(repeat_len, &mut rng);
    let contigs = build_synthetic_contigs(num_contigs, contig_len, &repeat, &mut rng);

    let contigs_path = output_dir.join("contigs.fa");
    let reads1_path = output_dir.join("reads_1.fastq");
    let reads2_path = paired.then(|| output_dir.join("reads_2.fastq"));
    let truth_path = output_dir.join("truth.tsv");
    let metadata_path = output_dir.join("task.json");

    write_contigs_fasta(&contigs_path, &contigs)?;
    write_synthetic_reads(
        &contigs,
        &repeat,
        read_len,
        reads,
        paired,
        insert_size,
        error_rate,
        decoy_rate,
        ambiguous_repeat_decoy_rate,
        &reads1_path,
        reads2_path.as_deref(),
        &truth_path,
        &mut rng,
    )?;

    let metadata = ReadMappingTaskMetadata {
        component: "read_mapping".to_string(),
        name: output_dir
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("read_mapping_task")
            .to_string(),
        description: "Synthetic mapper task with repeated sequence and planted truth.".to_string(),
        paired,
        minimizer_k,
        minimizer_w,
        num_contigs,
        contig_len,
        read_len,
        reads_generated: total_generated_reads(
            reads,
            paired,
            decoy_rate,
            ambiguous_repeat_decoy_rate,
        ),
        decoy_rate,
        ambiguous_repeat_decoy_rate,
        repeat_len,
        error_rate,
        seed,
    };
    let metadata_json = serde_json::to_string_pretty(&metadata)
        .map_err(|err| io::Error::other(format!("serialize metadata: {err}")))?;
    fs::write(metadata_path, metadata_json)?;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn prepare_read_mapping_panel(
    output_dir: &str,
    tasks: usize,
    paired: bool,
    num_contigs: usize,
    contig_len: usize,
    read_len: usize,
    reads: usize,
    insert_size: usize,
    repeat_len: usize,
    error_rate: f64,
    decoy_rate: f64,
    ambiguous_repeat_decoy_rate: f64,
    seed: u64,
    seed_step: u64,
    minimizer_k: usize,
    minimizer_w: usize,
) -> io::Result<()> {
    if tasks == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "tasks must be greater than zero",
        ));
    }

    let output_dir = Path::new(output_dir);
    fs::create_dir_all(output_dir)?;

    for task_index in 0..tasks {
        let task_seed = seed.saturating_add(seed_step.saturating_mul(task_index as u64));
        let task_dir = output_dir.join(format!("task_{task_index:03}"));
        prepare_read_mapping_task(
            path_to_str(&task_dir)?,
            paired,
            num_contigs,
            contig_len,
            read_len,
            reads,
            insert_size,
            repeat_len,
            error_rate,
            decoy_rate,
            ambiguous_repeat_decoy_rate,
            task_seed,
            minimizer_k,
            minimizer_w,
        )?;
    }

    Ok(())
}

pub fn run_branch_resolution_bench(
    task_path: &str,
    branch_support_min_win: u32,
    branch_support_min_margin: u32,
    prefer_non_repeat_branches: bool,
    json: bool,
    output: Option<&str>,
) -> io::Result<()> {
    let report = evaluate_branch_resolution_bench(
        Path::new(task_path),
        branch_support_min_win,
        branch_support_min_margin,
        prefer_non_repeat_branches,
    )?;

    if let Some(output_path) = output {
        let serialized = serde_json::to_string_pretty(&report)
            .map_err(|err| io::Error::other(format!("serialize report: {err}")))?;
        if let Some(parent) = Path::new(output_path)
            .parent()
            .filter(|path| !path.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }
        fs::write(output_path, serialized)?;
    }

    if json {
        let serialized = serde_json::to_string_pretty(&report)
            .map_err(|err| io::Error::other(format!("serialize report: {err}")))?;
        println!("{serialized}");
    } else {
        print_branch_resolution_report(&report);
    }

    Ok(())
}

pub fn prepare_branch_resolution_task(output_dir: &str, cases: usize, seed: u64) -> io::Result<()> {
    prepare_branch_resolution_task_with_profile(output_dir, cases, seed, "balanced")
}

pub fn prepare_branch_resolution_task_with_profile(
    output_dir: &str,
    cases: usize,
    seed: u64,
    scenario_profile: &str,
) -> io::Result<()> {
    if cases == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "cases must be greater than zero",
        ));
    }

    let output_dir = Path::new(output_dir);
    fs::create_dir_all(output_dir)?;
    let scenario_profile = BranchScenarioProfile::parse(scenario_profile)?;

    let task_file = BranchResolutionTaskFile {
        cases: build_branch_resolution_cases(cases, seed, scenario_profile),
    };
    let metadata = BranchResolutionTaskMetadata {
        component: "branch_resolution".to_string(),
        name: output_dir
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("branch_resolution_task")
            .to_string(),
        description: "Synthetic branch-choice task panel targeting branch-selection heuristics."
            .to_string(),
        kmer_len: DEFAULT_BRANCH_KMER_LEN,
        cases_generated: task_file.cases.len(),
        seed,
        scenario_profile,
    };

    fs::write(
        output_dir.join("cases.json"),
        serde_json::to_string_pretty(&task_file)
            .map_err(|err| io::Error::other(format!("serialize cases: {err}")))?,
    )?;
    fs::write(
        output_dir.join("task.json"),
        serde_json::to_string_pretty(&metadata)
            .map_err(|err| io::Error::other(format!("serialize metadata: {err}")))?,
    )?;

    Ok(())
}

pub fn prepare_branch_resolution_panel(
    output_dir: &str,
    tasks: usize,
    cases_per_task: usize,
    seed: u64,
    seed_step: u64,
) -> io::Result<()> {
    prepare_branch_resolution_panel_with_profile(
        output_dir,
        tasks,
        cases_per_task,
        seed,
        seed_step,
        "balanced",
    )
}

pub fn prepare_branch_resolution_panel_with_profile(
    output_dir: &str,
    tasks: usize,
    cases_per_task: usize,
    seed: u64,
    seed_step: u64,
    scenario_profile: &str,
) -> io::Result<()> {
    if tasks == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "tasks must be greater than zero",
        ));
    }

    let output_dir = Path::new(output_dir);
    fs::create_dir_all(output_dir)?;

    for task_index in 0..tasks {
        let task_seed = seed.saturating_add(seed_step.saturating_mul(task_index as u64));
        let task_dir = output_dir.join(format!("task_{task_index:03}"));
        prepare_branch_resolution_task_with_profile(
            path_to_str(&task_dir)?,
            cases_per_task,
            task_seed,
            scenario_profile,
        )?;
    }

    Ok(())
}

pub fn run_scaffold_polish_bench(
    task_path: &str,
    min_scaffold_links: usize,
    min_primary_matches: usize,
    min_scaffold_matches: usize,
    json: bool,
    output: Option<&str>,
) -> io::Result<()> {
    let report = evaluate_scaffold_polish_bench(
        Path::new(task_path),
        min_scaffold_links,
        min_primary_matches,
        min_scaffold_matches,
    )?;

    if let Some(output_path) = output {
        let serialized = serde_json::to_string_pretty(&report)
            .map_err(|err| io::Error::other(format!("serialize report: {err}")))?;
        if let Some(parent) = Path::new(output_path)
            .parent()
            .filter(|path| !path.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }
        fs::write(output_path, serialized)?;
    }

    if json {
        let serialized = serde_json::to_string_pretty(&report)
            .map_err(|err| io::Error::other(format!("serialize report: {err}")))?;
        println!("{serialized}");
    } else {
        print_scaffold_polish_report(&report);
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn prepare_scaffold_polish_task(
    output_dir: &str,
    scaffold_groups: usize,
    contigs_per_scaffold: usize,
    contig_len: usize,
    read_len: usize,
    insert_size: usize,
    internal_pairs_per_contig: usize,
    true_link_pairs: usize,
    decoy_link_pairs: usize,
    mutations_per_contig: usize,
    seed: u64,
) -> io::Result<()> {
    validate_scaffold_polish_args(
        scaffold_groups,
        contigs_per_scaffold,
        contig_len,
        read_len,
        insert_size,
        internal_pairs_per_contig,
        true_link_pairs,
        decoy_link_pairs,
    )?;

    let output_dir = Path::new(output_dir);
    fs::create_dir_all(output_dir)?;

    let total_contigs = scaffold_groups * contigs_per_scaffold;
    let mut rng = StdRng::seed_from_u64(seed);
    let repeat = random_dna(read_len.clamp(24, 48), &mut rng);
    let correct_contigs = build_synthetic_contigs(total_contigs, contig_len, &repeat, &mut rng);
    let mutated_contigs =
        mutate_contigs_for_polishing(&correct_contigs, mutations_per_contig, &mut rng)?;
    let truth_scaffolds = build_truth_scaffolds(
        &correct_contigs,
        scaffold_groups,
        contigs_per_scaffold,
        &mut rng,
    );

    let contigs_path = output_dir.join("contigs.fa");
    let reads1_path = output_dir.join("reads_1.fastq");
    let reads2_path = output_dir.join("reads_2.fastq");
    let truth_path = output_dir.join("truth.json");
    let metadata_path = output_dir.join("task.json");

    write_contigs_fasta(&contigs_path, &mutated_contigs)?;
    let pairs_generated = write_scaffold_polish_reads(
        &correct_contigs,
        &truth_scaffolds,
        read_len,
        insert_size,
        internal_pairs_per_contig,
        true_link_pairs,
        decoy_link_pairs,
        &reads1_path,
        &reads2_path,
        &mut rng,
    )?;

    let truth = ScaffoldPolishTruthFile {
        scaffolds: truth_scaffolds,
        polished_contigs: correct_contigs
            .iter()
            .map(|(header, sequence)| {
                Ok(ScaffoldPolishTruthContig {
                    header: header.clone(),
                    sequence: std::str::from_utf8(sequence)
                        .map_err(|err| {
                            io::Error::new(
                                io::ErrorKind::InvalidData,
                                format!("truth contig sequence is not UTF-8 DNA: {err}"),
                            )
                        })?
                        .to_string(),
                })
            })
            .collect::<io::Result<Vec<_>>>()?,
    };
    fs::write(
        &truth_path,
        serde_json::to_string_pretty(&truth)
            .map_err(|err| io::Error::other(format!("serialize truth: {err}")))?,
    )?;

    let metadata = ScaffoldPolishTaskMetadata {
        component: "scaffold_polish".to_string(),
        name: output_dir
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("scaffold_polish_task")
            .to_string(),
        description: "Synthetic shared scaffold+polish task with true joins, decoy joins, and mutated contigs."
            .to_string(),
        scaffold_groups,
        contigs_per_scaffold,
        contigs_generated: total_contigs,
        contig_len,
        read_len,
        insert_size,
        internal_pairs_per_contig,
        true_link_pairs,
        decoy_link_pairs,
        mutations_per_contig,
        pairs_generated,
        seed,
    };
    fs::write(
        metadata_path,
        serde_json::to_string_pretty(&metadata)
            .map_err(|err| io::Error::other(format!("serialize metadata: {err}")))?,
    )?;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn prepare_scaffold_polish_panel(
    output_dir: &str,
    tasks: usize,
    scaffold_groups: usize,
    contigs_per_scaffold: usize,
    contig_len: usize,
    read_len: usize,
    insert_size: usize,
    internal_pairs_per_contig: usize,
    true_link_pairs: usize,
    decoy_link_pairs: usize,
    mutations_per_contig: usize,
    seed: u64,
    seed_step: u64,
) -> io::Result<()> {
    if tasks == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "tasks must be greater than zero",
        ));
    }

    let output_dir = Path::new(output_dir);
    fs::create_dir_all(output_dir)?;

    for task_index in 0..tasks {
        let task_seed = seed.saturating_add(seed_step.saturating_mul(task_index as u64));
        let task_dir = output_dir.join(format!("task_{task_index:03}"));
        prepare_scaffold_polish_task(
            path_to_str(&task_dir)?,
            scaffold_groups,
            contigs_per_scaffold,
            contig_len,
            read_len,
            insert_size,
            internal_pairs_per_contig,
            true_link_pairs,
            decoy_link_pairs,
            mutations_per_contig,
            task_seed,
        )?;
    }

    Ok(())
}

fn evaluate_scaffold_polish_bench(
    task_path: &Path,
    min_scaffold_links: usize,
    min_primary_matches: usize,
    min_scaffold_matches: usize,
) -> io::Result<ScaffoldPolishBenchmarkReport> {
    let task_dirs = discover_scaffold_polish_task_dirs(task_path)?;
    let tasks = task_dirs
        .par_iter()
        .map(|task_dir| {
            evaluate_single_scaffold_polish_task(
                task_dir,
                min_scaffold_links,
                min_primary_matches,
                min_scaffold_matches,
            )
        })
        .collect::<io::Result<Vec<_>>>()?;
    let summary = summarize_scaffold_polish_reports(&tasks);

    Ok(ScaffoldPolishBenchmarkReport {
        component: "scaffold_polish",
        min_scaffold_links,
        min_primary_matches,
        min_scaffold_matches,
        tasks,
        summary,
    })
}

fn discover_scaffold_polish_task_dirs(task_path: &Path) -> io::Result<Vec<PathBuf>> {
    if is_scaffold_polish_task_dir(task_path) {
        return Ok(vec![task_path.to_path_buf()]);
    }

    let mut discovered = Vec::new();
    for entry in fs::read_dir(task_path)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() && is_scaffold_polish_task_dir(&path) {
            discovered.push(path);
        }
    }
    discovered.sort();

    if discovered.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "no scaffold-polish tasks found under {}",
                task_path.display()
            ),
        ));
    }

    Ok(discovered)
}

fn is_scaffold_polish_task_dir(task_dir: &Path) -> bool {
    resolve_scaffold_polish_truth_path(task_dir).is_ok()
        && resolve_reads2_path(task_dir).is_some()
        && resolve_contigs_path(task_dir).is_ok()
}

fn evaluate_single_scaffold_polish_task(
    task_dir: &Path,
    min_scaffold_links: usize,
    min_primary_matches: usize,
    min_scaffold_matches: usize,
) -> io::Result<ScaffoldPolishTaskReport> {
    let task_name = read_scaffold_polish_task_metadata(task_dir)?
        .and_then(|task| {
            if task.name.trim().is_empty() {
                None
            } else {
                Some(task.name)
            }
        })
        .unwrap_or_else(|| {
            task_dir
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("scaffold_polish_task")
                .to_string()
        });
    let metadata = read_scaffold_polish_task_metadata(task_dir)?.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "missing scaffold-polish metadata under {}",
                task_dir.display()
            ),
        )
    })?;
    let truth = load_scaffold_polish_truth_file(task_dir)?;
    let input_contigs = read_contigs(path_to_str(&resolve_contigs_path(task_dir)?)?)?;

    let temp_dir = TempDir::new()?;
    let scaffold_output_path = temp_dir.path().join("scaffolds.fa");
    let polish_output_path = temp_dir.path().join("polished.fa");
    let mapping_config = ReadMappingConfig {
        min_primary_matches,
        min_scaffold_matches,
        ..ReadMappingConfig::default()
    };

    let start = Instant::now();
    scaffold_and_polish_contigs(
        path_to_str(&resolve_contigs_path(task_dir)?)?,
        path_to_str(&resolve_reads1_path(task_dir)?)?,
        path_to_str(&resolve_reads2_path(task_dir).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("{} is missing paired reads_2.fastq", task_dir.display()),
            )
        })?)?,
        path_to_str(&scaffold_output_path)?,
        path_to_str(&polish_output_path)?,
        min_scaffold_links,
        mapping_config,
    )?;
    let runtime_seconds = start.elapsed().as_secs_f64();

    let observed_scaffolds = decode_scaffold_output(&scaffold_output_path, &input_contigs)?;
    let polished_contigs = read_contigs(path_to_str(&polish_output_path)?)?;

    Ok(build_scaffold_polish_task_report(
        task_dir,
        task_name,
        metadata.pairs_generated as u64,
        &truth,
        &observed_scaffolds,
        &polished_contigs,
        runtime_seconds,
    ))
}

fn build_scaffold_polish_task_report(
    task_dir: &Path,
    task_name: String,
    pairs_evaluated: u64,
    truth: &ScaffoldPolishTruthFile,
    observed_scaffolds: &[Vec<ScaffoldPolishPlacementSpec>],
    polished_contigs: &[(String, Vec<u8>)],
    runtime_seconds: f64,
) -> ScaffoldPolishTaskReport {
    let expected_links = scaffold_link_set_from_truth(&truth.scaffolds);
    let observed_links = scaffold_link_set_from_observed(observed_scaffolds);
    let correct_links = expected_links.intersection(&observed_links).count() as u64;
    let expected_scaffold_set = canonical_truth_scaffold_set(&truth.scaffolds);
    let observed_scaffold_set = canonical_observed_scaffold_set(observed_scaffolds);
    let exact_scaffolds = expected_scaffold_set
        .intersection(&observed_scaffold_set)
        .count() as u64;
    let (polished_exact_contigs, polished_bases_correct, polished_bases_expected) =
        compare_polished_contigs(&truth.polished_contigs, polished_contigs);

    let scaffold_link_precision = ratio(correct_links, observed_links.len() as u64);
    let scaffold_link_recall = ratio(correct_links, expected_links.len() as u64);
    let scaffold_link_f1 = f1_score(
        scaffold_link_precision.unwrap_or(0.0),
        scaffold_link_recall.unwrap_or(0.0),
    );
    let scaffold_exact_rate = ratio(exact_scaffolds, truth.scaffolds.len() as u64);
    let polished_exact_rate = ratio(polished_exact_contigs, truth.polished_contigs.len() as u64);
    let polished_base_accuracy = ratio(polished_bases_correct, polished_bases_expected);
    let pairs_per_second = if runtime_seconds > 0.0 {
        pairs_evaluated as f64 / runtime_seconds
    } else {
        0.0
    };
    let score = scaffold_link_f1.unwrap_or(0.0) * 0.6 + polished_exact_rate.unwrap_or(0.0) * 0.4;

    ScaffoldPolishTaskReport {
        task_dir: task_dir.display().to_string(),
        task_name,
        pairs_evaluated,
        expected_links: expected_links.len() as u64,
        observed_links: observed_links.len() as u64,
        correct_links,
        expected_scaffolds: truth.scaffolds.len() as u64,
        observed_scaffolds: observed_scaffolds.len() as u64,
        exact_scaffolds,
        polished_contigs_expected: truth.polished_contigs.len() as u64,
        polished_contigs_observed: polished_contigs.len() as u64,
        polished_contigs_exact: polished_exact_contigs,
        polished_bases_expected,
        polished_bases_correct,
        scaffold_link_precision,
        scaffold_link_recall,
        scaffold_link_f1,
        scaffold_exact_rate,
        polished_exact_rate,
        polished_base_accuracy,
        runtime_seconds,
        pairs_per_second,
        score,
    }
}

fn summarize_scaffold_polish_reports(
    tasks: &[ScaffoldPolishTaskReport],
) -> ScaffoldPolishAggregateReport {
    let pairs_evaluated = tasks.iter().map(|task| task.pairs_evaluated).sum();
    let expected_links = tasks.iter().map(|task| task.expected_links).sum();
    let observed_links = tasks.iter().map(|task| task.observed_links).sum();
    let correct_links = tasks.iter().map(|task| task.correct_links).sum();
    let expected_scaffolds = tasks.iter().map(|task| task.expected_scaffolds).sum();
    let observed_scaffolds = tasks.iter().map(|task| task.observed_scaffolds).sum();
    let exact_scaffolds = tasks.iter().map(|task| task.exact_scaffolds).sum();
    let polished_contigs_expected = tasks
        .iter()
        .map(|task| task.polished_contigs_expected)
        .sum();
    let polished_contigs_observed = tasks
        .iter()
        .map(|task| task.polished_contigs_observed)
        .sum();
    let polished_contigs_exact = tasks.iter().map(|task| task.polished_contigs_exact).sum();
    let polished_bases_expected = tasks.iter().map(|task| task.polished_bases_expected).sum();
    let polished_bases_correct = tasks.iter().map(|task| task.polished_bases_correct).sum();
    let runtime_seconds = tasks.iter().map(|task| task.runtime_seconds).sum();

    let scaffold_link_precision = ratio(correct_links, observed_links);
    let scaffold_link_recall = ratio(correct_links, expected_links);
    let scaffold_link_f1 = f1_score(
        scaffold_link_precision.unwrap_or(0.0),
        scaffold_link_recall.unwrap_or(0.0),
    );
    let scaffold_exact_rate = ratio(exact_scaffolds, expected_scaffolds);
    let polished_exact_rate = ratio(polished_contigs_exact, polished_contigs_expected);
    let polished_base_accuracy = ratio(polished_bases_correct, polished_bases_expected);
    let pairs_per_second = if runtime_seconds > 0.0 {
        pairs_evaluated as f64 / runtime_seconds
    } else {
        0.0
    };
    let score = scaffold_link_f1.unwrap_or(0.0) * 0.6 + polished_exact_rate.unwrap_or(0.0) * 0.4;

    ScaffoldPolishAggregateReport {
        task_count: tasks.len(),
        pairs_evaluated,
        expected_links,
        observed_links,
        correct_links,
        expected_scaffolds,
        observed_scaffolds,
        exact_scaffolds,
        polished_contigs_expected,
        polished_contigs_observed,
        polished_contigs_exact,
        polished_bases_expected,
        polished_bases_correct,
        scaffold_link_precision,
        scaffold_link_recall,
        scaffold_link_f1,
        scaffold_exact_rate,
        polished_exact_rate,
        polished_base_accuracy,
        runtime_seconds,
        pairs_per_second,
        score,
    }
}

pub fn run_error_correction_bench(
    task_path: &str,
    min_count: u32,
    min_trusted_count: u32,
    json: bool,
    output: Option<&str>,
) -> io::Result<()> {
    let report =
        evaluate_error_correction_bench(Path::new(task_path), min_count, min_trusted_count)?;

    if let Some(output_path) = output {
        let serialized = serde_json::to_string_pretty(&report)
            .map_err(|err| io::Error::other(format!("serialize report: {err}")))?;
        if let Some(parent) = Path::new(output_path)
            .parent()
            .filter(|path| !path.as_os_str().is_empty())
        {
            fs::create_dir_all(parent)?;
        }
        fs::write(output_path, serialized)?;
    }

    if json {
        let serialized = serde_json::to_string_pretty(&report)
            .map_err(|err| io::Error::other(format!("serialize report: {err}")))?;
        println!("{serialized}");
    } else {
        print_error_correction_report(&report);
    }

    Ok(())
}

pub fn prepare_error_correction_task(
    output_dir: &str,
    k: usize,
    trusted_roots: usize,
    weak_roots: usize,
    correctable_singletons_per_trusted: usize,
    protected_singletons_per_weak: usize,
    seed: u64,
) -> io::Result<()> {
    validate_error_correction_args(
        k,
        trusted_roots,
        weak_roots,
        correctable_singletons_per_trusted,
        protected_singletons_per_weak,
    )?;

    let output_dir = Path::new(output_dir);
    fs::create_dir_all(output_dir)?;

    let mut rng = StdRng::seed_from_u64(seed);
    let mut used_canonical = BTreeSet::new();
    let trusted_templates =
        collect_unique_canonical_kmers(k, trusted_roots, &mut used_canonical, &mut rng)?;
    let weak_templates =
        collect_unique_canonical_kmers(k, weak_roots, &mut used_canonical, &mut rng)?;
    let (task_file, truth_file) = build_error_correction_task_files(
        &trusted_templates,
        &weak_templates,
        correctable_singletons_per_trusted,
        protected_singletons_per_weak,
        &mut rng,
    )?;
    let metadata = ErrorCorrectionTaskMetadata {
        component: "error_correction".to_string(),
        name: output_dir
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("error_correction_task")
            .to_string(),
        description:
            "Synthetic phase-3 error-correction task with correctable singleton errors and protected rare variants."
                .to_string(),
        k,
        trusted_roots,
        weak_roots,
        correctable_singletons_per_trusted,
        protected_singletons_per_weak,
        kmers_generated: task_file.counts.len(),
        seed,
    };

    fs::write(
        output_dir.join("counts.json"),
        serde_json::to_string_pretty(&task_file)
            .map_err(|err| io::Error::other(format!("serialize counts: {err}")))?,
    )?;
    fs::write(
        output_dir.join("truth.json"),
        serde_json::to_string_pretty(&truth_file)
            .map_err(|err| io::Error::other(format!("serialize truth: {err}")))?,
    )?;
    fs::write(
        output_dir.join("task.json"),
        serde_json::to_string_pretty(&metadata)
            .map_err(|err| io::Error::other(format!("serialize metadata: {err}")))?,
    )?;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn prepare_error_correction_panel(
    output_dir: &str,
    tasks: usize,
    k: usize,
    trusted_roots: usize,
    weak_roots: usize,
    correctable_singletons_per_trusted: usize,
    protected_singletons_per_weak: usize,
    seed: u64,
    seed_step: u64,
) -> io::Result<()> {
    if tasks == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "tasks must be greater than zero",
        ));
    }

    let output_dir = Path::new(output_dir);
    fs::create_dir_all(output_dir)?;

    for task_index in 0..tasks {
        let task_seed = seed.saturating_add(seed_step.saturating_mul(task_index as u64));
        let task_dir = output_dir.join(format!("task_{task_index:03}"));
        prepare_error_correction_task(
            path_to_str(&task_dir)?,
            k,
            trusted_roots,
            weak_roots,
            correctable_singletons_per_trusted,
            protected_singletons_per_weak,
            task_seed,
        )?;
    }

    Ok(())
}

fn evaluate_error_correction_bench(
    task_path: &Path,
    min_count: u32,
    min_trusted_count: u32,
) -> io::Result<ErrorCorrectionBenchmarkReport> {
    let task_dirs = discover_error_correction_task_dirs(task_path)?;
    let tasks = task_dirs
        .par_iter()
        .map(|task_dir| {
            evaluate_single_error_correction_task(task_dir, min_count, min_trusted_count)
        })
        .collect::<io::Result<Vec<_>>>()?;
    let summary = summarize_error_correction_reports(&tasks);

    Ok(ErrorCorrectionBenchmarkReport {
        component: "error_correction",
        min_count,
        min_trusted_count,
        tasks,
        summary,
    })
}

fn discover_error_correction_task_dirs(task_path: &Path) -> io::Result<Vec<PathBuf>> {
    if is_error_correction_task_dir(task_path) {
        return Ok(vec![task_path.to_path_buf()]);
    }

    let mut discovered = Vec::new();
    for entry in fs::read_dir(task_path)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() && is_error_correction_task_dir(&path) {
            discovered.push(path);
        }
    }
    discovered.sort();

    if discovered.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "no error-correction tasks found under {}",
                task_path.display()
            ),
        ));
    }

    Ok(discovered)
}

fn is_error_correction_task_dir(task_dir: &Path) -> bool {
    resolve_error_correction_counts_path(task_dir).is_ok()
        && resolve_error_correction_truth_path(task_dir).is_ok()
}

fn evaluate_single_error_correction_task(
    task_dir: &Path,
    min_count: u32,
    min_trusted_count: u32,
) -> io::Result<ErrorCorrectionTaskReport> {
    let task_name = read_error_correction_task_metadata(task_dir)?
        .and_then(|task| {
            if task.name.trim().is_empty() {
                None
            } else {
                Some(task.name)
            }
        })
        .unwrap_or_else(|| {
            task_dir
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("error_correction_task")
                .to_string()
        });
    let metadata = read_error_correction_task_metadata(task_dir)?.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "missing error-correction metadata under {}",
                task_dir.display()
            ),
        )
    })?;
    let task_file = load_error_correction_task_file(task_dir)?;
    let truth = load_error_correction_truth_file(task_dir)?;
    let initial_counts = decode_error_correction_counts(&task_file.counts, metadata.k)?;

    let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
        k: metadata.k,
        min_count,
        error_correction_min_trusted_count: min_trusted_count,
        min_contig_len: metadata.k,
        ..Default::default()
    });

    let initial_singletons = initial_counts.values().filter(|&&count| count == 1).count() as u64;
    let start = Instant::now();
    let (corrected_counts, stats) = assembler.run_error_correction_case(initial_counts, metadata.k);
    let runtime_seconds = start.elapsed().as_secs_f64();

    build_error_correction_task_report(
        task_dir,
        task_name,
        initial_singletons,
        &truth,
        &corrected_counts,
        stats,
        runtime_seconds,
        metadata.k,
    )
}

#[allow(clippy::too_many_arguments)]
fn build_error_correction_task_report(
    task_dir: &Path,
    task_name: String,
    input_singletons: u64,
    truth: &ErrorCorrectionTruthFile,
    corrected_counts: &AHashMap<u64, u32>,
    stats: crate::pipeline::large_genome_assembler::LargeGenomeBenchErrorCorrectionStats,
    runtime_seconds: f64,
    k: usize,
) -> io::Result<ErrorCorrectionTaskReport> {
    let expected_removed = truth.expected_removed_singletons.len() as u64;
    let expected_preserved = truth.expected_preserved_singletons.len() as u64;
    let trusted_targets = truth.expected_trusted_totals.len() as u64;

    let expected_removed_encoded = encode_kmer_list(&truth.expected_removed_singletons, k)?;
    let expected_preserved_encoded = encode_kmer_list(&truth.expected_preserved_singletons, k)?;
    let expected_trusted_encoded = encode_trusted_totals(&truth.expected_trusted_totals, k)?;

    let correctly_removed = expected_removed_encoded
        .iter()
        .filter(|kmer| !corrected_counts.contains_key(kmer))
        .count() as u64;
    let preserved_retained = expected_preserved_encoded
        .iter()
        .filter(|kmer| corrected_counts.contains_key(kmer))
        .count() as u64;
    let observed_removed = expected_removed_encoded
        .iter()
        .chain(expected_preserved_encoded.iter())
        .filter(|kmer| !corrected_counts.contains_key(kmer))
        .count() as u64;
    let trusted_targets_exact = expected_trusted_encoded
        .iter()
        .filter(|(kmer, expected_count)| {
            corrected_counts.get(kmer).copied() == Some(**expected_count)
        })
        .count() as u64;

    let correction_precision = ratio(correctly_removed, observed_removed);
    let correction_recall = ratio(correctly_removed, expected_removed);
    let correction_f1 = f1_score(
        correction_precision.unwrap_or(0.0),
        correction_recall.unwrap_or(0.0),
    );
    let preserved_retention_rate = ratio(preserved_retained, expected_preserved);
    let trusted_target_exact_rate = ratio(trusted_targets_exact, trusted_targets);
    let kmers_per_second = if runtime_seconds > 0.0 {
        stats.remaining_kmers as f64 / runtime_seconds
    } else {
        0.0
    };
    let score = error_correction_score(
        correction_f1,
        preserved_retention_rate,
        trusted_target_exact_rate,
    );

    Ok(ErrorCorrectionTaskReport {
        task_dir: task_dir.display().to_string(),
        task_name,
        input_kmers: (stats.remaining_kmers as u64).saturating_add(stats.corrected_kmers),
        input_singletons,
        corrected_kmers: stats.corrected_kmers,
        remaining_kmers: stats.remaining_kmers as u64,
        expected_removed_singletons: expected_removed,
        observed_removed_singletons: observed_removed,
        correctly_removed_singletons: correctly_removed,
        expected_preserved_singletons: expected_preserved,
        preserved_singletons_retained: preserved_retained,
        trusted_targets,
        trusted_targets_exact,
        correction_precision,
        correction_recall,
        correction_f1,
        preserved_retention_rate,
        trusted_target_exact_rate,
        runtime_seconds,
        kmers_per_second,
        score,
    })
}

fn summarize_error_correction_reports(
    tasks: &[ErrorCorrectionTaskReport],
) -> ErrorCorrectionAggregateReport {
    let input_kmers = tasks.iter().map(|task| task.input_kmers).sum();
    let input_singletons = tasks.iter().map(|task| task.input_singletons).sum();
    let corrected_kmers = tasks.iter().map(|task| task.corrected_kmers).sum();
    let remaining_kmers = tasks.iter().map(|task| task.remaining_kmers).sum();
    let expected_removed_singletons = tasks
        .iter()
        .map(|task| task.expected_removed_singletons)
        .sum();
    let observed_removed_singletons = tasks
        .iter()
        .map(|task| task.observed_removed_singletons)
        .sum();
    let correctly_removed_singletons = tasks
        .iter()
        .map(|task| task.correctly_removed_singletons)
        .sum();
    let expected_preserved_singletons = tasks
        .iter()
        .map(|task| task.expected_preserved_singletons)
        .sum();
    let preserved_singletons_retained = tasks
        .iter()
        .map(|task| task.preserved_singletons_retained)
        .sum();
    let trusted_targets = tasks.iter().map(|task| task.trusted_targets).sum();
    let trusted_targets_exact = tasks.iter().map(|task| task.trusted_targets_exact).sum();
    let runtime_seconds = tasks.iter().map(|task| task.runtime_seconds).sum();

    let correction_precision = ratio(correctly_removed_singletons, observed_removed_singletons);
    let correction_recall = ratio(correctly_removed_singletons, expected_removed_singletons);
    let correction_f1 = f1_score(
        correction_precision.unwrap_or(0.0),
        correction_recall.unwrap_or(0.0),
    );
    let preserved_retention_rate =
        ratio(preserved_singletons_retained, expected_preserved_singletons);
    let trusted_target_exact_rate = ratio(trusted_targets_exact, trusted_targets);
    let kmers_per_second = if runtime_seconds > 0.0 {
        remaining_kmers as f64 / runtime_seconds
    } else {
        0.0
    };
    let score = error_correction_score(
        correction_f1,
        preserved_retention_rate,
        trusted_target_exact_rate,
    );

    ErrorCorrectionAggregateReport {
        task_count: tasks.len(),
        input_kmers,
        input_singletons,
        corrected_kmers,
        remaining_kmers,
        expected_removed_singletons,
        observed_removed_singletons,
        correctly_removed_singletons,
        expected_preserved_singletons,
        preserved_singletons_retained,
        trusted_targets,
        trusted_targets_exact,
        correction_precision,
        correction_recall,
        correction_f1,
        preserved_retention_rate,
        trusted_target_exact_rate,
        runtime_seconds,
        kmers_per_second,
        score,
    }
}

pub fn prepare_contig_extraction_task(
    output_dir: &str,
    k: usize,
    component_count: usize,
    primary_reads_per_component: usize,
    alternate_reads_per_component: usize,
    seed: u64,
) -> io::Result<()> {
    prepare_contig_extraction_task_with_profile(
        output_dir,
        k,
        component_count,
        primary_reads_per_component,
        alternate_reads_per_component,
        seed,
        "branching",
    )
}

#[allow(clippy::too_many_arguments)]
pub fn prepare_contig_extraction_task_with_profile(
    output_dir: &str,
    k: usize,
    component_count: usize,
    primary_reads_per_component: usize,
    alternate_reads_per_component: usize,
    seed: u64,
    scenario_profile: &str,
) -> io::Result<()> {
    validate_contig_extraction_args(
        k,
        component_count,
        primary_reads_per_component,
        alternate_reads_per_component,
    )?;

    let profile = ContigExtractionScenarioProfile::parse(scenario_profile)?;
    let (counts, graph, truth, metadata) = build_contig_extraction_task_files(
        k,
        component_count,
        primary_reads_per_component,
        alternate_reads_per_component,
        seed,
        profile,
    )?;

    let output_dir = Path::new(output_dir);
    fs::create_dir_all(output_dir)?;
    fs::write(
        output_dir.join("counts.json"),
        serde_json::to_string_pretty(&counts)
            .map_err(|err| io::Error::other(format!("serialize counts: {err}")))?,
    )?;
    fs::write(
        output_dir.join("graph.json"),
        serde_json::to_string_pretty(&graph)
            .map_err(|err| io::Error::other(format!("serialize graph: {err}")))?,
    )?;
    fs::write(
        output_dir.join("truth.json"),
        serde_json::to_string_pretty(&truth)
            .map_err(|err| io::Error::other(format!("serialize truth: {err}")))?,
    )?;
    fs::write(
        output_dir.join("task.json"),
        serde_json::to_string_pretty(&metadata)
            .map_err(|err| io::Error::other(format!("serialize metadata: {err}")))?,
    )?;

    Ok(())
}

pub fn prepare_contig_extraction_panel(
    output_dir: &str,
    tasks: usize,
    k: usize,
    component_count: usize,
    primary_reads_per_component: usize,
    alternate_reads_per_component: usize,
    seed: u64,
    seed_step: u64,
) -> io::Result<()> {
    prepare_contig_extraction_panel_with_profile(
        output_dir,
        tasks,
        k,
        component_count,
        primary_reads_per_component,
        alternate_reads_per_component,
        seed,
        seed_step,
        "mixed",
    )
}

#[allow(clippy::too_many_arguments)]
pub fn prepare_contig_extraction_panel_with_profile(
    output_dir: &str,
    tasks: usize,
    k: usize,
    component_count: usize,
    primary_reads_per_component: usize,
    alternate_reads_per_component: usize,
    seed: u64,
    seed_step: u64,
    scenario_profile: &str,
) -> io::Result<()> {
    if tasks == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "tasks must be greater than zero",
        ));
    }

    let output_dir = Path::new(output_dir);
    fs::create_dir_all(output_dir)?;
    let profiles = match scenario_profile {
        "mixed" => vec![
            ContigExtractionScenarioProfile::Branching,
            ContigExtractionScenarioProfile::RepeatCompletion,
            ContigExtractionScenarioProfile::RepeatFallback,
        ],
        other => vec![ContigExtractionScenarioProfile::parse(other)?],
    };

    for task_index in 0..tasks {
        let task_seed = seed.saturating_add(seed_step.saturating_mul(task_index as u64));
        let profile = profiles[task_index % profiles.len()];
        let task_dir = output_dir.join(format!("task_{task_index:03}"));
        prepare_contig_extraction_task_with_profile(
            path_to_str(&task_dir)?,
            k,
            component_count,
            primary_reads_per_component,
            alternate_reads_per_component,
            task_seed,
            profile.as_str(),
        )?;
    }

    Ok(())
}

fn evaluate_contig_extraction_bench(
    task_path: &Path,
    prefer_high_count_seeds: bool,
    prefer_non_repeat_seeds: bool,
    enable_repeat_seed_completion: bool,
    suppress_redundant_contigs: bool,
) -> io::Result<ContigExtractionBenchmarkReport> {
    let task_dirs = discover_contig_extraction_task_dirs(task_path)?;
    let tasks = task_dirs
        .par_iter()
        .map(|task_dir| {
            evaluate_single_contig_extraction_task(
                task_dir,
                prefer_high_count_seeds,
                prefer_non_repeat_seeds,
                enable_repeat_seed_completion,
                suppress_redundant_contigs,
            )
        })
        .collect::<io::Result<Vec<_>>>()?;
    let summary = summarize_contig_extraction_reports(&tasks);

    Ok(ContigExtractionBenchmarkReport {
        component: "contig_extraction",
        prefer_high_count_seeds,
        prefer_non_repeat_seeds,
        enable_repeat_seed_completion,
        suppress_redundant_contigs,
        tasks,
        summary,
    })
}

fn discover_contig_extraction_task_dirs(task_path: &Path) -> io::Result<Vec<PathBuf>> {
    if is_contig_extraction_task_dir(task_path) {
        return Ok(vec![task_path.to_path_buf()]);
    }

    let mut discovered = Vec::new();
    for entry in fs::read_dir(task_path)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() && is_contig_extraction_task_dir(&path) {
            discovered.push(path);
        }
    }
    discovered.sort();

    if discovered.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "no contig-extraction tasks found under {}",
                task_path.display()
            ),
        ));
    }

    Ok(discovered)
}

fn is_contig_extraction_task_dir(task_dir: &Path) -> bool {
    resolve_contig_extraction_counts_path(task_dir).is_ok()
        && resolve_contig_extraction_graph_path(task_dir).is_ok()
        && resolve_contig_extraction_truth_path(task_dir).is_ok()
}

fn evaluate_single_contig_extraction_task(
    task_dir: &Path,
    prefer_high_count_seeds: bool,
    prefer_non_repeat_seeds: bool,
    enable_repeat_seed_completion: bool,
    suppress_redundant_contigs: bool,
) -> io::Result<ContigExtractionTaskReport> {
    let task_name = read_contig_extraction_task_metadata(task_dir)?
        .as_ref()
        .and_then(|task| {
            if task.name.trim().is_empty() {
                None
            } else {
                Some(task.name.clone())
            }
        })
        .unwrap_or_else(|| {
            task_dir
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("contig_extraction_task")
                .to_string()
        });
    let metadata = read_contig_extraction_task_metadata(task_dir)?.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "missing contig-extraction metadata under {}",
                task_dir.display()
            ),
        )
    })?;
    let counts_file = load_contig_extraction_counts_file(task_dir)?;
    let graph_file = load_contig_extraction_graph_file(task_dir)?;
    let truth_file = load_contig_extraction_truth_file(task_dir)?;
    let kmer_counts = decode_contig_extraction_counts(&counts_file, metadata.k)?;
    let valid_kmers = decode_contig_extraction_valid_kmers(&graph_file.valid_kmers)?;
    let branch_support = decode_contig_extraction_branch_support(&graph_file.branch_support)?;

    let assembler = LargeGenomeAssembler::new(LargeGenomeConfig {
        k: metadata.k,
        min_count: 1,
        min_contig_len: metadata.k,
        prefer_high_count_seeds,
        prefer_non_repeat_seeds,
        enable_repeat_seed_completion,
        suppress_redundant_contigs,
        ..Default::default()
    });

    let graph_nodes = valid_kmers.len() as u64;
    let start = Instant::now();
    let contigs = assembler.run_contig_extraction_case(
        kmer_counts,
        &valid_kmers,
        &branch_support,
        metadata.k,
    );
    let runtime_seconds = start.elapsed().as_secs_f64();

    build_contig_extraction_task_report(
        task_dir,
        task_name,
        metadata.profile.as_str().to_string(),
        graph_nodes,
        &truth_file,
        &contigs,
        runtime_seconds,
        metadata.k,
    )
}

#[allow(clippy::too_many_arguments)]
fn build_contig_extraction_task_report(
    task_dir: &Path,
    task_name: String,
    profile: String,
    graph_nodes: u64,
    truth: &ContigExtractionTruthFile,
    observed_contigs: &[String],
    runtime_seconds: f64,
    k: usize,
) -> io::Result<ContigExtractionTaskReport> {
    let expected_contigs = truth.expected_contigs.len() as u64;
    let observed_contig_count = observed_contigs.len() as u64;
    let expected_canonical = canonical_contig_set(&truth.expected_contigs);
    let observed_canonical = canonical_contig_set(observed_contigs);
    let exact_contigs = expected_canonical.intersection(&observed_canonical).count() as u64;

    let expected_kmers = canonical_kmer_set_from_sequences(&truth.expected_contigs, k)?;
    let observed_kmers = canonical_kmer_set_from_sequences(observed_contigs, k)?;
    let truth_kmers_correct = expected_kmers.intersection(&observed_kmers).count() as u64;
    let truth_kmers_expected = expected_kmers.len() as u64;
    let truth_kmers_observed = observed_kmers.len() as u64;

    let exact_contig_rate = ratio(exact_contigs, expected_contigs);
    let contig_count_agreement = ratio_f64(
        expected_contigs.min(observed_contig_count) as f64,
        expected_contigs.max(observed_contig_count) as f64,
    );
    let truth_kmer_precision = ratio(truth_kmers_correct, truth_kmers_observed);
    let truth_kmer_recall = ratio(truth_kmers_correct, truth_kmers_expected);
    let truth_kmer_f1 = f1_score(
        truth_kmer_precision.unwrap_or(0.0),
        truth_kmer_recall.unwrap_or(0.0),
    );
    let graph_nodes_per_second = if runtime_seconds > 0.0 {
        graph_nodes as f64 / runtime_seconds
    } else {
        0.0
    };
    let score = contig_extraction_score(truth_kmer_f1, exact_contig_rate, contig_count_agreement);

    Ok(ContigExtractionTaskReport {
        task_dir: task_dir.display().to_string(),
        task_name,
        profile,
        graph_nodes,
        expected_contigs,
        observed_contigs: observed_contig_count,
        exact_contigs,
        truth_kmers_expected,
        truth_kmers_observed,
        truth_kmers_correct,
        exact_contig_rate,
        contig_count_agreement,
        truth_kmer_precision,
        truth_kmer_recall,
        truth_kmer_f1,
        runtime_seconds,
        graph_nodes_per_second,
        score,
    })
}

fn summarize_contig_extraction_reports(
    tasks: &[ContigExtractionTaskReport],
) -> ContigExtractionAggregateReport {
    let graph_nodes = tasks.iter().map(|task| task.graph_nodes).sum();
    let expected_contigs = tasks.iter().map(|task| task.expected_contigs).sum();
    let observed_contigs = tasks.iter().map(|task| task.observed_contigs).sum();
    let exact_contigs = tasks.iter().map(|task| task.exact_contigs).sum();
    let truth_kmers_expected = tasks.iter().map(|task| task.truth_kmers_expected).sum();
    let truth_kmers_observed = tasks.iter().map(|task| task.truth_kmers_observed).sum();
    let truth_kmers_correct = tasks.iter().map(|task| task.truth_kmers_correct).sum();
    let runtime_seconds = tasks.iter().map(|task| task.runtime_seconds).sum();

    let exact_contig_rate = ratio(exact_contigs, expected_contigs);
    let contig_count_agreement = ratio_f64(
        expected_contigs.min(observed_contigs) as f64,
        expected_contigs.max(observed_contigs) as f64,
    );
    let truth_kmer_precision = ratio(truth_kmers_correct, truth_kmers_observed);
    let truth_kmer_recall = ratio(truth_kmers_correct, truth_kmers_expected);
    let truth_kmer_f1 = f1_score(
        truth_kmer_precision.unwrap_or(0.0),
        truth_kmer_recall.unwrap_or(0.0),
    );
    let graph_nodes_per_second = if runtime_seconds > 0.0 {
        graph_nodes as f64 / runtime_seconds
    } else {
        0.0
    };
    let score = contig_extraction_score(truth_kmer_f1, exact_contig_rate, contig_count_agreement);

    ContigExtractionAggregateReport {
        task_count: tasks.len(),
        graph_nodes,
        expected_contigs,
        observed_contigs,
        exact_contigs,
        truth_kmers_expected,
        truth_kmers_observed,
        truth_kmers_correct,
        exact_contig_rate,
        contig_count_agreement,
        truth_kmer_precision,
        truth_kmer_recall,
        truth_kmer_f1,
        runtime_seconds,
        graph_nodes_per_second,
        score,
    }
}

fn evaluate_branch_resolution_bench(
    task_path: &Path,
    branch_support_min_win: u32,
    branch_support_min_margin: u32,
    prefer_non_repeat_branches: bool,
) -> io::Result<BranchResolutionBenchmarkReport> {
    let task_dirs = discover_branch_resolution_task_dirs(task_path)?;
    let config = LargeGenomeConfig {
        branch_support_min_win,
        branch_support_min_margin,
        prefer_non_repeat_branches,
        ..Default::default()
    };
    let tasks = task_dirs
        .par_iter()
        .map(|task_dir| evaluate_single_branch_resolution_task(task_dir, &config))
        .collect::<io::Result<Vec<_>>>()?;

    let summary = summarize_branch_resolution_reports(&tasks);
    Ok(BranchResolutionBenchmarkReport {
        component: "branch_resolution",
        branch_support_min_win,
        branch_support_min_margin,
        prefer_non_repeat_branches,
        tasks,
        summary,
    })
}

fn discover_branch_resolution_task_dirs(task_path: &Path) -> io::Result<Vec<PathBuf>> {
    if is_branch_resolution_task_dir(task_path) {
        return Ok(vec![task_path.to_path_buf()]);
    }

    let mut discovered = Vec::new();
    for entry in fs::read_dir(task_path)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() && is_branch_resolution_task_dir(&path) {
            discovered.push(path);
        }
    }
    discovered.sort();

    if discovered.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "no branch-resolution tasks found under {}",
                task_path.display()
            ),
        ));
    }

    Ok(discovered)
}

fn is_branch_resolution_task_dir(task_dir: &Path) -> bool {
    resolve_branch_cases_path(task_dir).is_ok()
}

fn evaluate_single_branch_resolution_task(
    task_dir: &Path,
    config: &LargeGenomeConfig,
) -> io::Result<BranchResolutionTaskReport> {
    let task_name = read_branch_task_metadata(task_dir)?
        .and_then(|task| {
            if task.name.trim().is_empty() {
                None
            } else {
                Some(task.name)
            }
        })
        .unwrap_or_else(|| {
            task_dir
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("branch_resolution_task")
                .to_string()
        });
    let task_file = load_branch_resolution_task_file(task_dir)?;
    let assembler = LargeGenomeAssembler::new(config.clone());

    let start = Instant::now();
    let mut counters = BranchResolutionCounters::default();

    for case in &task_file.cases {
        let expected_next = encode_branch_kmer(&case.expected_next_kmer)?;
        let choice = assembler.run_branch_resolution_case(&LargeGenomeBenchBranchResolutionCase {
            current_count: case.current_count,
            candidates: case
                .candidates
                .iter()
                .map(|candidate| {
                    Ok(LargeGenomeBenchBranchCandidate {
                        base_idx: candidate.base_idx,
                        next: encode_branch_kmer(&candidate.next_kmer)?,
                        count: candidate.count,
                        is_repeat: candidate.is_repeat,
                        read_support: candidate.read_support,
                    })
                })
                .collect::<io::Result<Vec<_>>>()?,
        });

        let matched = choice
            .map(|selected| {
                selected.base_idx == case.expected_base_idx && selected.next == expected_next
            })
            .unwrap_or(false);
        counters.observe(&case.scenario, matched, choice.is_some());
    }

    let runtime_seconds = start.elapsed().as_secs_f64();
    Ok(build_branch_resolution_task_report(
        task_dir,
        task_name,
        counters,
        runtime_seconds,
    ))
}

fn build_branch_resolution_task_report(
    task_dir: &Path,
    task_name: String,
    counters: BranchResolutionCounters,
    runtime_seconds: f64,
) -> BranchResolutionTaskReport {
    let cases_per_second = if runtime_seconds > 0.0 {
        counters.cases_evaluated as f64 / runtime_seconds
    } else {
        0.0
    };
    let choice_rate = ratio(counters.cases_with_choice, counters.cases_evaluated);
    let exact_rate = ratio(counters.exact_matches, counters.cases_evaluated);

    BranchResolutionTaskReport {
        task_dir: task_dir.display().to_string(),
        task_name,
        cases_evaluated: counters.cases_evaluated,
        cases_with_choice: counters.cases_with_choice,
        exact_matches: counters.exact_matches,
        choice_rate,
        exact_rate,
        runtime_seconds,
        cases_per_second,
        score: exact_rate.unwrap_or(0.0),
        scenario_metrics: branch_resolution_scenario_metrics(&counters),
    }
}

fn summarize_branch_resolution_reports(
    tasks: &[BranchResolutionTaskReport],
) -> BranchResolutionAggregateReport {
    let mut counters = BranchResolutionCounters::default();
    let mut runtime_seconds = 0.0;

    for task in tasks {
        counters.merge(&BranchResolutionCounters {
            cases_evaluated: task.cases_evaluated,
            cases_with_choice: task.cases_with_choice,
            exact_matches: task.exact_matches,
            scenario_counts: task
                .scenario_metrics
                .iter()
                .map(|(scenario, metrics)| (scenario.clone(), metrics.cases_evaluated))
                .collect(),
            scenario_correct: task
                .scenario_metrics
                .iter()
                .map(|(scenario, metrics)| (scenario.clone(), metrics.exact_matches))
                .collect(),
        });
        runtime_seconds += task.runtime_seconds;
    }

    let cases_per_second = if runtime_seconds > 0.0 {
        counters.cases_evaluated as f64 / runtime_seconds
    } else {
        0.0
    };
    let choice_rate = ratio(counters.cases_with_choice, counters.cases_evaluated);
    let exact_rate = ratio(counters.exact_matches, counters.cases_evaluated);

    BranchResolutionAggregateReport {
        task_count: tasks.len(),
        cases_evaluated: counters.cases_evaluated,
        cases_with_choice: counters.cases_with_choice,
        exact_matches: counters.exact_matches,
        choice_rate,
        exact_rate,
        runtime_seconds,
        cases_per_second,
        score: exact_rate.unwrap_or(0.0),
        scenario_metrics: branch_resolution_scenario_metrics(&counters),
    }
}

fn branch_resolution_scenario_metrics(
    counters: &BranchResolutionCounters,
) -> BTreeMap<String, BranchResolutionScenarioMetric> {
    let mut metrics = BTreeMap::new();
    for (scenario, cases_evaluated) in &counters.scenario_counts {
        let exact_matches = counters
            .scenario_correct
            .get(scenario)
            .copied()
            .unwrap_or(0);
        metrics.insert(
            scenario.clone(),
            BranchResolutionScenarioMetric {
                cases_evaluated: *cases_evaluated,
                exact_matches,
                exact_rate: ratio(exact_matches, *cases_evaluated),
            },
        );
    }
    metrics
}

fn load_branch_resolution_task_file(task_dir: &Path) -> io::Result<BranchResolutionTaskFile> {
    let cases_path = resolve_branch_cases_path(task_dir)?;
    serde_json::from_reader(BufReader::new(File::open(&cases_path)?)).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("parse {}: {}", cases_path.display(), err),
        )
    })
}

fn read_branch_task_metadata(task_dir: &Path) -> io::Result<Option<BranchResolutionTaskMetadata>> {
    let metadata_path = task_dir.join("task.json");
    if !metadata_path.exists() {
        return Ok(None);
    }

    let metadata =
        serde_json::from_reader(BufReader::new(File::open(&metadata_path)?)).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("parse {}: {}", metadata_path.display(), err),
            )
        })?;
    Ok(Some(metadata))
}

fn resolve_branch_cases_path(task_dir: &Path) -> io::Result<PathBuf> {
    resolve_required_path(task_dir, &["cases.json"])
}

fn encode_branch_kmer(sequence: &str) -> io::Result<u64> {
    KmerU64::from_str(sequence)
        .map(|encoded| encoded.encoded)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid branch k-mer '{}'", sequence),
            )
        })
}

fn build_branch_resolution_cases(
    cases: usize,
    seed: u64,
    scenario_profile: BranchScenarioProfile,
) -> Vec<BranchResolutionCaseSpec> {
    let mut rng = StdRng::seed_from_u64(seed);
    (0..cases)
        .map(
            |case_idx| match branch_scenario_for_case(case_idx, scenario_profile) {
                "support_dominates" => branch_support_dominates_case(case_idx, &mut rng),
                "support_floor_gate" => branch_support_floor_gate_case(case_idx, &mut rng),
                "support_margin_gate" => branch_support_margin_gate_case(case_idx, &mut rng),
                "coverage_closeness" => branch_coverage_closeness_case(case_idx, &mut rng),
                "non_repeat_preferred" => branch_non_repeat_case(case_idx, &mut rng),
                _ => branch_repeat_tiebreak_case(case_idx, &mut rng),
            },
        )
        .collect()
}

fn branch_scenario_for_case(
    case_idx: usize,
    scenario_profile: BranchScenarioProfile,
) -> &'static str {
    match scenario_profile {
        BranchScenarioProfile::Balanced => match case_idx % 6 {
            0 => "support_dominates",
            1 => "support_floor_gate",
            2 => "support_margin_gate",
            3 => "coverage_closeness",
            4 => "non_repeat_preferred",
            _ => "repeat_tiebreak",
        },
        BranchScenarioProfile::SupportStress => match case_idx % 8 {
            0 => "support_floor_gate",
            1 => "support_margin_gate",
            2 => "coverage_closeness",
            3 => "support_margin_gate",
            4 => "non_repeat_preferred",
            5 => "support_floor_gate",
            6 => "coverage_closeness",
            _ => "repeat_tiebreak",
        },
    }
}

fn branch_support_dominates_case(case_idx: usize, rng: &mut StdRng) -> BranchResolutionCaseSpec {
    let current_count: u32 = rng.gen_range(18..=42);
    let mut kmers = unique_kmers(DEFAULT_BRANCH_KMER_LEN, 3, rng);
    let winning_support: u32 = rng.gen_range(2..=3);
    let losing_support = winning_support.saturating_sub(2);

    BranchResolutionCaseSpec {
        scenario: "support_dominates".to_string(),
        current_count,
        expected_base_idx: case_idx % 4,
        expected_next_kmer: kmers[0].clone(),
        candidates: vec![
            BranchResolutionCandidateSpec {
                base_idx: case_idx % 4,
                next_kmer: kmers.remove(0),
                count: current_count.saturating_add(rng.gen_range(8..=16)),
                is_repeat: false,
                read_support: winning_support,
            },
            BranchResolutionCandidateSpec {
                base_idx: (case_idx + 1) % 4,
                next_kmer: kmers.remove(0),
                count: current_count.saturating_sub(rng.gen_range(0..=1)),
                is_repeat: false,
                read_support: losing_support,
            },
            BranchResolutionCandidateSpec {
                base_idx: (case_idx + 2) % 4,
                next_kmer: kmers.remove(0),
                count: current_count.saturating_add(rng.gen_range(2..=6)),
                is_repeat: true,
                read_support: rng.gen_range(0..=1),
            },
        ],
    }
}

fn branch_support_floor_gate_case(case_idx: usize, rng: &mut StdRng) -> BranchResolutionCaseSpec {
    let current_count: u32 = rng.gen_range(18..=36);
    let mut kmers = unique_kmers(DEFAULT_BRANCH_KMER_LEN, 3, rng);
    let expected_count = current_count.saturating_sub(rng.gen_range(0..=1));

    BranchResolutionCaseSpec {
        scenario: "support_floor_gate".to_string(),
        current_count,
        expected_base_idx: case_idx % 4,
        expected_next_kmer: kmers[0].clone(),
        candidates: vec![
            BranchResolutionCandidateSpec {
                base_idx: case_idx % 4,
                next_kmer: kmers.remove(0),
                count: expected_count,
                is_repeat: false,
                read_support: 0,
            },
            BranchResolutionCandidateSpec {
                base_idx: (case_idx + 1) % 4,
                next_kmer: kmers.remove(0),
                count: current_count.saturating_add(rng.gen_range(8..=16)),
                is_repeat: false,
                read_support: 1,
            },
            BranchResolutionCandidateSpec {
                base_idx: (case_idx + 2) % 4,
                next_kmer: kmers.remove(0),
                count: current_count.saturating_add(rng.gen_range(2..=6)),
                is_repeat: true,
                read_support: 0,
            },
        ],
    }
}

fn branch_support_margin_gate_case(case_idx: usize, rng: &mut StdRng) -> BranchResolutionCaseSpec {
    let current_count: u32 = rng.gen_range(20..=40);
    let mut kmers = unique_kmers(DEFAULT_BRANCH_KMER_LEN, 3, rng);

    BranchResolutionCaseSpec {
        scenario: "support_margin_gate".to_string(),
        current_count,
        expected_base_idx: case_idx % 4,
        expected_next_kmer: kmers[0].clone(),
        candidates: vec![
            BranchResolutionCandidateSpec {
                base_idx: case_idx % 4,
                next_kmer: kmers.remove(0),
                count: current_count.saturating_sub(rng.gen_range(0..=1)),
                is_repeat: false,
                read_support: 2,
            },
            BranchResolutionCandidateSpec {
                base_idx: (case_idx + 1) % 4,
                next_kmer: kmers.remove(0),
                count: current_count.saturating_add(rng.gen_range(10..=16)),
                is_repeat: false,
                read_support: 3,
            },
            BranchResolutionCandidateSpec {
                base_idx: (case_idx + 2) % 4,
                next_kmer: kmers.remove(0),
                count: current_count.saturating_add(rng.gen_range(2..=6)),
                is_repeat: true,
                read_support: 1,
            },
        ],
    }
}

fn branch_coverage_closeness_case(case_idx: usize, rng: &mut StdRng) -> BranchResolutionCaseSpec {
    let current_count: u32 = rng.gen_range(20..=48);
    let mut kmers = unique_kmers(DEFAULT_BRANCH_KMER_LEN, 3, rng);
    let winning_delta = rng.gen_range(0..=2);
    let losing_delta = rng.gen_range(9..=16);

    BranchResolutionCaseSpec {
        scenario: "coverage_closeness".to_string(),
        current_count,
        expected_base_idx: case_idx % 4,
        expected_next_kmer: kmers[0].clone(),
        candidates: vec![
            BranchResolutionCandidateSpec {
                base_idx: case_idx % 4,
                next_kmer: kmers.remove(0),
                count: current_count.saturating_sub(winning_delta),
                is_repeat: false,
                read_support: rng.gen_range(0..=1),
            },
            BranchResolutionCandidateSpec {
                base_idx: (case_idx + 1) % 4,
                next_kmer: kmers.remove(0),
                count: current_count.saturating_add(losing_delta),
                is_repeat: false,
                read_support: rng.gen_range(0..=1),
            },
            BranchResolutionCandidateSpec {
                base_idx: (case_idx + 2) % 4,
                next_kmer: kmers.remove(0),
                count: current_count.saturating_sub(rng.gen_range(7..=12)),
                is_repeat: true,
                read_support: 0,
            },
        ],
    }
}

fn branch_non_repeat_case(case_idx: usize, rng: &mut StdRng) -> BranchResolutionCaseSpec {
    let current_count: u32 = rng.gen_range(16..=36);
    let mut kmers = unique_kmers(DEFAULT_BRANCH_KMER_LEN, 3, rng);
    let winning_count = current_count.saturating_sub(rng.gen_range(0..=2));

    BranchResolutionCaseSpec {
        scenario: "non_repeat_preferred".to_string(),
        current_count,
        expected_base_idx: case_idx % 4,
        expected_next_kmer: kmers[0].clone(),
        candidates: vec![
            BranchResolutionCandidateSpec {
                base_idx: case_idx % 4,
                next_kmer: kmers.remove(0),
                count: winning_count,
                is_repeat: false,
                read_support: rng.gen_range(0..=1),
            },
            BranchResolutionCandidateSpec {
                base_idx: (case_idx + 1) % 4,
                next_kmer: kmers.remove(0),
                count: winning_count.saturating_add(rng.gen_range(0..=2)),
                is_repeat: true,
                read_support: rng.gen_range(0..=1),
            },
            BranchResolutionCandidateSpec {
                base_idx: (case_idx + 2) % 4,
                next_kmer: kmers.remove(0),
                count: winning_count.saturating_add(rng.gen_range(6..=10)),
                is_repeat: true,
                read_support: 0,
            },
        ],
    }
}

fn branch_repeat_tiebreak_case(case_idx: usize, rng: &mut StdRng) -> BranchResolutionCaseSpec {
    let current_count: u32 = rng.gen_range(12..=28);
    let mut kmers = unique_kmers(DEFAULT_BRANCH_KMER_LEN, 3, rng);
    kmers.sort();
    let winning = kmers[0].clone();

    BranchResolutionCaseSpec {
        scenario: "repeat_tiebreak".to_string(),
        current_count,
        expected_base_idx: case_idx % 4,
        expected_next_kmer: winning,
        candidates: vec![
            BranchResolutionCandidateSpec {
                base_idx: case_idx % 4,
                next_kmer: kmers.remove(0),
                count: current_count,
                is_repeat: true,
                read_support: 2,
            },
            BranchResolutionCandidateSpec {
                base_idx: (case_idx + 1) % 4,
                next_kmer: kmers.remove(0),
                count: current_count,
                is_repeat: true,
                read_support: 2,
            },
            BranchResolutionCandidateSpec {
                base_idx: (case_idx + 2) % 4,
                next_kmer: kmers.remove(0),
                count: current_count.saturating_sub(2),
                is_repeat: true,
                read_support: 1,
            },
        ],
    }
}

fn unique_kmers(k: usize, count: usize, rng: &mut StdRng) -> Vec<String> {
    let mut kmers = Vec::with_capacity(count);
    while kmers.len() < count {
        let candidate = String::from_utf8(random_dna(k, rng)).unwrap_or_default();
        if !kmers.contains(&candidate) {
            kmers.push(candidate);
        }
    }
    kmers
}

fn build_error_correction_task_files(
    trusted_templates: &[String],
    weak_templates: &[String],
    correctable_singletons_per_trusted: usize,
    protected_singletons_per_weak: usize,
    rng: &mut StdRng,
) -> io::Result<(ErrorCorrectionTaskFile, ErrorCorrectionTruthFile)> {
    let mut used_canonical = BTreeSet::new();
    let mut counts = Vec::new();
    let mut expected_removed_singletons = Vec::new();
    let mut expected_preserved_singletons = Vec::new();
    let mut expected_trusted_totals = BTreeMap::new();

    let trusted_bytes = trusted_templates
        .iter()
        .map(|template| template.as_bytes().to_vec())
        .collect::<Vec<_>>();

    for (idx, template) in trusted_templates.iter().enumerate() {
        let encoded = encode_canonical_kmer(template)?;
        used_canonical.insert(encoded);

        let root_count = match idx % 3 {
            0 => 4,
            1 => 6,
            _ => 9,
        };
        counts.push(ErrorCorrectionCountSpec {
            kmer: template.clone(),
            count: root_count,
        });

        let mut expected_total = root_count;
        for _ in 0..correctable_singletons_per_trusted {
            let singleton = generate_unique_single_edit_kmer(
                template.as_bytes(),
                &mut used_canonical,
                &[],
                rng,
            )?;
            counts.push(ErrorCorrectionCountSpec {
                kmer: singleton.clone(),
                count: 1,
            });
            expected_removed_singletons.push(singleton);
            expected_total += 1;
        }

        expected_trusted_totals.insert(template.clone(), expected_total);
    }

    for (idx, template) in weak_templates.iter().enumerate() {
        let encoded = encode_canonical_kmer(template)?;
        used_canonical.insert(encoded);

        counts.push(ErrorCorrectionCountSpec {
            kmer: template.clone(),
            count: if idx % 2 == 0 { 2 } else { 3 },
        });

        for _ in 0..protected_singletons_per_weak {
            let singleton = generate_unique_single_edit_kmer(
                template.as_bytes(),
                &mut used_canonical,
                &trusted_bytes,
                rng,
            )?;
            counts.push(ErrorCorrectionCountSpec {
                kmer: singleton.clone(),
                count: 1,
            });
            expected_preserved_singletons.push(singleton);
        }
    }

    counts.sort_by(|left, right| left.kmer.cmp(&right.kmer));
    expected_removed_singletons.sort();
    expected_preserved_singletons.sort();

    Ok((
        ErrorCorrectionTaskFile { counts },
        ErrorCorrectionTruthFile {
            expected_removed_singletons,
            expected_preserved_singletons,
            expected_trusted_totals,
        },
    ))
}

fn collect_unique_canonical_kmers(
    k: usize,
    count: usize,
    used_canonical: &mut BTreeSet<u64>,
    rng: &mut StdRng,
) -> io::Result<Vec<String>> {
    let mut kmers = Vec::with_capacity(count);
    while kmers.len() < count {
        let candidate = String::from_utf8(random_dna(k, rng)).map_err(invalid_utf8_error)?;
        let encoded = encode_canonical_kmer(&candidate)?;
        if used_canonical.insert(encoded) {
            kmers.push(candidate);
        }
    }
    Ok(kmers)
}

fn generate_unique_single_edit_kmer(
    template: &[u8],
    used_canonical: &mut BTreeSet<u64>,
    disallow_hamming_one_to: &[Vec<u8>],
    rng: &mut StdRng,
) -> io::Result<String> {
    for attempt in 0..4096usize {
        let mut mutated = template.to_vec();
        let pos = (attempt * 7 + rng.gen_range(0..template.len())) % template.len();
        mutated[pos] = random_base_excluding(mutated[pos], rng);

        if disallow_hamming_one_to
            .iter()
            .any(|other| hamming_distance_is_one(&mutated, other))
        {
            continue;
        }

        let candidate = String::from_utf8(mutated).map_err(invalid_utf8_error)?;
        let encoded = encode_canonical_kmer(&candidate)?;
        if used_canonical.insert(encoded) {
            return Ok(candidate);
        }
    }

    Err(io::Error::other(
        "failed to generate a unique single-edit k-mer for error-correction task",
    ))
}

fn hamming_distance_is_one(left: &[u8], right: &[u8]) -> bool {
    if left.len() != right.len() {
        return false;
    }

    let mut mismatches = 0usize;
    for (left_base, right_base) in left.iter().zip(right.iter()) {
        if left_base != right_base {
            mismatches += 1;
            if mismatches > 1 {
                return false;
            }
        }
    }

    mismatches == 1
}

fn encode_canonical_kmer(sequence: &str) -> io::Result<u64> {
    KmerU64::from_str(sequence)
        .map(|encoded| encoded.canonical().encoded)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid DNA k-mer '{}'", sequence),
            )
        })
}

fn encode_kmer_list(kmers: &[String], _k: usize) -> io::Result<Vec<u64>> {
    kmers
        .iter()
        .map(|kmer| encode_canonical_kmer(kmer))
        .collect()
}

fn encode_trusted_totals(
    trusted_totals: &BTreeMap<String, u32>,
    _k: usize,
) -> io::Result<AHashMap<u64, u32>> {
    trusted_totals
        .iter()
        .map(|(kmer, count)| Ok((encode_canonical_kmer(kmer)?, *count)))
        .collect()
}

fn decode_error_correction_counts(
    counts: &[ErrorCorrectionCountSpec],
    _k: usize,
) -> io::Result<AHashMap<u64, u32>> {
    let mut decoded = AHashMap::with_capacity(counts.len());
    for entry in counts {
        let encoded = encode_canonical_kmer(&entry.kmer)?;
        decoded.insert(encoded, entry.count);
    }
    Ok(decoded)
}

fn validate_contig_extraction_args(
    k: usize,
    component_count: usize,
    primary_reads_per_component: usize,
    alternate_reads_per_component: usize,
) -> io::Result<()> {
    if !(2..=32).contains(&k) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("contig-extraction k must be in 2..=32, got {k}"),
        ));
    }
    if component_count == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "component_count must be greater than zero",
        ));
    }
    if primary_reads_per_component == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "primary_reads_per_component must be greater than zero",
        ));
    }
    if alternate_reads_per_component == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "alternate_reads_per_component must be greater than zero",
        ));
    }
    Ok(())
}

fn build_contig_extraction_task_files(
    k: usize,
    component_count: usize,
    primary_reads_per_component: usize,
    alternate_reads_per_component: usize,
    seed: u64,
    profile: ContigExtractionScenarioProfile,
) -> io::Result<(
    Vec<ContigExtractionCountSpec>,
    ContigExtractionGraphFile,
    ContigExtractionTruthFile,
    ContigExtractionTaskMetadata,
)> {
    let (kmer_counts, valid_kmers, branch_support, expected_contigs, description) = match profile {
        ContigExtractionScenarioProfile::Branching => {
            let fixture = LargeGenomeStageBenchFixture::synthetic_branching(
                k,
                component_count,
                primary_reads_per_component,
                alternate_reads_per_component,
            );
            (
                fixture.cloned_kmer_counts(),
                fixture.valid_kmers(),
                fixture.branch_support_map(),
                fixture.expected_contigs().to_vec(),
                "Synthetic branching graph that rewards high-count seed selection and accurate branch-aware contig extraction."
                    .to_string(),
            )
        }
        ContigExtractionScenarioProfile::RepeatFallback => {
            let (counts, valid_kmers, branch_support, expected_contigs) =
                build_repeat_fallback_contig_task(k, seed)?;
            (
                counts,
                valid_kmers,
                branch_support,
                expected_contigs,
                "Synthetic graph where only repeat-classified seeds survive cleanup, forcing repeat-seed fallback."
                    .to_string(),
            )
        }
        ContigExtractionScenarioProfile::RepeatCompletion => {
            let (counts, valid_kmers, branch_support, expected_contigs) =
                build_repeat_completion_contig_task(k, seed)?;
            (
                counts,
                valid_kmers,
                branch_support,
                expected_contigs,
                "Synthetic graph with a non-repeat component plus a repeat-only component that needs the repeat completion pass."
                    .to_string(),
            )
        }
        ContigExtractionScenarioProfile::RepeatPriority => {
            let (counts, valid_kmers, branch_support, expected_contigs) =
                build_repeat_priority_contig_task(k, seed)?;
            (
                counts,
                valid_kmers,
                branch_support,
                expected_contigs,
                "Synthetic graph where a high-count repeat branch steals shared suffix nodes unless non-repeat seeds are processed first."
                    .to_string(),
            )
        }
    };

    let counts = encode_contig_extraction_counts(&kmer_counts, k);
    let graph = ContigExtractionGraphFile {
        valid_kmers: encode_canonical_kmer_set(&valid_kmers, k),
        branch_support: encode_contig_extraction_branch_support(&branch_support, k),
    };
    let truth = ContigExtractionTruthFile { expected_contigs };
    let metadata = ContigExtractionTaskMetadata {
        component: "contig_extraction".to_string(),
        name: profile.as_str().to_string(),
        description,
        k,
        profile,
        component_count,
        primary_reads_per_component,
        alternate_reads_per_component,
        expected_contigs: truth.expected_contigs.len(),
        seed,
    };

    Ok((counts, graph, truth, metadata))
}

fn build_repeat_fallback_contig_task(
    k: usize,
    seed: u64,
) -> io::Result<(
    AHashMap<u64, u32>,
    AHashSet<u64>,
    AHashMap<(u64, u64), u32>,
    Vec<String>,
)> {
    let root = "A".repeat(k);
    let mid = format!("{}C", "A".repeat(k.saturating_sub(1)));
    let end = format!("{}CC", "A".repeat(k.saturating_sub(2)));
    let expected_contig = format!("{root}CC");

    let mut counts = AHashMap::new();
    let mut valid_kmers = AHashSet::new();
    for sequence in [&root, &mid, &end] {
        let encoded = encode_canonical_kmer(sequence)?;
        counts.insert(encoded, 20);
        valid_kmers.insert(encoded);
    }

    let mut rng = StdRng::seed_from_u64(seed.max(1));
    let mut used_canonical = valid_kmers.iter().copied().collect::<BTreeSet<_>>();
    for decoy in collect_unique_canonical_kmers(k, 24, &mut used_canonical, &mut rng)? {
        counts.insert(encode_canonical_kmer(&decoy)?, 2);
    }

    Ok((counts, valid_kmers, AHashMap::new(), vec![expected_contig]))
}

fn build_repeat_completion_contig_task(
    k: usize,
    seed: u64,
) -> io::Result<(
    AHashMap<u64, u32>,
    AHashSet<u64>,
    AHashMap<(u64, u64), u32>,
    Vec<String>,
)> {
    let left_root = "A".repeat(k);
    let left_next = format!("{}T", "A".repeat(k.saturating_sub(1)));
    let right_root = "C".repeat(k);
    let right_next = format!("{}G", "C".repeat(k.saturating_sub(1)));

    let mut counts = AHashMap::new();
    let mut valid_kmers = AHashSet::new();
    for sequence in [&left_root, &left_next] {
        let encoded = encode_canonical_kmer(sequence)?;
        counts.insert(encoded, 6);
        valid_kmers.insert(encoded);
    }
    for sequence in [&right_root, &right_next] {
        let encoded = encode_canonical_kmer(sequence)?;
        counts.insert(encoded, 20);
        valid_kmers.insert(encoded);
    }

    let mut rng = StdRng::seed_from_u64(seed.max(1) ^ 0xC0DE_CAFE_u64);
    let mut used_canonical = valid_kmers.iter().copied().collect::<BTreeSet<_>>();
    for decoy in collect_unique_canonical_kmers(k, 24, &mut used_canonical, &mut rng)? {
        counts.insert(encode_canonical_kmer(&decoy)?, 2);
    }

    Ok((
        counts,
        valid_kmers,
        AHashMap::new(),
        vec![format!("{left_root}T"), format!("{right_root}G")],
    ))
}

fn build_repeat_priority_contig_task(
    k: usize,
    seed: u64,
) -> io::Result<(
    AHashMap<u64, u32>,
    AHashSet<u64>,
    AHashMap<(u64, u64), u32>,
    Vec<String>,
)> {
    if k < 5 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "repeat_priority profile requires k >= 5",
        ));
    }

    let root = "A".repeat(k);
    let common_tail = contig_repeat_priority_tail(k.saturating_add(3));
    let truth_sequence = format!("{root}T{common_tail}");
    let repeat_sequence = format!("{root}G{common_tail}");
    let truth_windows = contig_sequence_windows(&truth_sequence, k)?;
    let repeat_windows = contig_sequence_windows(&repeat_sequence, k)?;
    let truth_window_set = truth_windows.iter().cloned().collect::<BTreeSet<_>>();
    let repeat_window_set = repeat_windows.iter().cloned().collect::<BTreeSet<_>>();
    let shared_windows = truth_window_set
        .intersection(&repeat_window_set)
        .cloned()
        .collect::<BTreeSet<_>>();

    let mut counts = AHashMap::new();
    let mut valid_kmers = AHashSet::new();
    for window in &truth_windows {
        let count = if window == &root {
            7
        } else if shared_windows.contains(window) {
            5
        } else {
            6
        };
        insert_contig_task_kmer(window, count, &mut counts, &mut valid_kmers)?;
    }
    for window in &repeat_windows {
        let count = if window == &root {
            7
        } else if shared_windows.contains(window) {
            5
        } else {
            20
        };
        insert_contig_task_kmer(window, count, &mut counts, &mut valid_kmers)?;
    }

    let mut branch_support = AHashMap::new();
    // Make the true branch easy to recover when extraction starts from a shared
    // non-repeat seed, while still letting repeat-first seeding capture the
    // shared suffix component.
    branch_support.insert(
        (
            encode_oriented_kmer(&truth_windows[0])?,
            encode_oriented_kmer(&truth_windows[1])?,
        ),
        6,
    );
    branch_support.insert(
        (
            encode_oriented_kmer(&repeat_windows[0])?,
            encode_oriented_kmer(&repeat_windows[1])?,
        ),
        1,
    );

    let mut rng = StdRng::seed_from_u64(seed.max(1) ^ 0xA11C_E551_u64);
    let mut used_canonical = valid_kmers.iter().copied().collect::<BTreeSet<_>>();
    for decoy in collect_unique_canonical_kmers(k, 32, &mut used_canonical, &mut rng)? {
        counts.insert(encode_canonical_kmer(&decoy)?, 2);
    }

    Ok((counts, valid_kmers, branch_support, vec![truth_sequence]))
}

fn insert_contig_task_kmer(
    sequence: &str,
    count: u32,
    counts: &mut AHashMap<u64, u32>,
    valid_kmers: &mut AHashSet<u64>,
) -> io::Result<()> {
    let encoded = encode_canonical_kmer(sequence)?;
    counts.insert(encoded, count);
    valid_kmers.insert(encoded);
    Ok(())
}

fn contig_sequence_windows(sequence: &str, k: usize) -> io::Result<Vec<String>> {
    if sequence.len() < k {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("sequence length {} is smaller than k={k}", sequence.len()),
        ));
    }

    sequence
        .as_bytes()
        .windows(k)
        .map(|window| String::from_utf8(window.to_vec()).map_err(invalid_utf8_error))
        .collect()
}

fn contig_repeat_priority_tail(len: usize) -> String {
    const MOTIF: &[u8] = b"CGTAC";
    let mut tail = String::with_capacity(len);
    for idx in 0..len {
        tail.push(MOTIF[idx % MOTIF.len()] as char);
    }
    tail
}

fn encode_contig_extraction_counts(
    counts: &AHashMap<u64, u32>,
    k: usize,
) -> Vec<ContigExtractionCountSpec> {
    let mut encoded = counts
        .iter()
        .map(|(&kmer, &count)| ContigExtractionCountSpec {
            kmer: decode_encoded_kmer(kmer, k),
            count,
        })
        .collect::<Vec<_>>();
    encoded.sort_by(|left, right| left.kmer.cmp(&right.kmer));
    encoded
}

fn encode_canonical_kmer_set(valid_kmers: &AHashSet<u64>, k: usize) -> Vec<String> {
    let mut kmers = valid_kmers
        .iter()
        .copied()
        .map(|encoded| decode_encoded_kmer(encoded, k))
        .collect::<Vec<_>>();
    kmers.sort();
    kmers
}

fn encode_contig_extraction_branch_support(
    branch_support: &AHashMap<(u64, u64), u32>,
    k: usize,
) -> Vec<ContigExtractionBranchSupportSpec> {
    let mut encoded = branch_support
        .iter()
        .map(
            |(&(from, to), &support)| ContigExtractionBranchSupportSpec {
                from_kmer: decode_encoded_kmer(from, k),
                to_kmer: decode_encoded_kmer(to, k),
                support,
            },
        )
        .collect::<Vec<_>>();
    encoded.sort_by(|left, right| {
        left.from_kmer
            .cmp(&right.from_kmer)
            .then_with(|| left.to_kmer.cmp(&right.to_kmer))
    });
    encoded
}

fn decode_contig_extraction_counts(
    counts: &[ContigExtractionCountSpec],
    _k: usize,
) -> io::Result<AHashMap<u64, u32>> {
    let mut decoded = AHashMap::with_capacity(counts.len());
    for entry in counts {
        decoded.insert(encode_canonical_kmer(&entry.kmer)?, entry.count);
    }
    Ok(decoded)
}

fn decode_contig_extraction_valid_kmers(valid_kmers: &[String]) -> io::Result<AHashSet<u64>> {
    valid_kmers
        .iter()
        .map(|sequence| encode_canonical_kmer(sequence))
        .collect()
}

fn decode_contig_extraction_branch_support(
    branch_support: &[ContigExtractionBranchSupportSpec],
) -> io::Result<AHashMap<(u64, u64), u32>> {
    let mut decoded = AHashMap::with_capacity(branch_support.len());
    for entry in branch_support {
        decoded.insert(
            (
                encode_oriented_kmer(&entry.from_kmer)?,
                encode_oriented_kmer(&entry.to_kmer)?,
            ),
            entry.support,
        );
    }
    Ok(decoded)
}

fn decode_encoded_kmer(encoded: u64, k: usize) -> String {
    KmerU64 {
        encoded,
        len: k as u8,
    }
    .decode()
}

fn encode_oriented_kmer(sequence: &str) -> io::Result<u64> {
    KmerU64::from_str(sequence)
        .map(|encoded| encoded.encoded)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid oriented DNA k-mer '{}'", sequence),
            )
        })
}

fn canonical_contig_set(sequences: &[String]) -> BTreeSet<String> {
    sequences
        .iter()
        .map(|sequence| canonicalize_sequence(sequence))
        .collect()
}

fn canonicalize_sequence(sequence: &str) -> String {
    let reverse = crate::kmer::kmer::reverse_complement(sequence);
    if sequence <= reverse.as_str() {
        sequence.to_string()
    } else {
        reverse
    }
}

fn canonical_kmer_set_from_sequences(sequences: &[String], k: usize) -> io::Result<AHashSet<u64>> {
    let mut kmers = AHashSet::new();
    for sequence in sequences {
        if sequence.len() < k {
            continue;
        }
        for window in sequence.as_bytes().windows(k) {
            let encoded = KmerU64::from_slice(window).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid DNA sequence '{}'", sequence),
                )
            })?;
            kmers.insert(encoded.canonical().encoded);
        }
    }
    Ok(kmers)
}

fn evaluate_read_mapping_bench(
    task_path: &Path,
    minimizer_k: usize,
    minimizer_w: usize,
    min_primary_matches: usize,
    min_scaffold_matches: usize,
    position_tolerance: usize,
) -> io::Result<ReadMappingBenchmarkReport> {
    let task_dirs = discover_task_dirs(task_path)?;
    let tasks = task_dirs
        .par_iter()
        .map(|task_dir| {
            evaluate_single_read_mapping_task(
                task_dir,
                minimizer_k,
                minimizer_w,
                min_primary_matches,
                min_scaffold_matches,
                position_tolerance,
            )
        })
        .collect::<io::Result<Vec<_>>>()?;

    let summary = summarize_reports(&tasks);
    Ok(ReadMappingBenchmarkReport {
        component: "read_mapping",
        minimizer_k,
        minimizer_w,
        min_primary_matches,
        min_scaffold_matches,
        position_tolerance,
        tasks,
        summary,
    })
}

fn discover_task_dirs(task_path: &Path) -> io::Result<Vec<PathBuf>> {
    if is_task_dir(task_path) {
        return Ok(vec![task_path.to_path_buf()]);
    }

    let mut discovered = Vec::new();
    for entry in fs::read_dir(task_path)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() && is_task_dir(&path) {
            discovered.push(path);
        }
    }
    discovered.sort();

    if discovered.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("no read-mapping tasks found under {}", task_path.display()),
        ));
    }

    Ok(discovered)
}

fn is_task_dir(task_dir: &Path) -> bool {
    resolve_contigs_path(task_dir).is_ok()
        && resolve_reads1_path(task_dir).is_ok()
        && resolve_truth_path(task_dir).is_ok()
}

fn evaluate_single_read_mapping_task(
    task_dir: &Path,
    minimizer_k: usize,
    minimizer_w: usize,
    min_primary_matches: usize,
    min_scaffold_matches: usize,
    position_tolerance: usize,
) -> io::Result<ReadMappingTaskReport> {
    let contigs_path = resolve_contigs_path(task_dir)?;
    let reads1_path = resolve_reads1_path(task_dir)?;
    let reads2_path = resolve_reads2_path(task_dir);
    let truth_path = resolve_truth_path(task_dir)?;
    let metadata = read_task_metadata(task_dir)?;
    let task_name = metadata
        .as_ref()
        .and_then(|task| {
            if task.name.trim().is_empty() {
                None
            } else {
                Some(task.name.clone())
            }
        })
        .unwrap_or_else(|| {
            task_dir
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("read_mapping_task")
                .to_string()
        });

    let contigs = read_contigs(path_to_str(&contigs_path)?)?;
    let mut header_to_id = AHashMap::with_capacity(contigs.len());
    for (contig_id, (header, _seq)) in contigs.iter().enumerate() {
        header_to_id.insert(header.clone(), contig_id);
    }

    let truth = load_truth(&truth_path, &header_to_id)?;

    let mut index = MinimizerIndex::new(minimizer_k, minimizer_w);
    for (contig_id, (_header, sequence)) in contigs.iter().enumerate() {
        index.add_contig(contig_id, sequence);
    }

    let start = Instant::now();
    let mut counters = MappingCounters::default();
    let mut observed = AHashMap::with_capacity(truth.len());
    evaluate_fastq_file(
        &reads1_path,
        &index,
        &truth,
        min_primary_matches,
        min_scaffold_matches,
        position_tolerance,
        &mut counters,
        &mut observed,
    )?;
    if let Some(reads2_path) = reads2_path {
        evaluate_fastq_file(
            &reads2_path,
            &index,
            &truth,
            min_primary_matches,
            min_scaffold_matches,
            position_tolerance,
            &mut counters,
            &mut observed,
        )?;
    }
    let runtime_seconds = start.elapsed().as_secs_f64();
    let pair_counters = paired_pair_mapping_counters(&truth, &observed, metadata.as_ref());

    if counters.reads_evaluated != truth.len() as u64 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "truth/read count mismatch for {}: evaluated {} reads but loaded {} truth entries",
                task_dir.display(),
                counters.reads_evaluated,
                truth.len()
            ),
        ));
    }

    Ok(build_task_report(
        task_dir,
        task_name,
        counters,
        pair_counters,
        runtime_seconds,
    ))
}

fn evaluate_fastq_file(
    reads_path: &Path,
    index: &MinimizerIndex,
    truth: &AHashMap<String, ReadExpectation>,
    min_primary_matches: usize,
    min_scaffold_matches: usize,
    position_tolerance: usize,
    counters: &mut MappingCounters,
    observed: &mut AHashMap<String, ObservedReadMapping>,
) -> io::Result<()> {
    let reader = try_open_fastq(path_to_str(reads_path)?)?;
    let mut hits = AHashMap::new();
    let mut contig_hits = AHashMap::new();
    let mut reverse_scratch = Vec::new();

    for record in stream_fastq_records_checked(reader) {
        let record = record?;
        let read_id = normalize_read_id(&record.header).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "FASTQ read header missing ID")
        })?;
        let expected = truth.get(read_id).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "read '{}' in {} is missing from truth.tsv",
                    read_id,
                    reads_path.display()
                ),
            )
        })?;

        let (primary_hit, scaffold_hit) = index.map_read_with_scaffold_support_thresholds(
            record.sequence.as_bytes(),
            &mut hits,
            &mut contig_hits,
            &mut reverse_scratch,
            min_primary_matches,
            min_scaffold_matches,
        );
        counters.observe(expected, primary_hit, scaffold_hit, position_tolerance);
        observed.insert(
            read_id.to_string(),
            ObservedReadMapping {
                scaffold_hit: scaffold_hit
                    .map(|(contig_id, _support, is_reverse)| (contig_id, is_reverse)),
            },
        );
    }

    Ok(())
}

fn build_task_report(
    task_dir: &Path,
    task_name: String,
    counters: MappingCounters,
    pair_counters: PairMappingCounters,
    runtime_seconds: f64,
) -> ReadMappingTaskReport {
    let reads_per_second = if runtime_seconds > 0.0 {
        counters.reads_evaluated as f64 / runtime_seconds
    } else {
        0.0
    };
    let primary_mapped_rate = ratio(counters.primary_mapped, counters.primary_expected);
    let primary_exact_rate = ratio(counters.primary_exact, counters.primary_expected);
    let primary_near_rate = ratio(counters.primary_near, counters.primary_expected);
    let primary_contig_rate = ratio(counters.primary_contig_correct, counters.primary_expected);
    let primary_orientation_rate = ratio(
        counters.primary_orientation_correct,
        counters.primary_expected,
    );
    let primary_unexpected_rate = ratio(
        counters.primary_unexpected_mapped,
        counters.primary_unexpected_total,
    );
    let primary_specificity_rate = primary_unexpected_rate.map(|rate| 1.0 - rate);
    let scaffold_mapped_rate = ratio(counters.scaffold_mapped, counters.scaffold_expected);
    let scaffold_exact_rate = ratio(counters.scaffold_exact, counters.scaffold_expected);
    let scaffold_unexpected_rate = ratio(
        counters.scaffold_unexpected_mapped,
        counters.scaffold_unexpected_total,
    );
    let scaffold_specificity_rate = scaffold_unexpected_rate.map(|rate| 1.0 - rate);
    let pair_mapped_rate = ratio(pair_counters.pair_mapped, pair_counters.pair_expected);
    let pair_cross_contig_rate = ratio(
        pair_counters.pair_cross_contig_correct,
        pair_counters.pair_expected,
    );
    let pair_unexpected_rate = ratio(
        pair_counters.pair_unexpected_mapped,
        pair_counters.pair_unexpected_total,
    );
    let pair_specificity_rate = pair_unexpected_rate.map(|rate| 1.0 - rate);
    let primary_position_mae = ratio_f64(
        counters.primary_position_abs_error_sum as f64,
        counters.primary_position_abs_error_count as f64,
    );
    let score = mapping_score(
        primary_near_rate,
        primary_exact_rate,
        scaffold_exact_rate,
        primary_mapped_rate,
        primary_specificity_rate,
        scaffold_specificity_rate,
        pair_cross_contig_rate,
        pair_specificity_rate,
    );

    ReadMappingTaskReport {
        task_dir: task_dir.display().to_string(),
        task_name,
        reads_evaluated: counters.reads_evaluated,
        primary_expected: counters.primary_expected,
        primary_mapped: counters.primary_mapped,
        primary_exact: counters.primary_exact,
        primary_near: counters.primary_near,
        primary_contig_correct: counters.primary_contig_correct,
        primary_orientation_correct: counters.primary_orientation_correct,
        primary_unexpected_total: counters.primary_unexpected_total,
        primary_unexpected_mapped: counters.primary_unexpected_mapped,
        scaffold_expected: counters.scaffold_expected,
        scaffold_mapped: counters.scaffold_mapped,
        scaffold_exact: counters.scaffold_exact,
        scaffold_unexpected_total: counters.scaffold_unexpected_total,
        scaffold_unexpected_mapped: counters.scaffold_unexpected_mapped,
        pair_expected: pair_counters.pair_expected,
        pair_mapped: pair_counters.pair_mapped,
        pair_cross_contig_correct: pair_counters.pair_cross_contig_correct,
        pair_unexpected_total: pair_counters.pair_unexpected_total,
        pair_unexpected_mapped: pair_counters.pair_unexpected_mapped,
        pair_cross_contig_errors: pair_counters.pair_cross_contig_errors,
        primary_mapped_rate,
        primary_exact_rate,
        primary_near_rate,
        primary_contig_rate,
        primary_orientation_rate,
        primary_unexpected_rate,
        primary_specificity_rate,
        scaffold_mapped_rate,
        scaffold_exact_rate,
        scaffold_unexpected_rate,
        scaffold_specificity_rate,
        pair_mapped_rate,
        pair_cross_contig_rate,
        pair_unexpected_rate,
        pair_specificity_rate,
        primary_position_mae,
        primary_position_samples: counters.primary_position_abs_error_count,
        runtime_seconds,
        reads_per_second,
        score,
    }
}

fn summarize_reports(tasks: &[ReadMappingTaskReport]) -> ReadMappingAggregateReport {
    let mut counters = MappingCounters::default();
    let mut pair_counters = PairMappingCounters::default();
    let mut runtime_seconds = 0.0;

    for task in tasks {
        counters.merge(&MappingCounters {
            reads_evaluated: task.reads_evaluated,
            primary_expected: task.primary_expected,
            primary_mapped: task.primary_mapped,
            primary_exact: task.primary_exact,
            primary_near: task.primary_near,
            primary_contig_correct: task.primary_contig_correct,
            primary_orientation_correct: task.primary_orientation_correct,
            primary_unexpected_total: task.primary_unexpected_total,
            primary_unexpected_mapped: task.primary_unexpected_mapped,
            scaffold_expected: task.scaffold_expected,
            scaffold_mapped: task.scaffold_mapped,
            scaffold_exact: task.scaffold_exact,
            scaffold_unexpected_total: task.scaffold_unexpected_total,
            scaffold_unexpected_mapped: task.scaffold_unexpected_mapped,
            primary_position_abs_error_sum: task
                .primary_position_mae
                .map(|mae| (mae * task.primary_position_samples as f64).round() as u64)
                .unwrap_or(0),
            primary_position_abs_error_count: task.primary_position_samples,
        });
        pair_counters.merge(&PairMappingCounters {
            pair_expected: task.pair_expected,
            pair_mapped: task.pair_mapped,
            pair_cross_contig_correct: task.pair_cross_contig_correct,
            pair_unexpected_total: task.pair_unexpected_total,
            pair_unexpected_mapped: task.pair_unexpected_mapped,
            pair_cross_contig_errors: task.pair_cross_contig_errors,
        });
        runtime_seconds += task.runtime_seconds;
    }

    let reads_per_second = if runtime_seconds > 0.0 {
        counters.reads_evaluated as f64 / runtime_seconds
    } else {
        0.0
    };
    let primary_mapped_rate = ratio(counters.primary_mapped, counters.primary_expected);
    let primary_exact_rate = ratio(counters.primary_exact, counters.primary_expected);
    let primary_near_rate = ratio(counters.primary_near, counters.primary_expected);
    let primary_contig_rate = ratio(counters.primary_contig_correct, counters.primary_expected);
    let primary_orientation_rate = ratio(
        counters.primary_orientation_correct,
        counters.primary_expected,
    );
    let primary_unexpected_rate = ratio(
        counters.primary_unexpected_mapped,
        counters.primary_unexpected_total,
    );
    let primary_specificity_rate = primary_unexpected_rate.map(|rate| 1.0 - rate);
    let scaffold_mapped_rate = ratio(counters.scaffold_mapped, counters.scaffold_expected);
    let scaffold_exact_rate = ratio(counters.scaffold_exact, counters.scaffold_expected);
    let scaffold_unexpected_rate = ratio(
        counters.scaffold_unexpected_mapped,
        counters.scaffold_unexpected_total,
    );
    let scaffold_specificity_rate = scaffold_unexpected_rate.map(|rate| 1.0 - rate);
    let pair_mapped_rate = ratio(pair_counters.pair_mapped, pair_counters.pair_expected);
    let pair_cross_contig_rate = ratio(
        pair_counters.pair_cross_contig_correct,
        pair_counters.pair_expected,
    );
    let pair_unexpected_rate = ratio(
        pair_counters.pair_unexpected_mapped,
        pair_counters.pair_unexpected_total,
    );
    let pair_specificity_rate = pair_unexpected_rate.map(|rate| 1.0 - rate);
    let primary_position_mae = ratio_f64(
        counters.primary_position_abs_error_sum as f64,
        counters.primary_position_abs_error_count as f64,
    );
    let score = mapping_score(
        primary_near_rate,
        primary_exact_rate,
        scaffold_exact_rate,
        primary_mapped_rate,
        primary_specificity_rate,
        scaffold_specificity_rate,
        pair_cross_contig_rate,
        pair_specificity_rate,
    );

    ReadMappingAggregateReport {
        task_count: tasks.len(),
        reads_evaluated: counters.reads_evaluated,
        primary_expected: counters.primary_expected,
        primary_mapped: counters.primary_mapped,
        primary_exact: counters.primary_exact,
        primary_near: counters.primary_near,
        primary_contig_correct: counters.primary_contig_correct,
        primary_orientation_correct: counters.primary_orientation_correct,
        primary_unexpected_total: counters.primary_unexpected_total,
        primary_unexpected_mapped: counters.primary_unexpected_mapped,
        scaffold_expected: counters.scaffold_expected,
        scaffold_mapped: counters.scaffold_mapped,
        scaffold_exact: counters.scaffold_exact,
        scaffold_unexpected_total: counters.scaffold_unexpected_total,
        scaffold_unexpected_mapped: counters.scaffold_unexpected_mapped,
        pair_expected: pair_counters.pair_expected,
        pair_mapped: pair_counters.pair_mapped,
        pair_cross_contig_correct: pair_counters.pair_cross_contig_correct,
        pair_unexpected_total: pair_counters.pair_unexpected_total,
        pair_unexpected_mapped: pair_counters.pair_unexpected_mapped,
        pair_cross_contig_errors: pair_counters.pair_cross_contig_errors,
        primary_mapped_rate,
        primary_exact_rate,
        primary_near_rate,
        primary_contig_rate,
        primary_orientation_rate,
        primary_unexpected_rate,
        primary_specificity_rate,
        scaffold_mapped_rate,
        scaffold_exact_rate,
        scaffold_unexpected_rate,
        scaffold_specificity_rate,
        pair_mapped_rate,
        pair_cross_contig_rate,
        pair_unexpected_rate,
        pair_specificity_rate,
        primary_position_mae,
        primary_position_samples: counters.primary_position_abs_error_count,
        runtime_seconds,
        reads_per_second,
        score,
    }
}

fn load_truth(
    truth_path: &Path,
    header_to_id: &AHashMap<String, usize>,
) -> io::Result<AHashMap<String, ReadExpectation>> {
    let reader = BufReader::new(File::open(truth_path)?);
    let mut lines = reader.lines();

    let header_line = lines.next().transpose()?.ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("truth file {} is empty", truth_path.display()),
        )
    })?;
    let header_map = parse_truth_header(&header_line)?;
    let mut truth = AHashMap::new();

    for (line_index, line) in lines.enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        let get = |name: &str| -> io::Result<&str> {
            let idx = *header_map.get(name).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("truth file missing required column '{name}'"),
                )
            })?;
            fields.get(idx).copied().ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "truth row {} in {} is missing field '{}'",
                        line_index + 2,
                        truth_path.display(),
                        name
                    ),
                )
            })
        };

        let read_id = get("read_id")?.trim();
        if read_id.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("truth row {} has empty read_id", line_index + 2),
            ));
        }

        let primary = parse_expected_primary_hit(
            get("contig")?,
            get("start")?,
            get("is_reverse")?,
            header_to_id,
            truth_path,
            line_index + 2,
        )?;
        let scaffold = parse_expected_scaffold_hit(
            get("scaffold_contig")?,
            get("scaffold_is_reverse")?,
            header_to_id,
            truth_path,
            line_index + 2,
        )?;

        if truth
            .insert(read_id.to_string(), ReadExpectation { primary, scaffold })
            .is_some()
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "duplicate truth row for read '{}' in {}",
                    read_id,
                    truth_path.display()
                ),
            ));
        }
    }

    Ok(truth)
}

fn parse_truth_header(header_line: &str) -> io::Result<AHashMap<String, usize>> {
    let mut header_map = AHashMap::new();
    for (idx, field) in header_line.split('\t').enumerate() {
        header_map.insert(field.trim().to_string(), idx);
    }

    for required in [
        "read_id",
        "contig",
        "start",
        "is_reverse",
        "scaffold_contig",
        "scaffold_is_reverse",
    ] {
        if !header_map.contains_key(required) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("truth header missing required column '{required}'"),
            ));
        }
    }

    Ok(header_map)
}

fn parse_expected_primary_hit(
    contig: &str,
    start: &str,
    is_reverse: &str,
    header_to_id: &AHashMap<String, usize>,
    truth_path: &Path,
    line_number: usize,
) -> io::Result<Option<ExpectedPrimaryHit>> {
    let contig = contig.trim();
    if contig.is_empty() {
        return Ok(None);
    }

    let contig_id = header_to_id.get(contig).copied().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "truth row {} references unknown contig '{}' in {}",
                line_number,
                contig,
                truth_path.display()
            ),
        )
    })?;
    let start = start.trim().parse::<usize>().map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "truth row {} has invalid start '{}' in {}: {}",
                line_number,
                start.trim(),
                truth_path.display(),
                err
            ),
        )
    })?;

    Ok(Some(ExpectedPrimaryHit {
        contig_id,
        start,
        is_reverse: parse_boolish(is_reverse.trim(), truth_path, line_number)?,
    }))
}

fn parse_expected_scaffold_hit(
    contig: &str,
    is_reverse: &str,
    header_to_id: &AHashMap<String, usize>,
    truth_path: &Path,
    line_number: usize,
) -> io::Result<Option<ExpectedScaffoldHit>> {
    let contig = contig.trim();
    if contig.is_empty() {
        return Ok(None);
    }

    let contig_id = header_to_id.get(contig).copied().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "truth row {} references unknown scaffold contig '{}' in {}",
                line_number,
                contig,
                truth_path.display()
            ),
        )
    })?;

    Ok(Some(ExpectedScaffoldHit {
        contig_id,
        is_reverse: parse_boolish(is_reverse.trim(), truth_path, line_number)?,
    }))
}

fn parse_boolish(raw: &str, truth_path: &Path, line_number: usize) -> io::Result<bool> {
    match raw {
        "1" | "true" | "True" | "TRUE" | "forward_rc" => Ok(true),
        "0" | "false" | "False" | "FALSE" | "" => Ok(false),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "truth row {} in {} has invalid boolean '{}'",
                line_number,
                truth_path.display(),
                raw
            ),
        )),
    }
}

fn paired_pair_mapping_counters(
    truth: &AHashMap<String, ReadExpectation>,
    observed: &AHashMap<String, ObservedReadMapping>,
    metadata: Option<&ReadMappingTaskMetadata>,
) -> PairMappingCounters {
    if !metadata.map(|task| task.paired).unwrap_or(false) {
        return PairMappingCounters::default();
    }

    let mut pairs = AHashMap::<String, PairObservation>::default();
    for (read_id, expectation) in truth {
        let Some((template_id, mate_index)) = split_pair_read_id(read_id) else {
            continue;
        };
        let pair = pairs.entry(template_id.to_string()).or_default();
        let observed_scaffold_hit = observed.get(read_id).and_then(|hit| hit.scaffold_hit);
        match mate_index {
            0 => {
                pair.mate1_expected = expectation.scaffold.clone();
                pair.mate1_scaffold_hit = observed_scaffold_hit;
            }
            1 => {
                pair.mate2_expected = expectation.scaffold.clone();
                pair.mate2_scaffold_hit = observed_scaffold_hit;
            }
            _ => {}
        }
    }

    let mut counters = PairMappingCounters::default();
    for observation in pairs.values() {
        counters.observe(observation);
    }
    counters
}

fn split_pair_read_id(read_id: &str) -> Option<(&str, usize)> {
    if let Some(template_id) = read_id.strip_suffix("/1") {
        return Some((template_id, 0));
    }
    if let Some(template_id) = read_id.strip_suffix("/2") {
        return Some((template_id, 1));
    }
    None
}

fn write_contigs_fasta(path: &Path, contigs: &[(String, Vec<u8>)]) -> io::Result<()> {
    let mut writer = FastaWriter::try_new(path_to_str(path)?)?;
    for (header, sequence) in contigs {
        let sequence = std::str::from_utf8(sequence).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("generated contig sequence is not UTF-8 DNA: {err}"),
            )
        })?;
        writer.write_record(header, sequence)?;
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn write_synthetic_reads(
    contigs: &[(String, Vec<u8>)],
    repeat: &[u8],
    read_len: usize,
    reads: usize,
    paired: bool,
    insert_size: usize,
    error_rate: f64,
    decoy_rate: f64,
    ambiguous_repeat_decoy_rate: f64,
    reads1_path: &Path,
    reads2_path: Option<&Path>,
    truth_path: &Path,
    rng: &mut StdRng,
) -> io::Result<()> {
    let mut writer1 = FastqWriter::try_new(path_to_str(reads1_path)?)?;
    let mut writer2 = reads2_path
        .map(|path| FastqWriter::try_new(path_to_str(path)?))
        .transpose()?;
    let mut truth_writer = File::create(truth_path)?;

    writeln!(
        truth_writer,
        "read_id\tcontig\tstart\tis_reverse\tscaffold_contig\tscaffold_is_reverse"
    )?;

    let repeat_positions = contigs
        .iter()
        .map(|(_, sequence)| {
            find_subsequence(sequence, repeat).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "generated repeat was not found in synthetic contig",
                )
            })
        })
        .collect::<io::Result<Vec<_>>>()?;

    for read_index in 0..reads {
        let contig_id = rng.gen_range(0..contigs.len());
        let (contig_name, contig_sequence) = &contigs[contig_id];
        let repeat_start = repeat_positions[contig_id];

        if paired {
            let fragment_len = insert_size.max(read_len * 2);
            let start1 = sample_fragment_start(
                contig_sequence.len(),
                fragment_len,
                repeat_start,
                repeat.len(),
                rng,
            )?;
            let start2 = start1 + fragment_len - read_len;

            let read1_id = format!("read_{read_index:06}/1");
            let read2_id = format!("read_{read_index:06}/2");
            let mut read1_seq = contig_sequence[start1..start1 + read_len].to_vec();
            let mut read2_seq = reverse_complement(&contig_sequence[start2..start2 + read_len]);
            inject_errors(&mut read1_seq, error_rate, rng);
            inject_errors(&mut read2_seq, error_rate, rng);

            writer1.write_record(&FastqRecord {
                header: format!("@{read1_id}"),
                sequence: String::from_utf8(read1_seq).map_err(invalid_utf8_error)?,
                plus: "+".to_string(),
                quality: "I".repeat(read_len),
            })?;
            if let Some(writer2) = writer2.as_mut() {
                writer2.write_record(&FastqRecord {
                    header: format!("@{read2_id}"),
                    sequence: String::from_utf8(read2_seq).map_err(invalid_utf8_error)?,
                    plus: "+".to_string(),
                    quality: "I".repeat(read_len),
                })?;
            }

            write_truth_row(&mut truth_writer, &read1_id, contig_name, start1, false)?;
            write_truth_row(&mut truth_writer, &read2_id, contig_name, start2, true)?;
        } else {
            let start = sample_fragment_start(
                contig_sequence.len(),
                read_len,
                repeat_start,
                repeat.len(),
                rng,
            )?;
            let is_reverse = rng.gen_bool(0.5);
            let mut read_seq = contig_sequence[start..start + read_len].to_vec();
            if is_reverse {
                read_seq = reverse_complement(&read_seq);
            }
            inject_errors(&mut read_seq, error_rate, rng);
            let read_id = format!("read_{read_index:06}");
            writer1.write_record(&FastqRecord {
                header: format!("@{read_id}"),
                sequence: String::from_utf8(read_seq).map_err(invalid_utf8_error)?,
                plus: "+".to_string(),
                quality: "I".repeat(read_len),
            })?;
            write_truth_row(&mut truth_writer, &read_id, contig_name, start, is_reverse)?;
        }
    }

    let decoy_templates = scaled_extra_templates(reads, decoy_rate);
    for decoy_index in 0..decoy_templates {
        if paired {
            let read1_id = format!("decoy_{decoy_index:06}/1");
            let read2_id = format!("decoy_{decoy_index:06}/2");
            let mut read1_seq = random_dna(read_len, rng);
            let mut read2_seq = random_dna(read_len, rng);
            inject_errors(&mut read1_seq, error_rate, rng);
            inject_errors(&mut read2_seq, error_rate, rng);

            writer1.write_record(&FastqRecord {
                header: format!("@{read1_id}"),
                sequence: String::from_utf8(read1_seq).map_err(invalid_utf8_error)?,
                plus: "+".to_string(),
                quality: "I".repeat(read_len),
            })?;
            if let Some(writer2) = writer2.as_mut() {
                writer2.write_record(&FastqRecord {
                    header: format!("@{read2_id}"),
                    sequence: String::from_utf8(read2_seq).map_err(invalid_utf8_error)?,
                    plus: "+".to_string(),
                    quality: "I".repeat(read_len),
                })?;
            }

            write_unmapped_truth_row(&mut truth_writer, &read1_id)?;
            write_unmapped_truth_row(&mut truth_writer, &read2_id)?;
        } else {
            let read_id = format!("decoy_{decoy_index:06}");
            let mut read_seq = random_dna(read_len, rng);
            inject_errors(&mut read_seq, error_rate, rng);
            writer1.write_record(&FastqRecord {
                header: format!("@{read_id}"),
                sequence: String::from_utf8(read_seq).map_err(invalid_utf8_error)?,
                plus: "+".to_string(),
                quality: "I".repeat(read_len),
            })?;
            write_unmapped_truth_row(&mut truth_writer, &read_id)?;
        }
    }

    let ambiguous_repeat_templates = scaled_extra_templates(reads, ambiguous_repeat_decoy_rate);
    for decoy_index in 0..ambiguous_repeat_templates {
        if paired {
            let read1_id = format!("repeat_decoy_{decoy_index:06}/1");
            let read2_id = format!("repeat_decoy_{decoy_index:06}/2");
            let mut read1_seq = sample_ambiguous_repeat_read(repeat, read_len, rng);
            let mut read2_seq =
                reverse_complement(&sample_ambiguous_repeat_read(repeat, read_len, rng));
            inject_errors(&mut read1_seq, error_rate, rng);
            inject_errors(&mut read2_seq, error_rate, rng);

            writer1.write_record(&FastqRecord {
                header: format!("@{read1_id}"),
                sequence: String::from_utf8(read1_seq).map_err(invalid_utf8_error)?,
                plus: "+".to_string(),
                quality: "I".repeat(read_len),
            })?;
            if let Some(writer2) = writer2.as_mut() {
                writer2.write_record(&FastqRecord {
                    header: format!("@{read2_id}"),
                    sequence: String::from_utf8(read2_seq).map_err(invalid_utf8_error)?,
                    plus: "+".to_string(),
                    quality: "I".repeat(read_len),
                })?;
            }

            write_unmapped_truth_row(&mut truth_writer, &read1_id)?;
            write_unmapped_truth_row(&mut truth_writer, &read2_id)?;
        } else {
            let read_id = format!("repeat_decoy_{decoy_index:06}");
            let mut read_seq = sample_ambiguous_repeat_read(repeat, read_len, rng);
            if rng.gen_bool(0.5) {
                read_seq = reverse_complement(&read_seq);
            }
            inject_errors(&mut read_seq, error_rate, rng);
            writer1.write_record(&FastqRecord {
                header: format!("@{read_id}"),
                sequence: String::from_utf8(read_seq).map_err(invalid_utf8_error)?,
                plus: "+".to_string(),
                quality: "I".repeat(read_len),
            })?;
            write_unmapped_truth_row(&mut truth_writer, &read_id)?;
        }
    }

    Ok(())
}

fn write_truth_row(
    writer: &mut File,
    read_id: &str,
    contig_name: &str,
    start: usize,
    is_reverse: bool,
) -> io::Result<()> {
    writeln!(
        writer,
        "{read_id}\t{contig_name}\t{start}\t{is_reverse}\t{contig_name}\t{is_reverse}"
    )
}

fn write_unmapped_truth_row(writer: &mut File, read_id: &str) -> io::Result<()> {
    writeln!(writer, "{read_id}\t\t\t\t\t")
}

fn scaled_extra_templates(reads: usize, rate: f64) -> usize {
    if rate <= 0.0 {
        0
    } else {
        ((reads as f64) * rate).round() as usize
    }
}

fn total_generated_reads(
    reads: usize,
    paired: bool,
    decoy_rate: f64,
    ambiguous_repeat_decoy_rate: f64,
) -> usize {
    let templates = reads
        + scaled_extra_templates(reads, decoy_rate)
        + scaled_extra_templates(reads, ambiguous_repeat_decoy_rate);
    if paired {
        templates * 2
    } else {
        templates
    }
}

fn sample_ambiguous_repeat_read(repeat: &[u8], read_len: usize, rng: &mut StdRng) -> Vec<u8> {
    if repeat.is_empty() || read_len == 0 {
        return Vec::new();
    }

    let start = rng.gen_range(0..repeat.len());
    (0..read_len)
        .map(|offset| repeat[(start + offset) % repeat.len()])
        .collect()
}

fn build_synthetic_contigs(
    num_contigs: usize,
    contig_len: usize,
    repeat: &[u8],
    rng: &mut StdRng,
) -> Vec<(String, Vec<u8>)> {
    let mut contigs = Vec::with_capacity(num_contigs);
    let max_repeat_start = contig_len - repeat.len();
    let base_repeat_start = contig_len / 3;

    for contig_id in 0..num_contigs {
        let mut sequence = random_dna(contig_len, rng);
        let repeat_start = (base_repeat_start + contig_id * 17).min(max_repeat_start);
        sequence[repeat_start..repeat_start + repeat.len()].copy_from_slice(repeat);
        contigs.push((format!("contig_{contig_id}"), sequence));
    }

    contigs
}

fn sample_fragment_start(
    contig_len: usize,
    fragment_len: usize,
    repeat_start: usize,
    repeat_len: usize,
    rng: &mut StdRng,
) -> io::Result<usize> {
    if contig_len < fragment_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("fragment length {fragment_len} exceeds contig length {contig_len}"),
        ));
    }

    let max_start = contig_len - fragment_len;
    let unique_anchor = fragment_len.min(24) / 2;
    let repeat_end = repeat_start + repeat_len;
    let left_min = repeat_start.saturating_sub(fragment_len.saturating_sub(unique_anchor));
    let left_max = repeat_start.saturating_sub(unique_anchor).min(max_start);
    let right_min = repeat_end.saturating_sub(fragment_len.saturating_sub(unique_anchor));
    let right_max = repeat_end.saturating_sub(unique_anchor).min(max_start);

    let preferred = [(left_min, left_max), (right_min.min(max_start), right_max)];
    let preferred_ranges = preferred
        .iter()
        .copied()
        .filter(|(start, end)| start <= end)
        .collect::<Vec<_>>();

    if !preferred_ranges.is_empty() && rng.gen_bool(0.65) {
        let (start, end) = preferred_ranges[rng.gen_range(0..preferred_ranges.len())];
        Ok(rng.gen_range(start..=end))
    } else {
        let unique_ranges = [
            (0, repeat_start.saturating_sub(fragment_len)),
            (repeat_end.min(max_start), max_start),
        ]
        .into_iter()
        .filter(|(start, end)| start <= end)
        .collect::<Vec<_>>();

        if !unique_ranges.is_empty() {
            let (start, end) = unique_ranges[rng.gen_range(0..unique_ranges.len())];
            Ok(rng.gen_range(start..=end))
        } else if !preferred_ranges.is_empty() {
            let (start, end) = preferred_ranges[rng.gen_range(0..preferred_ranges.len())];
            Ok(rng.gen_range(start..=end))
        } else {
            Ok(rng.gen_range(0..=max_start))
        }
    }
}

fn validate_prepare_args(
    paired: bool,
    num_contigs: usize,
    contig_len: usize,
    read_len: usize,
    reads: usize,
    insert_size: usize,
    repeat_len: usize,
    error_rate: f64,
    decoy_rate: f64,
    ambiguous_repeat_decoy_rate: f64,
) -> io::Result<()> {
    if num_contigs == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "num_contigs must be greater than zero",
        ));
    }
    if reads == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "reads must be greater than zero",
        ));
    }
    if read_len < DEFAULT_MINIMIZER_K {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("read_len must be at least {DEFAULT_MINIMIZER_K}"),
        ));
    }
    if repeat_len == 0 || repeat_len >= contig_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "repeat_len must be greater than zero and smaller than contig_len",
        ));
    }
    if contig_len <= read_len + repeat_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "contig_len must leave non-repeat flanks for the requested read_len and repeat_len",
        ));
    }
    if paired && insert_size < read_len * 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "insert_size must be at least 2 * read_len for paired tasks",
        ));
    }
    if error_rate.is_sign_negative() || error_rate >= 1.0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "error_rate must be in [0.0, 1.0)",
        ));
    }
    if decoy_rate.is_sign_negative() || decoy_rate > 1.0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "decoy_rate must be in [0.0, 1.0]",
        ));
    }
    if ambiguous_repeat_decoy_rate.is_sign_negative() || ambiguous_repeat_decoy_rate > 1.0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "ambiguous_repeat_decoy_rate must be in [0.0, 1.0]",
        ));
    }
    Ok(())
}

fn validate_scaffold_polish_args(
    scaffold_groups: usize,
    contigs_per_scaffold: usize,
    contig_len: usize,
    read_len: usize,
    insert_size: usize,
    internal_pairs_per_contig: usize,
    true_link_pairs: usize,
    decoy_link_pairs: usize,
) -> io::Result<()> {
    if scaffold_groups == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "scaffold_groups must be greater than zero",
        ));
    }
    if contigs_per_scaffold < 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "contigs_per_scaffold must be at least two",
        ));
    }
    if contig_len <= read_len * 2 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "contig_len must be greater than 2 * read_len",
        ));
    }
    if read_len < DEFAULT_MINIMIZER_K {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("read_len must be at least {DEFAULT_MINIMIZER_K}"),
        ));
    }
    if insert_size < read_len * 2 || insert_size >= contig_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "insert_size must be at least 2 * read_len and smaller than contig_len",
        ));
    }
    if internal_pairs_per_contig == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "internal_pairs_per_contig must be greater than zero",
        ));
    }
    if true_link_pairs == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "true_link_pairs must be greater than zero",
        ));
    }
    if scaffold_groups < 2 && decoy_link_pairs > 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "decoy_link_pairs require at least two scaffold groups",
        ));
    }
    Ok(())
}

fn validate_error_correction_args(
    k: usize,
    trusted_roots: usize,
    weak_roots: usize,
    correctable_singletons_per_trusted: usize,
    protected_singletons_per_weak: usize,
) -> io::Result<()> {
    if k == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "k must be greater than zero",
        ));
    }
    if trusted_roots == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "trusted_roots must be greater than zero",
        ));
    }
    if weak_roots == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "weak_roots must be greater than zero",
        ));
    }
    if correctable_singletons_per_trusted == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "correctable_singletons_per_trusted must be greater than zero",
        ));
    }
    if protected_singletons_per_weak == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "protected_singletons_per_weak must be greater than zero",
        ));
    }
    Ok(())
}

fn mutate_contigs_for_polishing(
    correct_contigs: &[(String, Vec<u8>)],
    mutations_per_contig: usize,
    rng: &mut StdRng,
) -> io::Result<Vec<(String, Vec<u8>)>> {
    correct_contigs
        .iter()
        .map(|(header, sequence)| {
            if sequence.len() < 3 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "synthetic contigs must be at least 3 bp long",
                ));
            }
            let mut mutated = sequence.clone();
            let mutations = mutations_per_contig.max(1);
            for mutation_idx in 0..mutations {
                let pos = 1
                    + ((mutation_idx * 37 + rng.gen_range(0..sequence.len() - 2))
                        % (sequence.len() - 2));
                let original = mutated[pos];
                mutated[pos] = random_base_excluding(original, rng);
            }
            Ok((header.clone(), mutated))
        })
        .collect()
}

fn random_base_excluding(base: u8, rng: &mut StdRng) -> u8 {
    let mut replacement = base;
    while replacement == base {
        replacement = match rng.gen_range(0..4) {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        };
    }
    replacement
}

fn build_truth_scaffolds(
    contigs: &[(String, Vec<u8>)],
    scaffold_groups: usize,
    contigs_per_scaffold: usize,
    rng: &mut StdRng,
) -> Vec<ScaffoldPolishTruthScaffold> {
    let mut scaffolds = Vec::with_capacity(scaffold_groups);
    for group_idx in 0..scaffold_groups {
        let start = group_idx * contigs_per_scaffold;
        let mut placements = Vec::with_capacity(contigs_per_scaffold);
        for contig_offset in 0..contigs_per_scaffold {
            let contig_idx = start + contig_offset;
            placements.push(ScaffoldPolishPlacementSpec {
                contig: contigs[contig_idx].0.clone(),
                is_reverse: rng.gen_bool(0.5),
            });
        }
        scaffolds.push(ScaffoldPolishTruthScaffold {
            contigs: placements,
        });
    }
    scaffolds
}

#[allow(clippy::too_many_arguments)]
fn write_scaffold_polish_reads(
    correct_contigs: &[(String, Vec<u8>)],
    truth_scaffolds: &[ScaffoldPolishTruthScaffold],
    read_len: usize,
    insert_size: usize,
    internal_pairs_per_contig: usize,
    true_link_pairs: usize,
    decoy_link_pairs: usize,
    reads1_path: &Path,
    reads2_path: &Path,
    rng: &mut StdRng,
) -> io::Result<usize> {
    let mut writer1 = FastqWriter::try_new(path_to_str(reads1_path)?)?;
    let mut writer2 = FastqWriter::try_new(path_to_str(reads2_path)?)?;
    let contig_lookup = correct_contigs
        .iter()
        .enumerate()
        .map(|(idx, (header, sequence))| (header.clone(), (idx, sequence.as_slice())))
        .collect::<AHashMap<_, _>>();
    let mut pair_index = 0usize;

    for (header, sequence) in correct_contigs {
        let max_start = sequence.len().checked_sub(insert_size).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "insert_size exceeds contig_len",
            )
        })?;
        for offset in 0..internal_pairs_per_contig {
            let start1 = (offset * 17 + rng.gen_range(0..=max_start)) % (max_start + 1);
            let start2 = start1 + insert_size - read_len;
            let read1 = sequence[start1..start1 + read_len].to_vec();
            let read2 = reverse_complement(&sequence[start2..start2 + read_len]);
            write_pair(
                &mut writer1,
                &mut writer2,
                pair_index,
                &read1,
                &read2,
                "self",
                header,
                header,
            )?;
            pair_index += 1;
        }
    }

    for scaffold in truth_scaffolds {
        for adjacency in scaffold.contigs.windows(2) {
            let left = &adjacency[0];
            let right = &adjacency[1];
            let left_sequence = contig_lookup
                .get(&left.contig)
                .map(|(_, sequence)| *sequence)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("missing contig '{}' in truth scaffold", left.contig),
                    )
                })?;
            let right_sequence = contig_lookup
                .get(&right.contig)
                .map(|(_, sequence)| *sequence)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("missing contig '{}' in truth scaffold", right.contig),
                    )
                })?;
            let flank = read_len.min(24);
            for pair_offset in 0..true_link_pairs {
                let left_start = left_sequence.len() - read_len - (pair_offset % flank);
                let right_start = pair_offset % flank;
                let read1 =
                    oriented_read_slice(left_sequence, left_start, read_len, left.is_reverse);
                let read2 =
                    oriented_read_slice(right_sequence, right_start, read_len, right.is_reverse);
                write_pair(
                    &mut writer1,
                    &mut writer2,
                    pair_index,
                    &read1,
                    &read2,
                    "link",
                    &left.contig,
                    &right.contig,
                )?;
                pair_index += 1;
            }
        }
    }

    if decoy_link_pairs > 0 && truth_scaffolds.len() >= 2 {
        for pair_offset in 0..decoy_link_pairs {
            let left = &truth_scaffolds[0].contigs[pair_offset % truth_scaffolds[0].contigs.len()];
            let right = &truth_scaffolds[1].contigs[pair_offset % truth_scaffolds[1].contigs.len()];
            let left_sequence = contig_lookup
                .get(&left.contig)
                .map(|(_, sequence)| *sequence)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("missing decoy contig '{}'", left.contig),
                    )
                })?;
            let right_sequence = contig_lookup
                .get(&right.contig)
                .map(|(_, sequence)| *sequence)
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("missing decoy contig '{}'", right.contig),
                    )
                })?;
            let left_start = left_sequence.len() - read_len - (pair_offset % read_len.min(24));
            let right_start = pair_offset % read_len.min(24);
            let read1 = oriented_read_slice(left_sequence, left_start, read_len, left.is_reverse);
            let read2 =
                oriented_read_slice(right_sequence, right_start, read_len, right.is_reverse);
            write_pair(
                &mut writer1,
                &mut writer2,
                pair_index,
                &read1,
                &read2,
                "decoy",
                &left.contig,
                &right.contig,
            )?;
            pair_index += 1;
        }
    }

    Ok(pair_index)
}

fn oriented_read_slice(
    sequence: &[u8],
    start: usize,
    read_len: usize,
    is_reverse: bool,
) -> Vec<u8> {
    let slice = &sequence[start..start + read_len];
    if is_reverse {
        reverse_complement(slice)
    } else {
        slice.to_vec()
    }
}

fn write_pair(
    writer1: &mut FastqWriter,
    writer2: &mut FastqWriter,
    pair_index: usize,
    read1: &[u8],
    read2: &[u8],
    kind: &str,
    left_label: &str,
    right_label: &str,
) -> io::Result<()> {
    let read1_id = format!("{kind}_{pair_index:06}_{left_label}_{right_label}/1");
    let read2_id = format!("{kind}_{pair_index:06}_{left_label}_{right_label}/2");
    let quality = "I".repeat(read1.len());
    writer1.write_record(&FastqRecord {
        header: format!("@{read1_id}"),
        sequence: String::from_utf8(read1.to_vec()).map_err(invalid_utf8_error)?,
        plus: "+".to_string(),
        quality: quality.clone(),
    })?;
    writer2.write_record(&FastqRecord {
        header: format!("@{read2_id}"),
        sequence: String::from_utf8(read2.to_vec()).map_err(invalid_utf8_error)?,
        plus: "+".to_string(),
        quality,
    })?;
    Ok(())
}

fn decode_scaffold_output(
    scaffold_path: &Path,
    input_contigs: &[(String, Vec<u8>)],
) -> io::Result<Vec<Vec<ScaffoldPolishPlacementSpec>>> {
    let scaffolds = read_contigs(path_to_str(scaffold_path)?)?;
    let mut contig_lookup = AHashMap::new();
    for (header, sequence) in input_contigs {
        contig_lookup.insert(
            sequence.clone(),
            ScaffoldPolishPlacementSpec {
                contig: header.clone(),
                is_reverse: false,
            },
        );
        contig_lookup.insert(
            reverse_complement(sequence),
            ScaffoldPolishPlacementSpec {
                contig: header.clone(),
                is_reverse: true,
            },
        );
    }

    scaffolds
        .into_iter()
        .map(|(_header, sequence)| {
            split_scaffold_segments(&sequence)
                .into_iter()
                .map(|segment| {
                    contig_lookup.get(&segment).cloned().ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "observed scaffold segment did not match any input contig orientation",
                        )
                    })
                })
                .collect::<io::Result<Vec<_>>>()
        })
        .collect()
}

fn split_scaffold_segments(sequence: &[u8]) -> Vec<Vec<u8>> {
    let mut segments = Vec::new();
    let mut start = 0usize;
    while start < sequence.len() {
        while start < sequence.len() && matches!(sequence[start], b'N' | b'n') {
            start += 1;
        }
        if start >= sequence.len() {
            break;
        }
        let mut end = start;
        while end < sequence.len() && !matches!(sequence[end], b'N' | b'n') {
            end += 1;
        }
        segments.push(sequence[start..end].to_vec());
        start = end;
    }
    segments
}

fn scaffold_link_set_from_truth(
    scaffolds: &[ScaffoldPolishTruthScaffold],
) -> BTreeSet<(ScaffoldPolishPlacementSpec, ScaffoldPolishPlacementSpec)> {
    scaffolds
        .iter()
        .flat_map(|scaffold| scaffold.contigs.windows(2))
        .map(|adjacency| canonical_link(&adjacency[0], &adjacency[1]))
        .collect()
}

fn scaffold_link_set_from_observed(
    scaffolds: &[Vec<ScaffoldPolishPlacementSpec>],
) -> BTreeSet<(ScaffoldPolishPlacementSpec, ScaffoldPolishPlacementSpec)> {
    scaffolds
        .iter()
        .flat_map(|scaffold| scaffold.windows(2))
        .map(|adjacency| canonical_link(&adjacency[0], &adjacency[1]))
        .collect()
}

fn canonical_link(
    left: &ScaffoldPolishPlacementSpec,
    right: &ScaffoldPolishPlacementSpec,
) -> (ScaffoldPolishPlacementSpec, ScaffoldPolishPlacementSpec) {
    let forward = (left.clone(), right.clone());
    let reverse = (
        ScaffoldPolishPlacementSpec {
            contig: right.contig.clone(),
            is_reverse: !right.is_reverse,
        },
        ScaffoldPolishPlacementSpec {
            contig: left.contig.clone(),
            is_reverse: !left.is_reverse,
        },
    );
    if reverse < forward {
        reverse
    } else {
        forward
    }
}

fn canonical_truth_scaffold_set(
    scaffolds: &[ScaffoldPolishTruthScaffold],
) -> BTreeSet<Vec<ScaffoldPolishPlacementSpec>> {
    scaffolds
        .iter()
        .map(|scaffold| canonical_scaffold(&scaffold.contigs))
        .collect()
}

fn canonical_observed_scaffold_set(
    scaffolds: &[Vec<ScaffoldPolishPlacementSpec>],
) -> BTreeSet<Vec<ScaffoldPolishPlacementSpec>> {
    scaffolds
        .iter()
        .map(|scaffold| canonical_scaffold(scaffold))
        .collect()
}

fn canonical_scaffold(
    placements: &[ScaffoldPolishPlacementSpec],
) -> Vec<ScaffoldPolishPlacementSpec> {
    let forward = placements.to_vec();
    let reverse = placements
        .iter()
        .rev()
        .map(|placement| ScaffoldPolishPlacementSpec {
            contig: placement.contig.clone(),
            is_reverse: !placement.is_reverse,
        })
        .collect::<Vec<_>>();
    if reverse < forward {
        reverse
    } else {
        forward
    }
}

fn compare_polished_contigs(
    expected_contigs: &[ScaffoldPolishTruthContig],
    observed_contigs: &[(String, Vec<u8>)],
) -> (u64, u64, u64) {
    let observed = observed_contigs
        .iter()
        .map(|(header, sequence)| (header.as_str(), sequence.as_slice()))
        .collect::<AHashMap<_, _>>();
    let mut exact = 0u64;
    let mut bases_correct = 0u64;
    let mut bases_expected = 0u64;

    for truth in expected_contigs {
        let expected = truth.sequence.as_bytes();
        let observed_sequence = observed
            .get(truth.header.as_str())
            .copied()
            .unwrap_or_default();
        if observed_sequence == expected {
            exact += 1;
        }
        let overlap = expected.len().min(observed_sequence.len());
        bases_correct += expected[..overlap]
            .iter()
            .zip(&observed_sequence[..overlap])
            .filter(|(left, right)| left == right)
            .count() as u64;
        bases_expected += expected.len() as u64;
    }

    (exact, bases_correct, bases_expected)
}

fn f1_score(precision: f64, recall: f64) -> Option<f64> {
    if precision <= 0.0 || recall <= 0.0 {
        None
    } else {
        Some((2.0 * precision * recall) / (precision + recall))
    }
}

fn read_task_metadata(task_dir: &Path) -> io::Result<Option<ReadMappingTaskMetadata>> {
    let metadata_path = task_dir.join("task.json");
    if !metadata_path.exists() {
        return Ok(None);
    }

    let metadata =
        serde_json::from_reader(BufReader::new(File::open(&metadata_path)?)).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("parse {}: {}", metadata_path.display(), err),
            )
        })?;
    Ok(Some(metadata))
}

fn read_scaffold_polish_task_metadata(
    task_dir: &Path,
) -> io::Result<Option<ScaffoldPolishTaskMetadata>> {
    let metadata_path = task_dir.join("task.json");
    if !metadata_path.exists() {
        return Ok(None);
    }

    let metadata =
        serde_json::from_reader(BufReader::new(File::open(&metadata_path)?)).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("parse {}: {}", metadata_path.display(), err),
            )
        })?;
    Ok(Some(metadata))
}

fn read_error_correction_task_metadata(
    task_dir: &Path,
) -> io::Result<Option<ErrorCorrectionTaskMetadata>> {
    let metadata_path = task_dir.join("task.json");
    if !metadata_path.exists() {
        return Ok(None);
    }

    let metadata =
        serde_json::from_reader(BufReader::new(File::open(&metadata_path)?)).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("parse {}: {}", metadata_path.display(), err),
            )
        })?;
    Ok(Some(metadata))
}

fn read_contig_extraction_task_metadata(
    task_dir: &Path,
) -> io::Result<Option<ContigExtractionTaskMetadata>> {
    let metadata_path = task_dir.join("task.json");
    if !metadata_path.exists() {
        return Ok(None);
    }

    let metadata =
        serde_json::from_reader(BufReader::new(File::open(&metadata_path)?)).map_err(|err| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("parse {}: {}", metadata_path.display(), err),
            )
        })?;
    Ok(Some(metadata))
}

fn resolve_contigs_path(task_dir: &Path) -> io::Result<PathBuf> {
    resolve_required_path(task_dir, &["contigs.fa", "contigs.fasta", "contigs.fa.gz"])
}

fn resolve_reads1_path(task_dir: &Path) -> io::Result<PathBuf> {
    resolve_required_path(
        task_dir,
        &[
            "reads_1.fastq",
            "reads_1.fastq.gz",
            "reads_1.fq",
            "reads_1.fq.gz",
        ],
    )
}

fn resolve_reads2_path(task_dir: &Path) -> Option<PathBuf> {
    resolve_optional_path(
        task_dir,
        &[
            "reads_2.fastq",
            "reads_2.fastq.gz",
            "reads_2.fq",
            "reads_2.fq.gz",
        ],
    )
}

fn resolve_truth_path(task_dir: &Path) -> io::Result<PathBuf> {
    resolve_required_path(task_dir, &["truth.tsv"])
}

fn resolve_scaffold_polish_truth_path(task_dir: &Path) -> io::Result<PathBuf> {
    resolve_required_path(task_dir, &["truth.json"])
}

fn resolve_error_correction_counts_path(task_dir: &Path) -> io::Result<PathBuf> {
    resolve_required_path(task_dir, &["counts.json"])
}

fn resolve_error_correction_truth_path(task_dir: &Path) -> io::Result<PathBuf> {
    resolve_required_path(task_dir, &["truth.json"])
}

fn resolve_contig_extraction_counts_path(task_dir: &Path) -> io::Result<PathBuf> {
    resolve_required_path(task_dir, &["counts.json"])
}

fn resolve_contig_extraction_graph_path(task_dir: &Path) -> io::Result<PathBuf> {
    resolve_required_path(task_dir, &["graph.json"])
}

fn resolve_contig_extraction_truth_path(task_dir: &Path) -> io::Result<PathBuf> {
    resolve_required_path(task_dir, &["truth.json"])
}

fn load_scaffold_polish_truth_file(task_dir: &Path) -> io::Result<ScaffoldPolishTruthFile> {
    let truth_path = resolve_scaffold_polish_truth_path(task_dir)?;
    serde_json::from_reader(BufReader::new(File::open(&truth_path)?)).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("parse {}: {}", truth_path.display(), err),
        )
    })
}

fn load_error_correction_task_file(task_dir: &Path) -> io::Result<ErrorCorrectionTaskFile> {
    let task_path = resolve_error_correction_counts_path(task_dir)?;
    serde_json::from_reader(BufReader::new(File::open(&task_path)?)).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("parse {}: {}", task_path.display(), err),
        )
    })
}

fn load_error_correction_truth_file(task_dir: &Path) -> io::Result<ErrorCorrectionTruthFile> {
    let truth_path = resolve_error_correction_truth_path(task_dir)?;
    serde_json::from_reader(BufReader::new(File::open(&truth_path)?)).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("parse {}: {}", truth_path.display(), err),
        )
    })
}

fn load_contig_extraction_counts_file(
    task_dir: &Path,
) -> io::Result<Vec<ContigExtractionCountSpec>> {
    let counts_path = resolve_contig_extraction_counts_path(task_dir)?;
    serde_json::from_reader(BufReader::new(File::open(&counts_path)?)).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("parse {}: {}", counts_path.display(), err),
        )
    })
}

fn load_contig_extraction_graph_file(task_dir: &Path) -> io::Result<ContigExtractionGraphFile> {
    let graph_path = resolve_contig_extraction_graph_path(task_dir)?;
    serde_json::from_reader(BufReader::new(File::open(&graph_path)?)).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("parse {}: {}", graph_path.display(), err),
        )
    })
}

fn load_contig_extraction_truth_file(task_dir: &Path) -> io::Result<ContigExtractionTruthFile> {
    let truth_path = resolve_contig_extraction_truth_path(task_dir)?;
    serde_json::from_reader(BufReader::new(File::open(&truth_path)?)).map_err(|err| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("parse {}: {}", truth_path.display(), err),
        )
    })
}

fn resolve_required_path(task_dir: &Path, candidates: &[&str]) -> io::Result<PathBuf> {
    resolve_optional_path(task_dir, candidates).ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::NotFound,
            format!(
                "{} is missing one of [{}]",
                task_dir.display(),
                candidates.join(", ")
            ),
        )
    })
}

fn resolve_optional_path(task_dir: &Path, candidates: &[&str]) -> Option<PathBuf> {
    candidates
        .iter()
        .map(|candidate| task_dir.join(candidate))
        .find(|path| path.exists())
}

fn normalize_read_id(header: &str) -> Option<&str> {
    header
        .trim()
        .trim_start_matches('@')
        .split_whitespace()
        .next()
        .filter(|value| !value.is_empty())
}

fn ratio(numerator: u64, denominator: u64) -> Option<f64> {
    if denominator == 0 {
        None
    } else {
        Some(numerator as f64 / denominator as f64)
    }
}

fn ratio_f64(numerator: f64, denominator: f64) -> Option<f64> {
    if denominator == 0.0 {
        None
    } else {
        Some(numerator / denominator)
    }
}

fn mapping_score(
    primary_near_rate: Option<f64>,
    primary_exact_rate: Option<f64>,
    scaffold_exact_rate: Option<f64>,
    primary_mapped_rate: Option<f64>,
    primary_specificity_rate: Option<f64>,
    scaffold_specificity_rate: Option<f64>,
    pair_cross_contig_rate: Option<f64>,
    pair_specificity_rate: Option<f64>,
) -> f64 {
    let mut weighted_sum = 0.0;
    let mut total_weight = 0.0;

    for (metric, weight) in [
        (primary_near_rate, 0.30),
        (primary_exact_rate, 0.15),
        (scaffold_exact_rate, 0.10),
        (primary_mapped_rate, 0.05),
        (primary_specificity_rate, 0.10),
        (scaffold_specificity_rate, 0.10),
        (pair_cross_contig_rate, 0.10),
        (pair_specificity_rate, 0.10),
    ] {
        if let Some(metric) = metric {
            weighted_sum += metric * weight;
            total_weight += weight;
        }
    }

    if total_weight == 0.0 {
        0.0
    } else {
        weighted_sum / total_weight
    }
}

fn error_correction_score(
    correction_f1: Option<f64>,
    preserved_retention_rate: Option<f64>,
    trusted_target_exact_rate: Option<f64>,
) -> f64 {
    let mut weighted_sum = 0.0;
    let mut total_weight = 0.0;

    for (metric, weight) in [
        (correction_f1, 0.50),
        (preserved_retention_rate, 0.25),
        (trusted_target_exact_rate, 0.25),
    ] {
        if let Some(metric) = metric {
            weighted_sum += metric * weight;
            total_weight += weight;
        }
    }

    if total_weight == 0.0 {
        0.0
    } else {
        weighted_sum / total_weight
    }
}

fn contig_extraction_score(
    truth_kmer_f1: Option<f64>,
    exact_contig_rate: Option<f64>,
    contig_count_agreement: Option<f64>,
) -> f64 {
    let mut weighted_sum = 0.0;
    let mut total_weight = 0.0;

    for (metric, weight) in [
        (truth_kmer_f1, 0.50),
        (exact_contig_rate, 0.35),
        (contig_count_agreement, 0.15),
    ] {
        if let Some(metric) = metric {
            weighted_sum += metric * weight;
            total_weight += weight;
        }
    }

    if total_weight == 0.0 {
        0.0
    } else {
        weighted_sum / total_weight
    }
}

fn print_read_mapping_report(report: &ReadMappingBenchmarkReport) {
    println!("Read-mapping component benchmark");
    println!(
        "  Tasks: {} | k/w: {}/{} | min matches: {}/{} | tolerance: {} bp",
        report.summary.task_count,
        report.minimizer_k,
        report.minimizer_w,
        report.min_primary_matches,
        report.min_scaffold_matches,
        report.position_tolerance
    );
    println!(
        "  Score: {:.4} | Reads: {} | Throughput: {:.1} reads/s",
        report.summary.score, report.summary.reads_evaluated, report.summary.reads_per_second
    );
    println!(
        "  Primary exact/near/mapped: {} / {} / {}",
        format_rate(report.summary.primary_exact_rate),
        format_rate(report.summary.primary_near_rate),
        format_rate(report.summary.primary_mapped_rate)
    );
    println!(
        "  Primary specificity: {} | Scaffold specificity: {} | Pair specificity: {} | Pair contig: {} | Scaffold exact: {} | Position MAE: {}",
        format_rate(report.summary.primary_specificity_rate),
        format_rate(report.summary.scaffold_specificity_rate),
        format_rate(report.summary.pair_specificity_rate),
        format_rate(report.summary.pair_cross_contig_rate),
        format_rate(report.summary.scaffold_exact_rate),
        format_metric(report.summary.primary_position_mae, "bp")
    );

    if report.tasks.len() > 1 {
        println!("Per task:");
        for task in &report.tasks {
            println!(
                "  {} | score {:.4} | exact {} | near {} | primary spec {} | scaffold spec {} | pair spec {} | pair contig {} | scaffold {} | {:.1} reads/s",
                task.task_name,
                task.score,
                format_rate(task.primary_exact_rate),
                format_rate(task.primary_near_rate),
                format_rate(task.primary_specificity_rate),
                format_rate(task.scaffold_specificity_rate),
                format_rate(task.pair_specificity_rate),
                format_rate(task.pair_cross_contig_rate),
                format_rate(task.scaffold_exact_rate),
                task.reads_per_second
            );
        }
    }
}

fn print_contig_extraction_report(report: &ContigExtractionBenchmarkReport) {
    println!("Component: {}", report.component);
    println!(
        "Config: prefer_high_count_seeds={} prefer_non_repeat_seeds={} enable_repeat_seed_completion={} suppress_redundant_contigs={}",
        report.prefer_high_count_seeds,
        report.prefer_non_repeat_seeds,
        report.enable_repeat_seed_completion,
        report.suppress_redundant_contigs
    );
    println!(
        "Summary: score {:.4} | exact contigs {} | truth k-mer F1 {} | count agreement {} | {:.1} nodes/s",
        report.summary.score,
        format_rate(report.summary.exact_contig_rate),
        format_rate(report.summary.truth_kmer_f1),
        format_rate(report.summary.contig_count_agreement),
        report.summary.graph_nodes_per_second
    );
    println!(
        "  Exact contigs: {}/{} | truth k-mers correct/observed/expected: {}/{}/{}",
        report.summary.exact_contigs,
        report.summary.expected_contigs,
        report.summary.truth_kmers_correct,
        report.summary.truth_kmers_observed,
        report.summary.truth_kmers_expected
    );

    if report.tasks.len() > 1 {
        println!("Per task:");
        for task in &report.tasks {
            println!(
                "  {} [{}] | score {:.4} | exact {} | k-mer F1 {} | count {} | {:.1} nodes/s",
                task.task_name,
                task.profile,
                task.score,
                format_rate(task.exact_contig_rate),
                format_rate(task.truth_kmer_f1),
                format_rate(task.contig_count_agreement),
                task.graph_nodes_per_second
            );
        }
    }
}

fn print_branch_resolution_report(report: &BranchResolutionBenchmarkReport) {
    println!("Branch-resolution component benchmark");
    println!(
        "  Tasks: {} | Cases: {} | support floor/margin: {}/{} | prefer non-repeat: {} | Throughput: {:.1} cases/s",
        report.summary.task_count,
        report.summary.cases_evaluated,
        report.branch_support_min_win,
        report.branch_support_min_margin,
        report.prefer_non_repeat_branches,
        report.summary.cases_per_second
    );
    println!(
        "  Score: {:.4} | Exact: {} | Any choice: {}",
        report.summary.score,
        format_rate(report.summary.exact_rate),
        format_rate(report.summary.choice_rate)
    );
    if !report.summary.scenario_metrics.is_empty() {
        let scenario_summary = report
            .summary
            .scenario_metrics
            .iter()
            .map(|(scenario, metrics)| format!("{scenario} {}", format_rate(metrics.exact_rate)))
            .collect::<Vec<_>>()
            .join(" | ");
        println!("  Scenario exact: {}", scenario_summary);
    }

    if report.tasks.len() > 1 {
        println!("Per task:");
        for task in &report.tasks {
            println!(
                "  {} | score {:.4} | exact {} | choice {} | {:.1} cases/s",
                task.task_name,
                task.score,
                format_rate(task.exact_rate),
                format_rate(task.choice_rate),
                task.cases_per_second
            );
        }
    }
}

fn print_scaffold_polish_report(report: &ScaffoldPolishBenchmarkReport) {
    println!("Component: {}", report.component);
    println!(
        "Config: min_scaffold_links={} min_primary_matches={} min_scaffold_matches={}",
        report.min_scaffold_links, report.min_primary_matches, report.min_scaffold_matches
    );
    println!(
        "Summary: score {:.4} | scaffold F1 {} | polished exact {} | polished bases {} | {:.1} pairs/s",
        report.summary.score,
        format_rate(report.summary.scaffold_link_f1),
        format_rate(report.summary.polished_exact_rate),
        format_rate(report.summary.polished_base_accuracy),
        report.summary.pairs_per_second
    );
    println!(
        "  Links correct/observed/expected: {}/{}/{}",
        report.summary.correct_links, report.summary.observed_links, report.summary.expected_links
    );
    println!(
        "  Exact scaffolds: {}/{} | polished exact contigs: {}/{}",
        report.summary.exact_scaffolds,
        report.summary.expected_scaffolds,
        report.summary.polished_contigs_exact,
        report.summary.polished_contigs_expected
    );

    if report.tasks.len() > 1 {
        println!("Per task:");
        for task in &report.tasks {
            println!(
                "  {} | score {:.4} | scaffold F1 {} | polished exact {} | {:.1} pairs/s",
                task.task_name,
                task.score,
                format_rate(task.scaffold_link_f1),
                format_rate(task.polished_exact_rate),
                task.pairs_per_second
            );
        }
    }
}

fn print_error_correction_report(report: &ErrorCorrectionBenchmarkReport) {
    println!("Component: {}", report.component);
    println!(
        "Config: min_count={}, min_trusted_count={}",
        report.min_count, report.min_trusted_count
    );
    println!(
        "Summary: score {:.4} | correction F1 {} | preserved {} | trusted exact {} | {:.1} kmers/s",
        report.summary.score,
        format_rate(report.summary.correction_f1),
        format_rate(report.summary.preserved_retention_rate),
        format_rate(report.summary.trusted_target_exact_rate),
        report.summary.kmers_per_second
    );
    println!(
        "  Correctly removed/observed/expected: {}/{}/{}",
        report.summary.correctly_removed_singletons,
        report.summary.observed_removed_singletons,
        report.summary.expected_removed_singletons
    );
    println!(
        "  Preserved retained: {}/{} | trusted targets exact: {}/{}",
        report.summary.preserved_singletons_retained,
        report.summary.expected_preserved_singletons,
        report.summary.trusted_targets_exact,
        report.summary.trusted_targets
    );

    if report.tasks.len() > 1 {
        println!("Per task:");
        for task in &report.tasks {
            println!(
                "  {} | score {:.4} | F1 {} | preserved {} | trusted {} | {:.1} kmers/s",
                task.task_name,
                task.score,
                format_rate(task.correction_f1),
                format_rate(task.preserved_retention_rate),
                format_rate(task.trusted_target_exact_rate),
                task.kmers_per_second
            );
        }
    }
}

fn format_rate(rate: Option<f64>) -> String {
    match rate {
        Some(rate) => format!("{:.2}%", rate * 100.0),
        None => "n/a".to_string(),
    }
}

fn format_metric(value: Option<f64>, suffix: &str) -> String {
    match value {
        Some(value) => format!("{value:.2} {suffix}"),
        None => "n/a".to_string(),
    }
}

fn random_dna(len: usize, rng: &mut StdRng) -> Vec<u8> {
    let alphabet = [b'A', b'C', b'G', b'T'];
    (0..len)
        .map(|_| alphabet[rng.gen_range(0..alphabet.len())])
        .collect()
}

fn inject_errors(sequence: &mut [u8], error_rate: f64, rng: &mut StdRng) {
    if error_rate == 0.0 {
        return;
    }

    for base in sequence {
        if rng.gen_bool(error_rate) {
            *base = random_alt_base(*base, rng);
        }
    }
}

fn random_alt_base(base: u8, rng: &mut StdRng) -> u8 {
    let alphabet = [b'A', b'C', b'G', b'T'];
    loop {
        let candidate = alphabet[rng.gen_range(0..alphabet.len())];
        if candidate != base {
            return candidate;
        }
    }
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

fn find_subsequence(sequence: &[u8], needle: &[u8]) -> Option<usize> {
    sequence
        .windows(needle.len())
        .position(|window| window == needle)
}

fn invalid_utf8_error(err: std::string::FromUtf8Error) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("generated DNA sequence is not valid UTF-8: {err}"),
    )
}

fn path_to_str(path: &Path) -> io::Result<&str> {
    path.to_str().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("path '{}' is not valid UTF-8", path.display()),
        )
    })
}

#[cfg(test)]
mod tests {
    use super::{
        evaluate_branch_resolution_bench, evaluate_contig_extraction_bench,
        evaluate_error_correction_bench, evaluate_read_mapping_bench,
        evaluate_scaffold_polish_bench, prepare_branch_resolution_panel,
        prepare_branch_resolution_task, prepare_contig_extraction_panel,
        prepare_contig_extraction_task, prepare_contig_extraction_task_with_profile,
        prepare_error_correction_task, prepare_read_mapping_task, prepare_scaffold_polish_task,
    };
    use std::collections::BTreeMap;
    use tempfile::TempDir;

    #[test]
    fn generated_single_task_round_trips_through_benchmark() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir.path().join("mapper_task");

        prepare_read_mapping_task(
            task_dir.to_str().expect("utf8 path"),
            true,
            4,
            1200,
            75,
            32,
            220,
            45,
            0.0,
            0.0,
            0.0,
            7,
            11,
            5,
        )
        .expect("prepare task");

        let report = evaluate_read_mapping_bench(&task_dir, 11, 5, 3, 2, 8).expect("evaluate task");
        assert_eq!(report.tasks.len(), 1);
        let task = &report.tasks[0];
        assert_eq!(task.reads_evaluated, 64);
        assert_eq!(task.primary_exact_rate, Some(1.0));
        assert_eq!(task.scaffold_exact_rate, Some(1.0));
        assert_eq!(task.pair_expected, 32);
        assert_eq!(task.pair_mapped_rate, Some(1.0));
        assert_eq!(task.pair_cross_contig_rate, Some(1.0));
        assert!(task.reads_per_second > 0.0);
    }

    #[test]
    fn decoy_reads_contribute_specificity_counts() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir.path().join("mapper_task_with_decoys");

        prepare_read_mapping_task(
            task_dir.to_str().expect("utf8 path"),
            false,
            4,
            1200,
            75,
            32,
            220,
            45,
            0.0,
            0.25,
            0.0,
            13,
            11,
            5,
        )
        .expect("prepare task");

        let report = evaluate_read_mapping_bench(&task_dir, 11, 5, 3, 2, 8).expect("evaluate task");
        let task = &report.tasks[0];
        assert_eq!(task.reads_evaluated, 40);
        assert_eq!(task.primary_unexpected_total, 8);
        assert_eq!(task.primary_specificity_rate, Some(1.0));
    }

    #[test]
    fn ambiguous_repeat_decoys_penalize_scaffold_specificity() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir.path().join("mapper_task_with_repeat_decoys");

        prepare_read_mapping_task(
            task_dir.to_str().expect("utf8 path"),
            false,
            5,
            1400,
            75,
            32,
            220,
            90,
            0.0,
            0.0,
            0.25,
            29,
            11,
            5,
        )
        .expect("prepare task");

        let report = evaluate_read_mapping_bench(&task_dir, 9, 4, 2, 1, 8).expect("evaluate task");
        let task = &report.tasks[0];
        assert_eq!(task.reads_evaluated, 40);
        assert_eq!(task.scaffold_unexpected_total, 8);
        assert!(task.scaffold_unexpected_mapped > 0);
        assert!(task.scaffold_specificity_rate.unwrap_or(1.0) < 1.0);
    }

    #[test]
    fn paired_repeat_decoys_penalize_pair_specificity() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir
            .path()
            .join("paired_mapper_task_with_repeat_decoys");

        prepare_read_mapping_task(
            task_dir.to_str().expect("utf8 path"),
            true,
            5,
            1400,
            75,
            32,
            220,
            90,
            0.0,
            0.0,
            0.25,
            31,
            11,
            5,
        )
        .expect("prepare task");

        let report = evaluate_read_mapping_bench(&task_dir, 9, 4, 2, 1, 8).expect("evaluate task");
        let task = &report.tasks[0];
        assert_eq!(task.pair_unexpected_total, 8);
        assert!(task.pair_unexpected_mapped > 0);
        assert!(task.pair_specificity_rate.unwrap_or(1.0) < 1.0);
        assert!(task.pair_cross_contig_errors > 0);
    }

    #[test]
    fn task_root_discovers_multiple_child_tasks() {
        let temp_dir = TempDir::new().expect("temp dir");
        for task_name in ["task_a", "task_b"] {
            let task_dir = temp_dir.path().join(task_name);
            prepare_read_mapping_task(
                task_dir.to_str().expect("utf8 path"),
                false,
                3,
                900,
                60,
                24,
                160,
                35,
                0.0,
                0.0,
                0.0,
                11,
                11,
                5,
            )
            .expect("prepare task");
        }

        let report =
            evaluate_read_mapping_bench(temp_dir.path(), 11, 5, 3, 2, 8).expect("evaluate root");
        assert_eq!(report.tasks.len(), 2);
        assert_eq!(report.summary.task_count, 2);
        assert_eq!(report.summary.reads_evaluated, 48);
    }

    #[test]
    fn prepare_panel_creates_discoverable_task_root() {
        let temp_dir = TempDir::new().expect("temp dir");
        let panel_dir = temp_dir.path().join("panel");

        super::prepare_read_mapping_panel(
            panel_dir.to_str().expect("utf8 path"),
            3,
            false,
            3,
            900,
            60,
            12,
            160,
            35,
            0.0,
            0.0,
            0.0,
            17,
            1,
            11,
            5,
        )
        .expect("prepare panel");

        let report =
            evaluate_read_mapping_bench(&panel_dir, 11, 5, 3, 2, 8).expect("evaluate panel");
        assert_eq!(report.tasks.len(), 3);
        assert_eq!(report.summary.task_count, 3);
        assert_eq!(report.summary.reads_evaluated, 36);
    }

    #[test]
    fn support_stress_branch_profile_biases_margin_and_floor_cases() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir.path().join("branch_task_support_stress");

        super::prepare_branch_resolution_task_with_profile(
            task_dir.to_str().expect("utf8 path"),
            24,
            19,
            "support_stress",
        )
        .expect("prepare branch task");

        let task_file = super::load_branch_resolution_task_file(&task_dir).expect("task file");
        let mut scenario_counts = BTreeMap::new();
        for case in task_file.cases {
            *scenario_counts.entry(case.scenario).or_insert(0usize) += 1;
        }

        assert_eq!(scenario_counts.get("support_margin_gate"), Some(&6));
        assert_eq!(scenario_counts.get("support_floor_gate"), Some(&6));
        assert_eq!(scenario_counts.get("coverage_closeness"), Some(&6));
        assert_eq!(scenario_counts.get("repeat_tiebreak"), Some(&3));
        assert_eq!(scenario_counts.get("non_repeat_preferred"), Some(&3));
        assert!(!scenario_counts.contains_key("support_dominates"));
    }

    #[test]
    fn generated_branch_task_round_trips_through_benchmark() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir.path().join("branch_task");

        prepare_branch_resolution_task(task_dir.to_str().expect("utf8 path"), 24, 19)
            .expect("prepare branch task");

        let report =
            evaluate_branch_resolution_bench(&task_dir, 1, 2, true).expect("evaluate branch task");
        assert_eq!(report.tasks.len(), 1);
        let task = &report.tasks[0];
        assert_eq!(task.cases_evaluated, 24);
        assert_eq!(task.choice_rate, Some(1.0));
        assert!(task.score > 0.85);
    }

    #[test]
    fn branch_task_penalizes_overeager_support_defaults() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir.path().join("branch_task");

        prepare_branch_resolution_task(task_dir.to_str().expect("utf8 path"), 24, 19)
            .expect("prepare branch task");

        let conservative =
            evaluate_branch_resolution_bench(&task_dir, 1, 2, true).expect("evaluate branch task");
        let eager =
            evaluate_branch_resolution_bench(&task_dir, 1, 1, true).expect("evaluate branch task");
        assert!(eager.summary.score < conservative.summary.score);
    }

    #[test]
    fn prepare_branch_panel_creates_discoverable_task_root() {
        let temp_dir = TempDir::new().expect("temp dir");
        let panel_dir = temp_dir.path().join("branch_panel");

        prepare_branch_resolution_panel(panel_dir.to_str().expect("utf8 path"), 3, 12, 23, 1)
            .expect("prepare branch panel");

        let report = evaluate_branch_resolution_bench(&panel_dir, 1, 2, true)
            .expect("evaluate branch panel");
        assert_eq!(report.tasks.len(), 3);
        assert_eq!(report.summary.task_count, 3);
        assert_eq!(report.summary.cases_evaluated, 36);
        assert_eq!(report.summary.choice_rate, Some(1.0));
        assert!(report.summary.score > 0.85);
        assert!(report
            .summary
            .scenario_metrics
            .contains_key("support_margin_gate"));
    }

    #[test]
    fn generated_scaffold_polish_task_round_trips_through_benchmark() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir.path().join("scaffold_polish_task");

        prepare_scaffold_polish_task(
            task_dir.to_str().expect("utf8 path"),
            2,
            3,
            400,
            80,
            200,
            12,
            5,
            2,
            1,
            7,
        )
        .expect("prepare scaffold+polish task");

        let report = evaluate_scaffold_polish_bench(&task_dir, 3, 2, 4)
            .expect("evaluate scaffold+polish task");
        assert_eq!(report.tasks.len(), 1);
        let task = &report.tasks[0];
        assert_eq!(task.expected_scaffolds, 2);
        assert_eq!(task.expected_links, 4);
        assert_eq!(task.scaffold_link_f1, Some(1.0));
        assert_eq!(task.polished_exact_rate, Some(1.0));
        assert!(task.score > 0.99);
    }

    #[test]
    fn scaffold_polish_task_penalizes_weak_link_thresholds() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir.path().join("scaffold_polish_task");

        prepare_scaffold_polish_task(
            task_dir.to_str().expect("utf8 path"),
            2,
            3,
            400,
            80,
            200,
            12,
            5,
            2,
            1,
            9,
        )
        .expect("prepare scaffold+polish task");

        let conservative = evaluate_scaffold_polish_bench(&task_dir, 3, 2, 4)
            .expect("evaluate scaffold+polish task");
        let eager = evaluate_scaffold_polish_bench(&task_dir, 1, 2, 4)
            .expect("evaluate scaffold+polish task");
        assert!(eager.summary.score < conservative.summary.score);
        assert!(
            eager.summary.scaffold_link_precision.unwrap_or(1.0)
                < conservative.summary.scaffold_link_precision.unwrap_or(0.0)
        );
    }

    #[test]
    fn generated_error_correction_task_round_trips_through_benchmark() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir.path().join("error_correction_task");

        prepare_error_correction_task(task_dir.to_str().expect("utf8 path"), 21, 16, 8, 2, 1, 7)
            .expect("prepare error-correction task");

        let report = evaluate_error_correction_bench(&task_dir, 1, 4)
            .expect("evaluate error-correction task");
        assert_eq!(report.tasks.len(), 1);
        let task = &report.tasks[0];
        assert!(task.expected_removed_singletons > 0);
        assert_eq!(task.preserved_retention_rate, Some(1.0));
        assert!(task.correction_recall.unwrap_or(0.0) > 0.95);
        assert!(task.score > 0.95);
    }

    #[test]
    fn error_correction_task_penalizes_overly_strict_min_count() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir.path().join("error_correction_task");

        prepare_error_correction_task(task_dir.to_str().expect("utf8 path"), 21, 24, 12, 2, 1, 11)
            .expect("prepare error-correction task");

        let permissive =
            evaluate_error_correction_bench(&task_dir, 1, 4).expect("evaluate permissive config");
        let strict =
            evaluate_error_correction_bench(&task_dir, 3, 4).expect("evaluate strict config");

        assert!(strict.summary.score < permissive.summary.score);
        assert!(
            strict.summary.correction_recall.unwrap_or(0.0)
                < permissive.summary.correction_recall.unwrap_or(0.0)
        );
    }

    #[test]
    fn error_correction_task_penalizes_overly_lenient_trusted_floor() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir.path().join("error_correction_task");

        prepare_error_correction_task(task_dir.to_str().expect("utf8 path"), 21, 24, 12, 2, 1, 17)
            .expect("prepare error-correction task");

        let baseline =
            evaluate_error_correction_bench(&task_dir, 1, 4).expect("evaluate baseline config");
        let overly_lenient = evaluate_error_correction_bench(&task_dir, 1, 3)
            .expect("evaluate overly lenient config");

        assert!(overly_lenient.summary.score < baseline.summary.score);
        assert!(
            overly_lenient
                .summary
                .preserved_retention_rate
                .unwrap_or(0.0)
                < baseline.summary.preserved_retention_rate.unwrap_or(0.0)
        );
    }

    #[test]
    fn generated_contig_extraction_task_round_trips_through_benchmark() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir.path().join("contig_extraction_task");

        prepare_contig_extraction_task(task_dir.to_str().expect("utf8 path"), 11, 4, 4, 3, 7)
            .expect("prepare contig-extraction task");

        let report = evaluate_contig_extraction_bench(&task_dir, true, true, true, false)
            .expect("evaluate task");
        assert_eq!(report.tasks.len(), 1);
        let task = &report.tasks[0];
        assert_eq!(task.expected_contigs, 4);
        assert_eq!(task.exact_contig_rate, Some(1.0));
        assert!(task.truth_kmer_f1.unwrap_or(0.0) > 0.85);
        assert!(task.score > 0.80);
    }

    #[test]
    fn contig_extraction_task_penalizes_low_count_seed_order() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir.path().join("contig_extraction_task");

        prepare_contig_extraction_task(task_dir.to_str().expect("utf8 path"), 11, 4, 4, 3, 11)
            .expect("prepare contig-extraction task");

        let baseline =
            evaluate_contig_extraction_bench(&task_dir, true, true, true, false).expect("baseline");
        let low_count_first =
            evaluate_contig_extraction_bench(&task_dir, false, true, true, false).expect("variant");

        assert!(low_count_first.summary.score < baseline.summary.score);
        assert!(
            low_count_first.summary.exact_contig_rate.unwrap_or(0.0)
                < baseline.summary.exact_contig_rate.unwrap_or(0.0)
        );
    }

    #[test]
    fn contig_extraction_task_penalizes_missing_repeat_completion() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir.path().join("contig_extraction_repeat_task");

        prepare_contig_extraction_task_with_profile(
            task_dir.to_str().expect("utf8 path"),
            11,
            4,
            4,
            3,
            13,
            "repeat_completion",
        )
        .expect("prepare contig-extraction repeat task");

        let baseline =
            evaluate_contig_extraction_bench(&task_dir, true, true, true, false).expect("baseline");
        let no_repeat_completion =
            evaluate_contig_extraction_bench(&task_dir, true, true, false, false).expect("variant");

        assert!(no_repeat_completion.summary.score < baseline.summary.score);
        assert!(
            no_repeat_completion
                .summary
                .exact_contig_rate
                .unwrap_or(0.0)
                < baseline.summary.exact_contig_rate.unwrap_or(0.0)
        );
    }

    #[test]
    fn contig_extraction_task_penalizes_repeat_first_seeding_under_suppression() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir
            .path()
            .join("contig_extraction_repeat_priority_task");

        prepare_contig_extraction_task_with_profile(
            task_dir.to_str().expect("utf8 path"),
            11,
            4,
            4,
            3,
            17,
            "repeat_priority",
        )
        .expect("prepare contig-extraction repeat-priority task");

        let baseline =
            evaluate_contig_extraction_bench(&task_dir, true, true, true, true).expect("baseline");
        let repeat_first_variant =
            evaluate_contig_extraction_bench(&task_dir, true, false, true, true)
                .expect("repeat-first variant");

        assert_eq!(baseline.summary.exact_contig_rate, Some(1.0));
        assert!(repeat_first_variant.summary.score < baseline.summary.score);
        assert!(
            repeat_first_variant
                .summary
                .exact_contig_rate
                .unwrap_or(0.0)
                < baseline.summary.exact_contig_rate.unwrap_or(0.0)
        );
        assert!(
            repeat_first_variant.summary.truth_kmer_f1.unwrap_or(0.0)
                < baseline.summary.truth_kmer_f1.unwrap_or(0.0)
        );
    }

    #[test]
    fn prepare_contig_extraction_panel_creates_discoverable_task_root() {
        let temp_dir = TempDir::new().expect("temp dir");
        let panel_dir = temp_dir.path().join("contig_panel");

        prepare_contig_extraction_panel(
            panel_dir.to_str().expect("utf8 path"),
            4,
            11,
            4,
            4,
            3,
            17,
            1,
        )
        .expect("prepare contig panel");

        let report = evaluate_contig_extraction_bench(&panel_dir, true, true, true, false)
            .expect("evaluate");
        assert_eq!(report.tasks.len(), 4);
        assert_eq!(report.summary.task_count, 4);
        assert!(report.summary.score > 0.70);
    }

    #[test]
    fn contig_extraction_task_rewards_redundant_contig_suppression() {
        let temp_dir = TempDir::new().expect("temp dir");
        let task_dir = temp_dir.path().join("contig_extraction_branch_task");

        prepare_contig_extraction_task_with_profile(
            task_dir.to_str().expect("utf8 path"),
            11,
            4,
            4,
            3,
            19,
            "branching",
        )
        .expect("prepare contig-extraction branching task");

        let baseline =
            evaluate_contig_extraction_bench(&task_dir, true, true, true, false).expect("baseline");
        let suppressed = evaluate_contig_extraction_bench(&task_dir, true, true, true, true)
            .expect("suppressed");

        assert!(suppressed.summary.score > baseline.summary.score);
        assert!(
            suppressed.summary.contig_count_agreement.unwrap_or(0.0)
                > baseline.summary.contig_count_agreement.unwrap_or(0.0)
        );
        assert!(
            suppressed.summary.truth_kmer_precision.unwrap_or(0.0)
                >= baseline.summary.truth_kmer_precision.unwrap_or(0.0)
        );
        assert_eq!(
            suppressed.summary.exact_contig_rate,
            baseline.summary.exact_contig_rate
        );
    }
}
