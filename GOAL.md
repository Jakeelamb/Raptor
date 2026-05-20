<goal>
Build Raptor into a full-fledged Rust replacement for Trinity RNA-seq de novo transcriptome assembly. Do not stop at feature resemblance: Raptor is complete only when every major Trinity pipeline stage has a mapped Raptor equivalent, reproducible benchmarks, and comparable biological outputs across a frozen benchmark panel.
</goal>

<context>
Start by reading:
- `SPEC.md`
- `README.md`
- `docs/rescue-baseline.md`
- `src/cli_main.rs`
- `src/main.rs`
- `src/pipeline/assemble.rs`
- `src/pipeline/normalize.rs`
- `src/pipeline/isoform_processor.rs`
- `src/pipeline/large_genome_assembler.rs`
- `src/graph/`
- `src/kmer/`
- `src/gpu/`
- `tests/`
- `benches/`

Use these discovery commands:
- `git status --short --branch`
- `rg -n "isoform|paired|normalize|gpu|de Bruijn|butterfly|chrysalis|inchworm|transcript" src tests docs README.md`
- `cargo run -- --help`
- `cargo run -- assemble --help`
- `cargo run -- assemble-large --help`
- `cargo test --features gpu`

Use Trinity's documented architecture as the external reference: normalization, Inchworm, Chrysalis, Butterfly, paired-end evidence, optional genome-guided mode, and downstream assembly quality assessment.
</context>

<constraints>
- Rust first. Python and shell are acceptable only for benchmark orchestration, plotting, download glue, and reports.
- No rewrite from scratch. Work stage by stage with a mapped Trinity-equivalent contract and evidence.
- Delete dead paths when evidence proves they are obsolete, but do not delete major functionality without a benchmark-backed replacement.
- Do not claim Trinity parity from unit tests, synthetic smoke tests, or README language.
- Do not tune by discarding hard biological cases.
- Do not add CUDA, OpenCL changes, or other accelerator dependencies unless CPU output equivalence and benchmark benefit are measured.
- Keep one production codepath per behavior. Tests should exercise the same codepath users run.
- Preserve unrelated user changes.
- Benchmark claims must include exact commands, versions, hardware, input data, metrics, and artifacts.
</constraints>

<scorecard>
Primary score: Trinity replacement readiness, computed as the minimum readiness level across all core stages. A single weak stage blocks completion.

Core stages:
- input preparation and in silico normalization
- Inchworm-equivalent contig construction
- Chrysalis-equivalent clustering and de Bruijn graph partitioning
- Butterfly-equivalent isoform reconstruction using read and read-pair support
- paired-end evidence handling
- output/evaluation/reporting
- end-to-end CLI workflow

Stage levels:
- 0: absent or only stubbed
- 1: compiles and has unit tests
- 2: passes representative synthetic fixtures
- 3: matches Trinity on controlled simulated RNA-seq fixtures
- 4: matches Trinity on at least three real public RNA-seq datasets by biological metrics
- 5: matches or beats Trinity on biological metrics and runtime/resource usage across the approved benchmark panel

Passing threshold: every core stage at level 5.

Regression checks:
- biological metrics must not regress on frozen fixtures without a documented, approved tradeoff
- CPU and GPU outputs must remain equivalent when acceleration is claimed
- existing rescue checks must stay green

Scoring artifacts:
- `docs/trinity-parity-map.md`
- `docs/trinity-parity-report.md`
- `bench/trinity_parity/results/`
- `ATTEMPTS.md`

Stop condition: only stop when every `done_when` item is true and `docs/trinity-parity-report.md` supports the conclusion that Raptor is a credible Rust Trinity replacement.
</scorecard>

<done_when>
The goal is complete only when all items below are true:

- `docs/trinity-parity-map.md` maps Trinity stages to Raptor modules and marks every core stage complete with evidence links.
- `bench/trinity_parity/` contains reproducible scripts for downloading/preparing fixtures, running Trinity, running Raptor, and comparing outputs.
- `docs/trinity-parity-report.md` records the frozen benchmark panel, exact commands, versions, hardware, biological metrics, resource metrics, and pass/fail conclusions.
- Every core stage reaches scorecard level 5.
- Raptor end-to-end outputs are biologically comparable to Trinity across the approved benchmark panel.
- Raptor runtime and peak memory are no worse than Trinity by more than the approved tolerance on the benchmark panel, and at least one major stage is measurably faster or lower-memory.
- GPU acceleration, if used in final claims, produces equivalent biological outputs to CPU and has measured speed/resource benefit.
- `cargo clippy --all-targets --all-features -- -D warnings` passes.
- `cargo test --all-targets --all-features` passes.
- `./scripts/rescue_smoke.sh` passes.
- No README or documentation claims Trinity replacement status without linking to `docs/trinity-parity-report.md`.
</done_when>

<feedback_loop>
Fast loop:
- Run after focused changes.
- Commands: `cargo check`, focused tests for touched module, and the smallest active synthetic RNA-seq fixture.
- Expected runtime: seconds to a few minutes.
- Proxy validity: catches compile breakage and stage-local behavior regressions early.

Medium loop:
- Run before phase transitions.
- Commands: `cargo test --features gpu`, `./scripts/rescue_smoke.sh`, `./bench/gpu_kmer_baseline.sh`, and active stage benchmark scripts.
- Expected runtime: minutes.
- Proxy validity: checks rescue baseline, GPU health, and active-stage biological proxies.

Final loop:
- Run before claiming parity.
- Commands: full Trinity-vs-Raptor benchmark panel under `bench/trinity_parity/`, `cargo clippy --all-targets --all-features -- -D warnings`, and `cargo test --all-targets --all-features`.
- Expected runtime: hours or longer.
- Proxy validity: this is the real completion gate.
</feedback_loop>

<workflow>
1. Establish the parity map.
   - Create `docs/trinity-parity-map.md`.
   - Inventory Trinity stages and current Raptor modules.
   - Mark each stage as complete, partial, missing, or misleading.

2. Freeze the benchmark panel.
   - Create `bench/trinity_parity/`.
   - Add scripts to prepare synthetic and public RNA-seq datasets.
   - Add scripts to run Trinity and Raptor with exact captured versions.
   - Require explicit approval before changing the frozen panel after it is set.

3. Build the comparison harness.
   - Add metrics extraction for transcript FASTA, read representation, BUSCO or equivalent completeness, full-length recovery on truth-known sets, isoform precision/recall where possible, fusion/paralog indicators, runtime, RSS, disk, and GPU usage.
   - Save raw outputs outside git and curated summaries in docs.

4. Close stages one at a time.
   - Normalization.
   - Inchworm-equivalent contig construction.
   - Chrysalis-equivalent clustering/partitioning.
   - Butterfly-equivalent isoform reconstruction.
   - Paired-end evidence handling.
   - Output/evaluation/reporting.
   - End-to-end CLI workflow.

5. Optimize only after correctness is measurable.
   - Use CPU baselines first.
   - Use OpenCL/CUDA only when output equivalence is demonstrated.
   - Keep before/after metrics in `ATTEMPTS.md` and reports.

6. Final parity review.
   - Run final checks.
   - Update README only after parity evidence exists.
   - Do not mark complete until every `done_when` item is satisfied.
</workflow>

<working_memory>
Maintain these files throughout the goal:

- `PLAN.md`: current phase, current strategy, open decisions, next actions.
- `ATTEMPTS.md`: every meaningful implementation attempt, benchmark run, failed approach, metric change, and result.
- `NOTES.md`: durable discoveries, Trinity/Raptor architecture notes, blockers, and context needed after compaction.
- `CONTROL.md`: human operator panel.

Update cadence:
- Update `PLAN.md` at phase changes and after strategic pivots.
- Update `ATTEMPTS.md` after each implementation attempt, benchmark run, or failed experiment.
- Update `NOTES.md` whenever a discovery should survive context compaction.
- Reread `CONTROL.md` before phase changes, expensive benchmark runs, dependency changes, accelerator backend work, or parity claims.
</working_memory>

<human_control_surface>
Create and maintain `CONTROL.md` as the compact human operator panel for this goal.

Before each phase change, strategic pivot, expensive benchmark run, dependency change, accelerator backend change, or final parity claim, reread `CONTROL.md`. If it changed, summarize the relevant change in `PLAN.md` and adapt before proceeding.

`CONTROL.md` may narrow scope, change priorities, pause work, or require approval. It cannot silently weaken `done_when`, benchmark thresholds, or scorecard requirements.
</human_control_surface>

<verification_loop>
Focused verification:
- `cargo check`
- focused tests for touched modules
- active stage benchmark or fixture

Broad verification:
- `cargo clippy --all-targets --all-features -- -D warnings`
- `cargo test --all-targets --all-features`
- `cargo test --features gpu`
- `./scripts/rescue_smoke.sh`
- `./bench/gpu_kmer_baseline.sh`

Final verification:
- full frozen Trinity-vs-Raptor benchmark panel
- `docs/trinity-parity-report.md` regenerated from current artifacts
- user review of final parity claim
</verification_loop>

<execution_rules>
- Check git status before edits.
- Preserve unrelated user changes.
- Prefer `rg` over `grep` when available.
- Use the runtime patch/edit tool for manual edits when available.
- Read context files before implementation.
- Batch independent file reads in parallel when the runtime supports it.
- Keep the goal scorecard current: know the primary metric, passing threshold, regression checks, scoring method, and stop condition.
- Use the fastest representative feedback check while iterating; reserve slower checks for escalation points and final verification.
- Maintain `PLAN.md`, `ATTEMPTS.md`, `NOTES.md`, and `CONTROL.md`.
- Update `ATTEMPTS.md` after each meaningful approach so future iterations do not repeat work without new evidence.
- Run focused tests before broad tests.
- Do not paper over failures.
- Do not widen scope.
- Keep final answers concise and evidence-based.
</execution_rules>

<output_contract>
Final deliverables:
- Rust implementation changes that make Raptor a credible Trinity replacement.
- `docs/trinity-parity-map.md`
- `docs/trinity-parity-report.md`
- reproducible scripts under `bench/trinity_parity/`
- maintained `PLAN.md`, `ATTEMPTS.md`, `NOTES.md`, and `CONTROL.md`
- passing final verification commands

Completion response must state:
- final commit hash
- benchmark panel used
- biological metric summary
- runtime/resource comparison summary
- commands run
- remaining caveats, if any
</output_contract>
