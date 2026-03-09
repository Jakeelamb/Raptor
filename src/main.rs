mod accel;
mod cli;
mod cli_main;
mod dist;
mod eval;
mod gpu;
mod graph;
mod io;
mod kmer;
mod pipeline;
mod polish;
mod quant;
mod stats;
mod visualize;

use clap::Parser;
use cli_main::{Cli, Commands};
use ndarray::Array2;
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
use tracing::info;
use tracing_subscriber::FmtSubscriber;

fn main() {
    let subscriber = FmtSubscriber::builder()
        .with_max_level(tracing::Level::INFO)
        .finish();
    tracing::subscriber::set_global_default(subscriber).expect("Setting tracing default failed");

    let cli = Cli::parse();

    match cli.command {
        Commands::Normalize {
            input1,
            input2,
            output,
            gpu,
            threads: _,
            streaming,
            coverage_target: _,
            max_reads: _,
        } => {
            // Run the normalization pipeline
            println!("Running normalization pipeline");
            let start = std::time::Instant::now();

            if let Some(input2_path) = input2 {
                // For paired-end reads
                pipeline::normalize::normalize_paired(
                    &input1,
                    &input2_path,
                    &output,
                    gpu,
                    streaming,
                );
            } else {
                // For single-end reads
                pipeline::normalize::normalize_single(&input1, &output, gpu, streaming);
            }

            println!(
                "Normalization completed in {:.2}s",
                start.elapsed().as_secs_f32()
            );
        }

        Commands::Assemble {
            input,
            output,
            min_len,
            threads: _,
            gfa,
            gpu,
            adaptive_k,
            rle,
            distributed: _,
            buckets: _,
            gfa2,
            collapse_repeats,
            min_repeat_len,
            polish,
            polish_window,
            streaming: _,
            export_metadata,
            json_metadata,
            tsv_metadata,
            isoforms,
            gtf,
            counts_matrix,
            gff3,
            max_path_depth,
            min_confidence,
            min_path_len: _,
            dev_mode: _,
            compute_tpm,
            polish_isoforms,
            samples,
            min_tpm,
            polish_reads,
        } => {
            // Run the assembly pipeline
            println!("Running assembly pipeline");
            if gpu {
                println!("GPU acceleration requested");
            }
            let start = std::time::Instant::now();

            pipeline::assemble::assemble_reads_with_gpu(
                &input,
                &output,
                min_len,
                gfa,
                gfa2,
                adaptive_k,
                rle,
                collapse_repeats,
                min_repeat_len,
                polish,
                polish_window,
                false, // streaming - now handled internally
                export_metadata,
                json_metadata,
                tsv_metadata,
                isoforms,
                gtf,
                gff3,
                max_path_depth,
                min_confidence,
                compute_tpm,
                polish_isoforms,
                samples,
                min_tpm,
                polish_reads,
                counts_matrix,
                gpu,
            );

            println!(
                "Assembly completed in {:.2}s",
                start.elapsed().as_secs_f32()
            );
        }

        Commands::Stats {
            input,
            format,
            graph,
        } => {
            info!("Calculating assembly statistics for: {}", input);
            use stats::{calculate_graph_stats, calculate_stats, update_with_graph_stats};

            let mut stats = calculate_stats(&input);

            // Add graph complexity stats if requested
            if graph {
                info!("Calculating graph complexity stats");
                if let Some(graph_stats) = calculate_graph_stats(&input) {
                    update_with_graph_stats(&mut stats, &graph_stats);
                    println!("Graph analysis:");
                    println!("  Total paths: {}", graph_stats.total_paths);
                    println!("  Average path length: {:.2}", graph_stats.average_length);
                    println!("  Shared segments (branches): {}", graph_stats.branch_count);
                }
            }

            match format.as_str() {
                "json" => {
                    println!("{}", serde_json::to_string_pretty(&stats).unwrap());
                }
                "tsv" => {
                    println!("contigs\ttotal_len\tavg_len\tn50");
                    println!(
                        "{}\t{}\t{:.2}\t{}",
                        stats.total_contigs, stats.total_length, stats.average_length, stats.n50
                    );
                }
                _ => eprintln!("Unsupported format: {}", format),
            }
        }

        Commands::Benchmark { input, k, threads } => {
            info!(
                "Starting benchmark: input = {}, k = {}, threads = {}",
                input, k, threads
            );

            ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .expect("Failed to build thread pool");

            crate::cli::benchmark::benchmark_kmer_counting(&input, k);
        }

        Commands::Isoform {
            input,
            expression,
            output,
            min_confidence,
            max_depth,
            formats,
            threads,
            stats,
            filter_similar,
            similarity_threshold,
            merge_similar,
        } => {
            info!(
                "Starting isoform reconstruction: input = {}, expression = {}, output = {}",
                input, expression, output
            );

            ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .expect("Failed to build thread pool");

            match io::isoform_runner::run_isoform_reconstruction(
                &input,
                &expression,
                &output,
                min_confidence,
                max_depth,
                &formats,
                stats,
                if filter_similar {
                    Some(similarity_threshold)
                } else {
                    None
                },
                merge_similar,
            ) {
                Ok(_) => info!("Isoform reconstruction completed successfully"),
                Err(e) => eprintln!("Error during isoform reconstruction: {}", e),
            }
        }

        Commands::DiffExp {
            matrix,
            group_a,
            group_b,
            output,
            p_value,
            fold_change,
        } => {
            info!("Performing differential expression analysis");
            info!("Reading counts matrix: {}", matrix);

            // Read the counts matrix
            let content = match std::fs::read_to_string(&matrix) {
                Ok(c) => c,
                Err(e) => {
                    eprintln!("Error reading counts matrix file: {}", e);
                    return;
                }
            };

            // Parse the matrix
            let mut lines = content.lines();
            let header = match lines.next() {
                Some(h) => h,
                None => {
                    eprintln!("Error: empty counts matrix file");
                    return;
                }
            };

            // Extract sample names from header
            let samples: Vec<&str> = header.split('\t').skip(1).collect();

            // Create a map to store TPM values for each sample
            let mut sample_tpms: HashMap<String, Vec<f64>> = HashMap::new();

            // Initialize empty vectors for each sample
            for &sample in &samples {
                sample_tpms.insert(sample.to_string(), Vec::new());
            }

            // Extract transcript IDs and TPM values
            let mut transcript_ids = Vec::new();
            for line in lines {
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() < 2 {
                    continue;
                }

                let transcript_id = parts[0].to_string();
                transcript_ids.push(transcript_id);

                for (i, &sample) in samples.iter().enumerate() {
                    if i + 1 < parts.len() {
                        if let Ok(tpm) = parts[i + 1].parse::<f64>() {
                            sample_tpms.get_mut(sample).unwrap().push(tpm);
                        } else {
                            sample_tpms.get_mut(sample).unwrap().push(0.0);
                        }
                    } else {
                        sample_tpms.get_mut(sample).unwrap().push(0.0);
                    }
                }
            }

            // Parse group sample names
            let group_a_samples: Vec<&str> = group_a.split(',').collect();
            let group_b_samples: Vec<&str> = group_b.split(',').collect();

            info!("Group A samples: {:?}", group_a_samples);
            info!("Group B samples: {:?}", group_b_samples);

            // Verify samples exist in the matrix
            for &sample in &group_a_samples {
                if !sample_tpms.contains_key(sample) {
                    eprintln!(
                        "Error: Sample '{}' from group A not found in matrix",
                        sample
                    );
                    return;
                }
            }

            for &sample in &group_b_samples {
                if !sample_tpms.contains_key(sample) {
                    eprintln!(
                        "Error: Sample '{}' from group B not found in matrix",
                        sample
                    );
                    return;
                }
            }

            // Define helper functions
            fn variance(x: &[f64]) -> f64 {
                let mean = x.iter().sum::<f64>() / x.len() as f64;
                let var = x.iter().map(|&v| (v - mean).powi(2)).sum::<f64>() / x.len() as f64;
                var
            }

            fn t_test(a: &[f64], b: &[f64]) -> f64 {
                if a.is_empty() || b.is_empty() {
                    return 1.0;
                }

                let mean_a = a.iter().sum::<f64>() / a.len() as f64;
                let mean_b = b.iter().sum::<f64>() / b.len() as f64;

                let var_a = variance(a);
                let var_b = variance(b);

                let pooled_var = (var_a / a.len() as f64) + (var_b / b.len() as f64);
                if pooled_var == 0.0 {
                    return 1.0;
                }

                let t = (mean_a - mean_b) / pooled_var.sqrt();

                // Simple p-value approximation based on t-statistic
                // In a real implementation, we would use a proper CDF function
                let p = 1.0 / (1.0 + t.abs());
                p
            }

            #[derive(Debug)]
            struct DiffExpResult {
                transcript: String,
                fold_change: f64,
                p_value: f64,
            }

            // Perform differential expression analysis
            let mut results = Vec::new();

            for (i, transcript) in transcript_ids.iter().enumerate() {
                let mut group_a_values = Vec::new();
                let mut group_b_values = Vec::new();

                for sample in &group_a_samples {
                    let value = sample_tpms.get(*sample).unwrap()[i];
                    group_a_values.push(value);
                }

                for sample in &group_b_samples {
                    let value = sample_tpms.get(*sample).unwrap()[i];
                    group_b_values.push(value);
                }

                let mean_a = group_a_values.iter().sum::<f64>() / group_a_values.len() as f64;
                let mean_b = group_b_values.iter().sum::<f64>() / group_b_values.len() as f64;

                // Handle potential zeros in means
                let fold_change = if mean_a > 0.0 && mean_b > 0.0 {
                    (mean_a / mean_b).log2()
                } else if mean_a > 0.0 {
                    f64::INFINITY
                } else if mean_b > 0.0 {
                    f64::NEG_INFINITY
                } else {
                    0.0 // Both are zero
                };

                let p = t_test(&group_a_values, &group_b_values);

                results.push(DiffExpResult {
                    transcript: transcript.clone(),
                    fold_change,
                    p_value: p,
                });
            }

            // Filter by significance thresholds
            let significant_results: Vec<_> = results
                .into_iter()
                .filter(|r| r.p_value < p_value && r.fold_change.abs() >= fold_change)
                .collect();

            // Output results
            let mut output_file = match std::fs::File::create(&output) {
                Ok(f) => f,
                Err(e) => {
                    eprintln!("Error creating output file: {}", e);
                    return;
                }
            };

            use std::io::Write;
            writeln!(output_file, "transcript\tlog2_fold_change\tp_value").unwrap();

            for result in significant_results {
                writeln!(
                    output_file,
                    "{}\t{:.4}\t{:.6}",
                    result.transcript, result.fold_change, result.p_value
                )
                .unwrap();
            }

            info!(
                "Differential expression analysis completed. Results written to: {}",
                output
            );
        }

        Commands::GtfCompare {
            truth,
            predicted,
            output,
        } => {
            info!("Comparing GTF files");

            // Compare GTF files
            match eval::compare_gtf(&truth, &predicted) {
                Ok((tp, fp, fn_)) => {
                    let (precision, recall, f1) = eval::calculate_metrics(tp, fp, fn_);

                    let results = format!(
                        "Exon-level metrics:\n\
                        True positives: {}\n\
                        False positives: {}\n\
                        False negatives: {}\n\
                        Precision: {:.4}\n\
                        Recall: {:.4}\n\
                        F1 score: {:.4}",
                        tp, fp, fn_, precision, recall, f1
                    );

                    println!("{}", results);

                    if let Some(output_path) = output {
                        if let Err(e) = std::fs::write(&output_path, results) {
                            eprintln!("Error writing output file: {}", e);
                        }
                    }
                }
                Err(e) => {
                    eprintln!("Error comparing GTF files: {}", e);
                }
            }
        }

        Commands::Eval {
            truth,
            pred,
            output,
        } => {
            info!("Evaluating assembly results");

            // This is a simplified evaluation
            if let Err(e) = eval::evaluate_assembly(&truth, &pred, output.as_deref()) {
                eprintln!("Error during evaluation: {}", e);
            }
        }

        Commands::Visualize {
            matrix,
            output: _,
            heatmap,
            pca,
            components,
        } => {
            info!("Visualizing expression data");

            // Simple validation
            if components < 2 {
                eprintln!("Error: Number of PCA components must be at least 2");
                return;
            }

            // Load the matrix file
            let content = match std::fs::read_to_string(&matrix) {
                Ok(c) => c,
                Err(e) => {
                    eprintln!("Error reading matrix file: {}", e);
                    return;
                }
            };

            // Parse the matrix
            let mut lines = content.lines();
            let header = match lines.next() {
                Some(h) => h,
                None => {
                    eprintln!("Error: empty matrix file");
                    return;
                }
            };

            let sample_names: Vec<String> =
                header.split('\t').skip(1).map(|s| s.to_string()).collect();
            let mut transcript_ids = Vec::new();
            let mut data = Vec::new();

            for line in lines {
                let parts: Vec<&str> = line.split('\t').collect();
                if parts.len() < 2 {
                    continue;
                }

                transcript_ids.push(parts[0].to_string());

                let mut row = Vec::new();
                for &value in parts.iter().skip(1) {
                    match value.parse::<f64>() {
                        Ok(v) => row.push(v),
                        Err(_) => row.push(0.0),
                    }
                }

                // Ensure the row has the right number of elements
                while row.len() < sample_names.len() {
                    row.push(0.0);
                }

                data.push(row);
            }

            // Create output
            if data.is_empty() {
                eprintln!("Error: No data found in matrix");
                return;
            }

            // Convert to ndarray format for our visualization functions
            let rows = data.len();
            let cols = if rows > 0 { data[0].len() } else { 0 };
            let mut matrix_data = Array2::zeros((rows, cols));

            for i in 0..rows {
                for j in 0..cols {
                    matrix_data[[i, j]] = data[i][j];
                }
            }

            if let Some(heatmap_path) = heatmap {
                match visualize::plot::plot_heatmap_with_labels(
                    &matrix_data,
                    &transcript_ids,
                    &sample_names,
                    &heatmap_path,
                ) {
                    Ok(()) => info!("Heatmap saved to {}", heatmap_path),
                    Err(e) => eprintln!("Error creating heatmap: {}", e),
                }
            }

            if let Some(pca_path) = pca {
                // Convert to HashMap format for our PCA function
                let mut sample_data = HashMap::new();
                for (j, sample) in sample_names.iter().enumerate() {
                    let mut values = Vec::new();
                    for i in 0..rows {
                        values.push(matrix_data[[i, j]]);
                    }
                    sample_data.insert(sample.clone(), values);
                }

                match visualize::pca::plot_pca(&sample_data, &pca_path) {
                    Ok(()) => info!("PCA plot saved to {}", pca_path),
                    Err(e) => eprintln!("Error creating PCA plot: {}", e),
                }
            }

            info!("Visualization completed");
        }

        Commands::Traverse {
            input,
            segments,
            output,
            formats,
            include_edges,
            visualize,
            metadata,
        } => {
            info!("Traversing paths in GFA file");

            match io::gfa_utils::traverse_paths(
                &input,
                &segments,
                &output,
                &formats,
                include_edges,
                visualize,
                metadata,
            ) {
                Ok(count) => info!("Successfully traversed {} paths", count),
                Err(e) => eprintln!("Error traversing paths: {}", e),
            }
        }

        Commands::AssembleLarge {
            input,
            input2,
            output,
            kmer,
            min_count,
            min_contig,
            threads,
            temp_dir,
            num_buckets,
            max_tip_len,
            max_bubble_len,
            scaffold,
            min_scaffold_links,
            polish,
            polish_iterations,
            long_reads,
            min_long_read_len,
            compress_buckets: _,
        } => {
            println!("Running large genome assembly pipeline");
            let start = std::time::Instant::now();

            // Configure thread pool
            ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .expect("Failed to build thread pool");

            // Build configuration
            let config = pipeline::large_genome_assembler::LargeGenomeConfig {
                k: kmer,
                min_count,
                min_contig_len: min_contig,
                num_buckets,
                temp_dir,
                max_tip_len,
                max_bubble_len,
            };

            let assembler = pipeline::large_genome_assembler::LargeGenomeAssembler::new(config);

            // Check for paired-end input
            let is_paired = input2.is_some();
            if is_paired {
                info!(
                    "Paired-end mode: {} and {}",
                    input,
                    input2.as_ref().unwrap()
                );
            }

            // Run assembly
            let assembly_result = if let Some(input2_path) = input2.as_deref() {
                assembler.assemble_paired(&input, input2_path, &output)
            } else {
                assembler.assemble(&input, &output)
            };

            match assembly_result {
                Ok(stats) => {
                    println!("{}", stats);

                    // Optional: scaffolding with paired-end reads
                    if scaffold && is_paired {
                        info!("Scaffolding with paired-end reads...");
                        let scaffold_output = output
                            .replace(".fa", ".scaffolds.fa")
                            .replace(".fasta", ".scaffolds.fasta");
                        match pipeline::scaffolder::scaffold_contigs(
                            &output,
                            &input,
                            input2.as_ref().unwrap(),
                            &scaffold_output,
                            min_scaffold_links,
                        ) {
                            Ok(scaffold_stats) => {
                                println!("Scaffolding completed:");
                                println!("  Scaffolds: {}", scaffold_stats.num_scaffolds);
                                println!("  N50: {} bp", scaffold_stats.scaffold_n50);
                            }
                            Err(e) => eprintln!("Scaffolding failed: {}", e),
                        }
                    } else if scaffold && !is_paired {
                        eprintln!("Warning: Scaffolding requires paired-end reads (--input2)");
                    }

                    // Optional: polishing
                    if polish {
                        info!("Polishing contigs...");
                        let polish_output = output
                            .replace(".fa", ".polished.fa")
                            .replace(".fasta", ".polished.fasta");
                        match pipeline::polisher::polish_contigs(
                            &output,
                            &input,
                            input2.as_deref(),
                            &polish_output,
                            polish_iterations,
                        ) {
                            Ok(polish_stats) => {
                                println!("Polishing completed:");
                                println!("  Corrections: {}", polish_stats.corrections);
                            }
                            Err(e) => eprintln!("Polishing failed: {}", e),
                        }
                    }

                    // Optional: long read integration
                    if let Some(long_reads_path) = long_reads {
                        info!("Integrating long reads from {}...", long_reads_path);
                        let long_read_output = output
                            .replace(".fa", ".hybrid.fa")
                            .replace(".fasta", ".hybrid.fasta");

                        let lr_config = pipeline::long_read_integration::LongReadConfig {
                            min_read_length: min_long_read_len,
                            ..Default::default()
                        };

                        match pipeline::long_read_integration::integrate_long_reads(
                            &output,
                            &long_reads_path,
                            &long_read_output,
                            lr_config,
                        ) {
                            Ok(lr_stats) => {
                                println!("Long read integration completed:");
                                println!("  Long reads mapped: {}", lr_stats.long_reads_mapped);
                                println!("  Contigs extended: {}", lr_stats.contigs_extended);
                                println!("  Contigs joined: {}", lr_stats.contigs_joined / 2);
                            }
                            Err(e) => eprintln!("Long read integration failed: {}", e),
                        }
                    }
                }
                Err(e) => eprintln!("Assembly failed: {}", e),
            }

            println!(
                "Large genome assembly completed in {:.2}s",
                start.elapsed().as_secs_f32()
            );
        }
    }
}
