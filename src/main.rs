mod cli_main;
mod io;
mod kmer;
mod graph;
mod accel;
mod dist;
mod pipeline;
mod gpu;
mod stats;
mod cli;
mod quant;
mod polish;
mod eval;
mod visualize;

use clap::Parser;
use tracing::info;
use tracing_subscriber::FmtSubscriber;
use rayon::ThreadPoolBuilder;
use cli_main::{Cli, Commands};

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
            threads, 
            streaming, 
            coverage_target, 
            max_reads 
        } => {
            // Run the normalization pipeline
            println!("Running normalization pipeline");
            let start = std::time::Instant::now();
            pipeline::normalize::normalize_paired_reads(
                input1, 
                input2, 
                output, 
                coverage_target,
                max_reads,
                *gpu, 
                *threads, 
                *streaming
            )?;
            println!("Normalization completed in {:.2}s", start.elapsed().as_secs_f32());
            Ok(())
        },

        Commands::Assemble { 
            input, 
            output, 
            min_len, 
            min_coverage, 
            min_confidence, 
            threads, 
            k_size, 
            num_buckets, 
            gfa, 
            streaming, 
            polish, 
            distributed, 
            isoforms, 
            json_metadata,
            export_metadata,
            alignment,
            debug,
            export_graph, 
            adaptive_k,
            use_rle,
            gtf,
            gff3,
            compute_tpm,
            min_kmers,
            min_abundance,
            min_pairs,
            output_gfa,
            counts_matrix,
            min_path_len,
            dev_mode
        } => {
            // Run the assembly pipeline
            println!("Running assembly pipeline");
            let start = std::time::Instant::now();
            pipeline::assemble::assemble_reads(
                input, 
                output, 
                *min_len, 
                *min_coverage, 
                *min_confidence,
                *threads, 
                *k_size, 
                *num_buckets, 
                *gfa, 
                *streaming, 
                *polish, 
                *distributed, 
                *isoforms, 
                json_metadata,
                *export_metadata,
                alignment,
                *debug,
                export_graph,
                *adaptive_k,
                *use_rle,
                gtf,
                gff3,
                *compute_tpm,
                *min_kmers,
                *min_abundance,
                *min_pairs,
                *min_path_len
            )?;
            println!("Assembly completed in {:.2}s", start.elapsed().as_secs_f32());
            Ok(())
        },
        
        Commands::Stats { input, format, graph } => {
            info!("Calculating assembly statistics for: {}", input);
            use stats::{calculate_stats, calculate_graph_stats, update_with_graph_stats};

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
                    println!("{}\t{}\t{:.2}\t{}", 
                             stats.total_contigs, 
                             stats.total_length, 
                             stats.average_length, 
                             stats.n50);
                }
                _ => eprintln!("Unsupported format: {}", format),
            }
        }

        Commands::Benchmark { input, k, threads } => {
            info!("Starting benchmark: input = {}, k = {}, threads = {}", input, k, threads);
            
            ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .expect("Failed to build thread pool");
                
            crate::cli::benchmark::benchmark_kmer_counting(&input, k);
        }
        
        Commands::Isoform { input, expression, output, min_confidence, max_depth, formats, threads, stats, filter_similar, similarity_threshold, merge_similar } => {
            info!("Starting isoform reconstruction: input = {}, expression = {}, output = {}", 
                  input, expression, output);
            
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
                if filter_similar { Some(similarity_threshold) } else { None },
                merge_similar
            ) {
                Ok(_) => info!("Isoform reconstruction completed successfully"),
                Err(e) => eprintln!("Error during isoform reconstruction: {}", e),
            }
        }
        
        Commands::DiffExp { matrix, group_a, group_b, output, p_value, fold_change } => {
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
            use std::collections::HashMap;
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
            
            // Parse group A and group B sample names
            let group_a_samples: Vec<&str> = group_a.split(',').collect();
            let group_b_samples: Vec<&str> = group_b.split(',').collect();
            
            info!("Group A samples: {:?}", group_a_samples);
            info!("Group B samples: {:?}", group_b_samples);
            
            // Validate that all group samples exist in the matrix
            for &sample in &group_a_samples {
                if !sample_tpms.contains_key(sample) {
                    eprintln!("Error: sample '{}' not found in counts matrix", sample);
                    return;
                }
            }
            
            for &sample in &group_b_samples {
                if !sample_tpms.contains_key(sample) {
                    eprintln!("Error: sample '{}' not found in counts matrix", sample);
                    return;
                }
            }
            
            // Perform differential expression analysis
            let len = sample_tpms.values().next().unwrap().len();
            let mut results = Vec::new();
            
            // Trait for log2 calculation
            trait Log2 {
                fn log2(self) -> f64;
            }
            
            impl Log2 for f64 {
                fn log2(self) -> f64 {
                    self.ln() / 2f64.ln()
                }
            }
            
            // Function to calculate variance
            fn variance(x: &[f64]) -> f64 {
                let mean = x.iter().sum::<f64>() / x.len() as f64;
                x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / x.len() as f64
            }
            
            // Simple t-test implementation
            fn t_test(a: &[f64], b: &[f64]) -> f64 {
                let var_a = variance(a);
                let var_b = variance(b);
                let mean_a = a.iter().sum::<f64>() / a.len() as f64;
                let mean_b = b.iter().sum::<f64>() / b.len() as f64;
                
                let t = (mean_a - mean_b).abs()
                    / ((var_a / a.len() as f64) + (var_b / b.len() as f64)).sqrt();
                
                // Simple approximation of p-value
                1.0 / (1.0 + t * t)
            }
            
            // DiffExp struct for storage
            struct DiffExpResult {
                transcript: String,
                fold_change: f64,
                p_value: f64,
            }
            
            for i in 0..len {
                let a_vals: Vec<f64> = group_a_samples.iter().map(|&s| sample_tpms[s][i]).collect();
                let b_vals: Vec<f64> = group_b_samples.iter().map(|&s| sample_tpms[s][i]).collect();
                
                let mean_a = a_vals.iter().sum::<f64>() / a_vals.len() as f64;
                let mean_b = b_vals.iter().sum::<f64>() / b_vals.len() as f64;
                let fold_change = (mean_b + 1.0) / (mean_a + 1.0); // add pseudo-count
                let p_value = t_test(&a_vals, &b_vals);
                
                results.push(DiffExpResult {
                    transcript: transcript_ids.get(i).unwrap_or(&format!("transcript_{}", i)).clone(),
                    fold_change,
                    p_value,
                });
            }
            
            // Filter results by p-value and fold change
            let filtered_results: Vec<_> = results.into_iter()
                .filter(|r| r.p_value < p_value && (r.fold_change.log2().abs() >= fold_change))
                .collect();
            
            // Write results to output file
            info!("Writing {} differentially expressed transcripts to {}", filtered_results.len(), output);
            
            // Write output file
            let mut file = match std::fs::File::create(&output) {
                Ok(f) => f,
                Err(e) => {
                    eprintln!("Error creating output file: {}", e);
                    return;
                }
            };
            
            use std::io::Write;
            writeln!(file, "transcript_id\tlog2FC\tp_value").unwrap();
            for r in &filtered_results {
                writeln!(file, "{}\t{:.3}\t{:.4}", r.transcript, r.fold_change.log2(), r.p_value).unwrap();
            }
            
            println!("Differential expression analysis complete.");
            println!("Found {} significantly differentially expressed transcripts.", filtered_results.len());
            
            // Generate PCA plot if sample TPMs are available
            if !sample_tpms.is_empty() {
                info!("Generating PCA visualization of samples");
                let pca_path = format!("{}_pca.svg", output);
                
                if let Err(e) = crate::visualize::pca::plot_pca(&sample_tpms, &pca_path) {
                    eprintln!("Error generating PCA plot: {}", e);
                } else {
                    info!("PCA visualization saved to {}", pca_path);
                }
            }
        }
        
        Commands::GtfCompare { truth, predicted, output } => {
            info!("Comparing predicted GTF with reference GTF");
            info!("Truth GTF: {}", truth);
            info!("Predicted GTF: {}", predicted);
            
            // Implementation without using unresolved imports
            
            // Helper function to load GTF ranges
            fn load_gtf_ranges(path: &str) -> std::io::Result<std::collections::HashSet<(String, String, usize, usize)>> {
                use std::collections::HashSet;
                use std::fs::File;
                use std::io::{BufRead, BufReader};
                
                let file = File::open(path)?;
                let reader = BufReader::new(file);
                let mut ranges = HashSet::new();
                
                for line_result in reader.lines() {
                    let line = line_result?;
                    if line.starts_with('#') {
                        continue; // Skip comment lines
                    }
                    
                    let cols: Vec<&str> = line.split('\t').collect();
                    if cols.len() < 9 {
                        continue; // Skip malformed lines
                    }
                    
                    // We only care about exon features for this comparison
                    if cols[2] != "exon" {
                        continue;
                    }
                    
                    // Parse the attributes field to extract transcript_id
                    let attrs = cols[8];
                    let transcript_id = if let Some(pos) = attrs.find("transcript_id") {
                        let remaining = &attrs[pos + 14..]; // Skip "transcript_id \""
                        if let Some(end_pos) = remaining.find('"') {
                            remaining[..end_pos].to_string()
                        } else {
                            continue; // Malformed attribute
                        }
                    } else {
                        continue; // No transcript ID found
                    };
                    
                    // Parse start and end positions
                    let start = match cols[3].parse::<usize>() {
                        Ok(s) => s,
                        Err(_) => continue,
                    };
                    
                    let end = match cols[4].parse::<usize>() {
                        Ok(e) => e,
                        Err(_) => continue,
                    };
                    
                    // Add to the ranges set
                    ranges.insert((cols[0].to_string(), transcript_id, start, end));
                }
                
                Ok(ranges)
            }
            
            // Function to calculate precision, recall, and F1 score
            fn calculate_metrics(tp: usize, fp: usize, fn_: usize) -> (f64, f64, f64) {
                let precision = if tp + fp > 0 {
                    tp as f64 / (tp + fp) as f64
                } else {
                    0.0
                };
                
                let recall = if tp + fn_ > 0 {
                    tp as f64 / (tp + fn_) as f64
                } else {
                    0.0
                };
                
                let f1 = if precision + recall > 0.0 {
                    2.0 * precision * recall / (precision + recall)
                } else {
                    0.0
                };
                
                (precision, recall, f1)
            }
            
            // Compare GTF files
            match (load_gtf_ranges(&truth), load_gtf_ranges(&predicted)) {
                (Ok(true_ranges), Ok(pred_ranges)) => {
                    // Count true positives, false positives, and false negatives
                    let tp = pred_ranges.intersection(&true_ranges).count();
                    let fp = pred_ranges.difference(&true_ranges).count();
                    let fn_ = true_ranges.difference(&pred_ranges).count();
                    
                    let (precision, recall, f1) = calculate_metrics(tp, fp, fn_);
                    
                    println!("GTF Comparison Results:");
                    println!("True Positives: {}", tp);
                    println!("False Positives: {}", fp);
                    println!("False Negatives: {}", fn_);
                    println!("Precision: {:.4}", precision);
                    println!("Recall: {:.4}", recall);
                    println!("F1 Score: {:.4}", f1);
                    
                    // Write results to output file if specified
                    if let Some(out_path) = output {
                        info!("Writing comparison results to {}", out_path);
                        let mut file = match std::fs::File::create(&out_path) {
                            Ok(f) => f,
                            Err(e) => {
                                eprintln!("Error creating output file: {}", e);
                                return;
                            }
                        };
                        
                        use std::io::Write;
                        writeln!(file, "Metric\tValue").unwrap();
                        writeln!(file, "true_positives\t{}", tp).unwrap();
                        writeln!(file, "false_positives\t{}", fp).unwrap();
                        writeln!(file, "false_negatives\t{}", fn_).unwrap();
                        writeln!(file, "precision\t{:.6}", precision).unwrap();
                        writeln!(file, "recall\t{:.6}", recall).unwrap();
                        writeln!(file, "f1_score\t{:.6}", f1).unwrap();
                    }
                }
                (Err(e), _) => {
                    eprintln!("Error reading truth GTF file: {}", e);
                }
                (_, Err(e)) => {
                    eprintln!("Error reading predicted GTF file: {}", e);
                }
            }
        }
        
        Commands::Eval { truth, pred, output } => {
            info!("Evaluating assembled transcripts against ground truth");
            info!("Truth file: {}", truth);
            info!("Predicted file: {}", pred);
            
            use crate::eval::gtf_compare::{compare_gtf, calculate_metrics};
            
            match compare_gtf(&truth, &pred) {
                Ok((tp, fp, fn_)) => {
                    let (precision, recall, f1) = calculate_metrics(tp, fp, fn_);
                    
                    println!("Evaluation Results:");
                    println!("  True Positives: {}", tp);
                    println!("  False Positives: {}", fp);
                    println!("  False Negatives: {}", fn_);
                    println!("  Precision: {:.4}", precision);
                    println!("  Recall: {:.4}", recall);
                    println!("  F1 Score: {:.4}", f1);
                    
                    // Write output if requested
                    if let Some(output_path) = output {
                        let content = format!(
                            "metric\tvalue\ntp\t{}\nfp\t{}\nfn\t{}\nprecision\t{:.4}\nrecall\t{:.4}\nf1\t{:.4}\n",
                            tp, fp, fn_, precision, recall, f1
                        );
                        
                        if let Err(e) = std::fs::write(&output_path, content) {
                            eprintln!("Error writing evaluation results to {}: {}", output_path, e);
                        } else {
                            info!("Evaluation results written to {}", output_path);
                        }
                    }
                }
                Err(e) => {
                    eprintln!("Error comparing GTF files: {}", e);
                }
            }
        }
        
        Commands::Visualize { matrix, output, heatmap, pca, components } => {
            info!("Generating visualization for TPM matrix: {}", matrix);
            
            use crate::quant::matrix::read_tpm_matrix;
            use crate::visualize::pca::{plot_pca, compute_pca, plot_pca_simple};
            use crate::visualize::plot::{plot_heatmap, plot_heatmap_with_labels};
            use ndarray::Array2;
            
            match read_tpm_matrix(&matrix) {
                Ok(samples) => {
                    info!("Loaded TPM data for {} samples", samples.len());
                    
                    // Legacy output parameter (backwards compatibility)
                    match plot_pca(&samples, &output) {
                        Ok(_) => {
                            info!("PCA visualization saved to {}", output);
                        },
                        Err(e) => {
                            eprintln!("Error generating PCA plot: {}", e);
                        }
                    }
                    
                    // Generate heatmap if requested
                    if let Some(heatmap_path) = heatmap {
                        // Convert HashMap to ndarray
                        let sample_names: Vec<String> = samples.keys().cloned().collect();
                        let n_samples = sample_names.len();
                        
                        if n_samples > 0 {
                            let first_sample = &samples[&sample_names[0]];
                            let n_features = first_sample.len();
                            
                            let mut matrix = Array2::zeros((n_samples, n_features));
                            for (i, name) in sample_names.iter().enumerate() {
                                for (j, val) in samples[name].iter().enumerate() {
                                    matrix[[i, j]] = *val;
                                }
                            }
                            
                            // Get transcript IDs as column labels
                            let transcript_labels: Vec<String> = (0..n_features)
                                .map(|i| format!("transcript_{}", i+1))
                                .collect();
                            
                            match plot_heatmap_with_labels(&matrix, &sample_names, &transcript_labels, &heatmap_path) {
                                Ok(_) => {
                                    info!("Heatmap visualization saved to {}", heatmap_path);
                                },
                                Err(e) => {
                                    eprintln!("Error generating heatmap: {}", e);
                                }
                            }
                        } else {
                            eprintln!("No samples found in the TPM matrix");
                        }
                    }
                    
                    // Generate PCA plot if requested
                    if let Some(pca_path) = pca {
                        // Convert HashMap to ndarray
                        let sample_names: Vec<String> = samples.keys().cloned().collect();
                        let n_samples = sample_names.len();
                        
                        if n_samples > 0 {
                            let first_sample = &samples[&sample_names[0]];
                            let n_features = first_sample.len();
                            
                            let mut matrix = Array2::zeros((n_samples, n_features));
                            for (i, name) in sample_names.iter().enumerate() {
                                for (j, val) in samples[name].iter().enumerate() {
                                    matrix[[i, j]] = *val;
                                }
                            }
                            
                            // Compute PCA and generate plot
                            let pca_result = compute_pca(&matrix, components);
                            match plot_pca_simple(&pca_result, &pca_path) {
                                Ok(_) => {
                                    info!("PCA visualization saved to {}", pca_path);
                                },
                                Err(e) => {
                                    eprintln!("Error generating PCA plot: {}", e);
                                }
                            }
                        } else {
                            eprintln!("No samples found in the TPM matrix");
                        }
                    }
                },
                Err(e) => {
                    eprintln!("Error reading TPM matrix: {}", e);
                }
            }
        }

        Commands::Traverse { input, segments, output, formats, include_edges: _, visualize, metadata: _ } => {
            info!("Traversing paths in GFA file: {}", input);
            
            // Load GFA file
            let gfa_content = match std::fs::read_to_string(&input) {
                Ok(content) => content,
                Err(e) => {
                    eprintln!("Error reading GFA file: {}", e);
                    return;
                }
            };
            
            // Load segment sequences
            let segment_map = match graph::navigation::load_segment_sequences(&segments) {
                Ok(map) => map,
                Err(e) => {
                    eprintln!("Error loading segment sequences: {}", e);
                    return;
                }
            };
            
            // Parse GFA paths
            let gfa_lines: Vec<String> = gfa_content.lines().map(String::from).collect();
            let paths = graph::navigation::parse_gfa_paths(&gfa_lines);
            
            // Reconstruct path sequences
            let reconstructed_paths = graph::navigation::reconstruct_paths(&paths, &segment_map);
            
            // Process and export based on requested formats
            let format_list: Vec<&str> = formats.split(',').collect();
            
            for format in format_list {
                match format.trim() {
                    "fasta" => {
                        let fasta_path = format!("{}.fasta", output);
                        if let Err(e) = io::export::export_paths_to_fasta(&reconstructed_paths, &fasta_path) {
                            eprintln!("Error exporting to FASTA: {}", e);
                        } else {
                            info!("Exported path sequences to FASTA: {}", fasta_path);
                        }
                    },
                    "json" => {
                        let json_path = format!("{}.json", output);
                        
                        #[derive(serde::Serialize)]
                        struct SerializablePathMeta {
                            path_id: String,
                            segment_count: usize,
                            unique_segments: usize,
                            total_length: usize,
                            has_inversions: bool,
                        }
                        
                        // Create metadata for each path
                        let mut path_metadatas = Vec::new();
                        for (id, sequence) in &reconstructed_paths {
                            // Retrieve the original path segments
                            let empty_vec = Vec::new();
                            let segments = paths.get(id).unwrap_or(&empty_vec);
                            let segment_count = segments.len();
                            
                            // Count unique segments
                            let mut unique_segment_ids = std::collections::HashSet::new();
                            for (seg_id, _) in segments {
                                unique_segment_ids.insert(seg_id);
                            }
                            
                            // Check if path has inversions
                            let has_inversions = segments.iter().any(|(_, dir)| *dir == '-');
                            
                            // Create metadata
                            let meta = SerializablePathMeta {
                                path_id: id.clone(),
                                segment_count,
                                unique_segments: unique_segment_ids.len(),
                                total_length: sequence.len(),
                                has_inversions,
                            };
                            
                            path_metadatas.push(meta);
                        }
                        
                        let json_data = serde_json::to_string_pretty(&path_metadatas)
                            .expect("Failed to serialize path metadata");
                        
                        if let Err(e) = std::fs::write(&json_path, json_data) {
                            eprintln!("Error writing JSON metadata: {}", e);
                        } else {
                            info!("Exported path metadata to JSON: {}", json_path);
                        }
                    },
                    "dot" => {
                        if visualize {
                            // Create a dot graph for visualization
                            let dot_path = format!("{}.dot", output);
                            // Using a simplified approach since export_paths_to_dot is not available
                            let dot_content = format!(
                                "digraph G {{\n  label=\"Paths in {}\";\n}}\n",
                                input
                            );
                            if let Err(e) = std::fs::write(&dot_path, dot_content) {
                                eprintln!("Error creating DOT file: {}", e);
                            } else {
                                info!("Created simple DOT graph at {}", dot_path);
                                info!("Note: Full DOT graph export is not implemented");
                            }
                        }
                    },
                    _ => {
                        eprintln!("Unsupported export format: {}", format);
                    }
                }
            }
            
            info!("Path traversal complete. Found {} paths, exported to requested formats.", paths.len());
        }
    }
}

