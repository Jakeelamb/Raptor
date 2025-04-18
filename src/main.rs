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
use std::collections::HashMap;
use ndarray::Array2;

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

            if let Some(input2_path) = input2 {
                // For paired-end reads
                pipeline::normalize::normalize_paired(&input1, &input2_path, &output, gpu, streaming);
            } else {
                // For single-end reads
                pipeline::normalize::normalize_single(&input1, &output, gpu, streaming);
        }
            
            println!("Normalization completed in {:.2}s", start.elapsed().as_secs_f32());
        },

        Commands::Assemble { 
            input, 
            output, 
            min_len, 
            threads, 
            gfa, 
            adaptive_k,
            rle,
            distributed,
            buckets,
            gfa2,
            collapse_repeats,
            min_repeat_len,
            polish,
            polish_window,
            streaming,
            export_metadata,
            json_metadata,
            tsv_metadata,
            isoforms, 
            gtf, 
            counts_matrix, 
            gff3,
            max_path_depth, 
            min_confidence,
            min_path_len,
            dev_mode, 
            compute_tpm,
            polish_isoforms,
            samples,
            min_tpm,
            polish_reads 
        } => {
            // Run the assembly pipeline
            println!("Running assembly pipeline");
            let start = std::time::Instant::now();
            
            pipeline::assemble::assemble_reads_old(
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
                streaming,
                export_metadata,
                json_metadata.map(String::from),
                tsv_metadata.map(String::from),
                isoforms,
                gtf.map(String::from),
                gff3.map(String::from),
                max_path_depth,
                min_confidence,
                compute_tpm,
                polish_isoforms,
                samples.map(String::from),
                min_tpm,
                polish_reads.map(String::from),
                    counts_matrix
                );
            
            println!("Assembly completed in {:.2}s", start.elapsed().as_secs_f32());
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
                    eprintln!("Error: Sample '{}' from group A not found in matrix", sample);
                    return;
                }
            }
            
            for &sample in &group_b_samples {
                if !sample_tpms.contains_key(sample) {
                    eprintln!("Error: Sample '{}' from group B not found in matrix", sample);
                    return;
                }
            }
            
            // Define helper functions
            trait Log2 {
                fn log2(self) -> f64;
            }
            
            impl Log2 for f64 {
                fn log2(self) -> f64 {
                    self.log2()
                }
            }
            
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
                    std::f64::INFINITY
                } else if mean_b > 0.0 {
                    std::f64::NEG_INFINITY
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
            let significant_results: Vec<_> = results.into_iter()
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
                    result.transcript, 
                    result.fold_change, 
                    result.p_value
                ).unwrap();
            }
            
            info!("Differential expression analysis completed. Results written to: {}", output);
        }
        
        Commands::GtfCompare { truth, predicted, output } => {
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
        
        Commands::Eval { truth, pred, output } => {
            info!("Evaluating assembly results");
            
            // This is a simplified evaluation
            if let Err(e) = eval::evaluate_assembly(&truth, &pred, output.as_deref()) {
                eprintln!("Error during evaluation: {}", e);
            }
        }
        
        Commands::Visualize { matrix, output, heatmap, pca, components } => {
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
            
            let sample_names: Vec<String> = header.split('\t').skip(1).map(|s| s.to_string()).collect();
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
                match visualize::plot::plot_heatmap_with_labels(&matrix_data, &transcript_ids, &sample_names, &heatmap_path) {
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
        
        Commands::Traverse { input, segments, output, formats, include_edges, visualize, metadata } => {
            info!("Traversing paths in GFA file");
            
            match io::gfa_utils::traverse_paths(&input, &segments, &output, &formats, include_edges, visualize, metadata) {
                Ok(count) => info!("Successfully traversed {} paths", count),
                Err(e) => eprintln!("Error traversing paths: {}", e),
            }
        }
    }
}

