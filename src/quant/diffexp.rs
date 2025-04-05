use std::collections::HashMap;

pub struct DiffExpResult {
    pub transcript: String,
    pub fold_change: f64,
    pub p_value: f64,
}

/// Trait extension to add log2 method to f64
pub trait Log2 {
    fn log2(self) -> f64;
}

impl Log2 for f64 {
    fn log2(self) -> f64 {
        self.ln() / 2f64.ln()
    }
}

pub fn compute_fold_change_matrix(
    samples: &HashMap<String, Vec<f64>>,
    group_a: &[&str],
    group_b: &[&str],
) -> Vec<DiffExpResult> {
    let len = samples.values().next().unwrap().len();
    let mut results = vec![];

    for i in 0..len {
        let a_vals: Vec<f64> = group_a.iter().map(|&s| samples[s][i]).collect();
        let b_vals: Vec<f64> = group_b.iter().map(|&s| samples[s][i]).collect();

        let mean_a = a_vals.iter().sum::<f64>() / a_vals.len() as f64;
        let mean_b = b_vals.iter().sum::<f64>() / b_vals.len() as f64;
        let fold_change = (mean_b + 1.0) / (mean_a + 1.0); // add pseudo-count
        let p_value = t_test(&a_vals, &b_vals);

        results.push(DiffExpResult {
            transcript: format!("transcript_{}", i),
            fold_change,
            p_value,
        });
    }

    results
}

fn t_test(a: &[f64], b: &[f64]) -> f64 {
    // Welch's t-test approximation
    let var_a = variance(a);
    let var_b = variance(b);
    let mean_a = a.iter().sum::<f64>() / a.len() as f64;
    let mean_b = b.iter().sum::<f64>() / b.len() as f64;

    let t = (mean_a - mean_b).abs()
        / ((var_a / a.len() as f64) + (var_b / b.len() as f64)).sqrt();

    let dof = a.len().min(b.len()) as f64 - 1.0;
    // Simple approximation of p-value
    // In a real implementation, use statrs crate with Student's T distribution
    1.0 / (1.0 + t * t)
}

fn variance(x: &[f64]) -> f64 {
    let mean = x.iter().sum::<f64>() / x.len() as f64;
    x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / x.len() as f64
}

/// Write differential expression results to a file
pub fn write_diffexp(results: &[DiffExpResult], output: &str) {
    use std::io::Write;
    let mut f = std::fs::File::create(output).unwrap();
    writeln!(f, "transcript_id\tlog2FC\tp_value").unwrap();
    for r in results {
        writeln!(f, "{}\t{:.3}\t{:.4}", r.transcript, r.fold_change.log2(), r.p_value).unwrap();
    }
} 