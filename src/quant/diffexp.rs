use std::collections::HashMap;
use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub struct DiffExpResult {
    pub transcript: String,
    pub fold_change: f64,
    pub p_value: f64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DiffExpError {
    EmptyGroup(&'static str),
    MissingSample(String),
}

impl fmt::Display for DiffExpError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DiffExpError::EmptyGroup(group) => {
                write!(f, "differential expression group '{}' is empty", group)
            }
            DiffExpError::MissingSample(sample) => {
                write!(f, "sample '{}' not found in matrix", sample)
            }
        }
    }
}

impl std::error::Error for DiffExpError {}

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
    compute_fold_change_matrix_checked(samples, group_a, group_b).unwrap_or_default()
}

pub fn compute_fold_change_matrix_checked(
    samples: &HashMap<String, Vec<f64>>,
    group_a: &[&str],
    group_b: &[&str],
) -> Result<Vec<DiffExpResult>, DiffExpError> {
    let group_a_values = resolve_group(samples, group_a, "A")?;
    let group_b_values = resolve_group(samples, group_b, "B")?;

    let len = group_a_values
        .iter()
        .chain(group_b_values.iter())
        .map(|values| values.len())
        .min()
        .unwrap_or(0);
    let mut results = Vec::with_capacity(len);
    let mut a_vals = Vec::with_capacity(group_a_values.len());
    let mut b_vals = Vec::with_capacity(group_b_values.len());

    for i in 0..len {
        a_vals.clear();
        b_vals.clear();
        a_vals.extend(
            group_a_values
                .iter()
                .map(|values| finite_or_zero(values[i])),
        );
        b_vals.extend(
            group_b_values
                .iter()
                .map(|values| finite_or_zero(values[i])),
        );

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

    Ok(results)
}

#[inline]
fn resolve_group<'a>(
    samples: &'a HashMap<String, Vec<f64>>,
    group: &[&str],
    group_name: &'static str,
) -> Result<Vec<&'a [f64]>, DiffExpError> {
    if group.is_empty() {
        return Err(DiffExpError::EmptyGroup(group_name));
    }

    let mut resolved = Vec::with_capacity(group.len());
    for &sample in group {
        let values = samples
            .get(sample)
            .ok_or_else(|| DiffExpError::MissingSample(sample.to_string()))?;
        resolved.push(values.as_slice());
    }
    Ok(resolved)
}

#[inline]
fn finite_or_zero(value: f64) -> f64 {
    if value.is_finite() {
        value
    } else {
        0.0
    }
}

fn t_test(a: &[f64], b: &[f64]) -> f64 {
    if a.is_empty() || b.is_empty() {
        return 1.0;
    }

    // Welch's t-test approximation
    let var_a = variance(a);
    let var_b = variance(b);
    let mean_a = a.iter().sum::<f64>() / a.len() as f64;
    let mean_b = b.iter().sum::<f64>() / b.len() as f64;

    let denom = (var_a / a.len() as f64) + (var_b / b.len() as f64);
    if !denom.is_finite() || denom <= f64::EPSILON {
        return 1.0;
    }

    let t = (mean_a - mean_b).abs() / denom.sqrt();
    if !t.is_finite() {
        return 1.0;
    }

    let _dof = a.len().min(b.len()) as f64 - 1.0;
    // Simple approximation of p-value
    // In a real implementation, use statrs crate with Student's T distribution
    (1.0 / (1.0 + t * t)).clamp(0.0, 1.0)
}

fn variance(x: &[f64]) -> f64 {
    if x.is_empty() {
        return 0.0;
    }
    let mean = x.iter().sum::<f64>() / x.len() as f64;
    let var = x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / x.len() as f64;
    if var.is_finite() {
        var
    } else {
        0.0
    }
}

/// Write differential expression results to a file
pub fn write_diffexp(results: &[DiffExpResult], output: &str) -> std::io::Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(output)?;
    writeln!(f, "transcript_id\tlog2FC\tp_value")?;
    for r in results {
        writeln!(
            f,
            "{}\t{:.3}\t{:.4}",
            r.transcript,
            r.fold_change.log2(),
            r.p_value
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use tempfile::NamedTempFile;

    #[test]
    fn compute_fold_change_is_deterministic_across_sample_insertion_orders() {
        let mut samples_a = HashMap::new();
        samples_a.insert("treated".to_string(), vec![5.0, 8.0, 13.0]);
        samples_a.insert("control".to_string(), vec![2.0, 3.0, 5.0]);

        let mut samples_b = HashMap::new();
        samples_b.insert("control".to_string(), vec![2.0, 3.0, 5.0]);
        samples_b.insert("treated".to_string(), vec![5.0, 8.0, 13.0]);

        let result_a =
            compute_fold_change_matrix_checked(&samples_a, &["control"], &["treated"]).unwrap();
        let result_b =
            compute_fold_change_matrix_checked(&samples_b, &["control"], &["treated"]).unwrap();

        assert_eq!(result_a, result_b);
    }

    #[test]
    fn compute_fold_change_checked_rejects_empty_or_missing_groups() {
        let mut samples = HashMap::new();
        samples.insert("s1".to_string(), vec![1.0, 2.0]);

        let err = compute_fold_change_matrix_checked(&samples, &[], &["s1"]).unwrap_err();
        assert_eq!(err, DiffExpError::EmptyGroup("A"));

        let err = compute_fold_change_matrix_checked(&samples, &["s1"], &["missing"]).unwrap_err();
        assert_eq!(err, DiffExpError::MissingSample("missing".to_string()));
    }

    #[test]
    fn compute_fold_change_uses_shortest_sample_length_without_panicking() {
        let mut samples = HashMap::new();
        samples.insert("a".to_string(), vec![1.0, 2.0, 3.0, 4.0]);
        samples.insert("b".to_string(), vec![2.0, 3.0]);

        let results = compute_fold_change_matrix_checked(&samples, &["a"], &["b"]).unwrap();
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].transcript, "transcript_0");
        assert_eq!(results[1].transcript, "transcript_1");
    }

    #[test]
    fn compute_fold_change_sanitizes_non_finite_inputs() {
        let mut samples = HashMap::new();
        samples.insert("a".to_string(), vec![f64::NAN, f64::INFINITY]);
        samples.insert("b".to_string(), vec![1.0, 1.0]);

        let results = compute_fold_change_matrix_checked(&samples, &["a"], &["b"]).unwrap();
        assert_eq!(results.len(), 2);
        for result in results {
            assert!(result.fold_change.is_finite());
            assert!((0.0..=1.0).contains(&result.p_value));
        }
    }

    #[test]
    fn write_diffexp_writes_tsv_output() {
        let file = NamedTempFile::new().unwrap();
        let results = vec![DiffExpResult {
            transcript: "transcript_5".to_string(),
            fold_change: 2.0,
            p_value: 0.125,
        }];

        write_diffexp(&results, file.path().to_str().unwrap()).unwrap();
        let content = std::fs::read_to_string(file.path()).unwrap();
        let lines: Vec<&str> = content.lines().collect();

        assert_eq!(lines[0], "transcript_id\tlog2FC\tp_value");
        assert_eq!(lines[1], "transcript_5\t1.000\t0.1250");
    }
}
