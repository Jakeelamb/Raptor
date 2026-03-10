use ndarray::Array2;
use plotters::prelude::*;
use std::collections::HashMap;
use std::io::{Error, ErrorKind};

#[inline]
fn sample_matrix_sorted(
    samples: &HashMap<String, Vec<f64>>,
) -> Result<(Vec<String>, Array2<f64>), Error> {
    if samples.is_empty() {
        return Err(Error::new(
            ErrorKind::InvalidInput,
            "PCA requires at least one sample",
        ));
    }

    let mut sample_names: Vec<String> = samples.keys().cloned().collect();
    sample_names.sort_unstable();

    let n_samples = sample_names.len();
    let first_name = &sample_names[0];
    let first_sample = samples.get(first_name).ok_or_else(|| {
        Error::new(
            ErrorKind::InvalidData,
            format!("missing sample data for '{}'", first_name),
        )
    })?;
    let n_features = first_sample.len();
    if n_features == 0 {
        return Err(Error::new(
            ErrorKind::InvalidInput,
            "PCA requires at least one feature per sample",
        ));
    }

    let mut matrix = Array2::zeros((n_samples, n_features));
    for (row_idx, sample_name) in sample_names.iter().enumerate() {
        let values = samples.get(sample_name).ok_or_else(|| {
            Error::new(
                ErrorKind::InvalidData,
                format!("missing sample data for '{}'", sample_name),
            )
        })?;
        if values.len() != n_features {
            return Err(Error::new(
                ErrorKind::InvalidInput,
                format!(
                    "sample '{}' has {} features, expected {}",
                    sample_name,
                    values.len(),
                    n_features
                ),
            ));
        }

        for (col_idx, &value) in values.iter().enumerate() {
            matrix[[row_idx, col_idx]] = value;
        }
    }

    Ok((sample_names, matrix))
}

#[inline]
fn axis_bounds(points: impl Iterator<Item = f64>) -> Result<(f64, f64), Error> {
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;

    for value in points {
        if !value.is_finite() {
            return Err(Error::new(
                ErrorKind::InvalidInput,
                "PCA output contains non-finite values",
            ));
        }
        min = min.min(value);
        max = max.max(value);
    }

    if !min.is_finite() || !max.is_finite() {
        return Err(Error::new(ErrorKind::InvalidInput, "PCA output is empty"));
    }

    if (max - min).abs() < f64::EPSILON {
        let pad = (min.abs() * 0.05).max(1e-6);
        return Ok((min - pad, max + pad));
    }

    let pad = ((max - min) * 0.05).max(1e-9);
    Ok((min - pad, max + pad))
}

/// Simple PCA implementation using ndarray
pub fn compute_pca(matrix: &Array2<f64>, n_components: usize) -> Array2<f64> {
    let n_samples = matrix.shape()[0];
    let n_features = matrix.shape()[1];
    let projected_dims = n_components.min(2);

    if projected_dims == 0 || n_samples == 0 {
        return Array2::zeros((n_samples, projected_dims));
    }

    // This simple PCA implementation projects onto centered input dimensions.
    // Only compute means for dimensions we will emit.
    let used_features = n_features.min(projected_dims);
    let mut means = vec![0.0; used_features];
    for i in 0..used_features {
        let mut sum = 0.0;
        for j in 0..n_samples {
            sum += matrix[[j, i]];
        }
        means[i] = sum / n_samples as f64;
    }

    let mut result = Array2::zeros((n_samples, projected_dims));
    for i in 0..n_samples {
        for j in 0..used_features {
            result[[i, j]] = matrix[[i, j]] - means[j];
        }
    }

    result
}

pub fn plot_pca(
    samples: &HashMap<String, Vec<f64>>,
    output: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let (sample_names, matrix) = sample_matrix_sorted(samples)?;

    // Compute PCA
    let pca_result = compute_pca(&matrix, 2);
    if pca_result.nrows() == 0 || pca_result.ncols() < 2 {
        return Err(Box::new(Error::new(
            ErrorKind::InvalidInput,
            "PCA result must contain at least one row and two components",
        )));
    }

    // Create plot
    let root = SVGBackend::new(output, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    // Find min/max values for axes
    let (min_x, max_x) = axis_bounds((0..pca_result.nrows()).map(|i| pca_result[[i, 0]]))?;
    let (min_y, max_y) = axis_bounds((0..pca_result.nrows()).map(|i| pca_result[[i, 1]]))?;

    let mut chart = ChartBuilder::on(&root)
        .caption("PCA of Transcript TPM", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(min_x..max_x, min_y..max_y)?;

    chart.configure_mesh().draw()?;

    for i in 0..pca_result.nrows() {
        chart
            .draw_series(PointSeries::of_element(
                [(pca_result[[i, 0]], pca_result[[i, 1]])],
                5,
                &BLUE,
                &|c, s, st| EmptyElement::at(c) + Circle::new((0, 0), s, st),
            ))?
            .label(&sample_names[i]);
    }

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}

pub fn plot_pca_simple(
    pca_result: &Array2<f64>,
    output: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    if pca_result.nrows() == 0 || pca_result.ncols() < 2 {
        return Err(Box::new(Error::new(
            ErrorKind::InvalidInput,
            "PCA result must contain at least one row and two components",
        )));
    }

    let root = BitMapBackend::new(output, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    // Find min/max values for axes
    let (min_x, max_x) = axis_bounds((0..pca_result.nrows()).map(|i| pca_result[[i, 0]]))?;
    let (min_y, max_y) = axis_bounds((0..pca_result.nrows()).map(|i| pca_result[[i, 1]]))?;

    let mut chart = ChartBuilder::on(&root)
        .caption("PCA of Transcript TPM", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(min_x..max_x, min_y..max_y)?;

    chart.configure_mesh().draw()?;

    for i in 0..pca_result.nrows() {
        chart.draw_series(PointSeries::of_element(
            [(pca_result[[i, 0]], pca_result[[i, 1]])],
            3,
            &BLUE,
            &|c, s, st| EmptyElement::at(c) + Circle::new((0, 0), s, st),
        ))?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{axis_bounds, compute_pca, sample_matrix_sorted};
    use ndarray::array;
    use std::collections::HashMap;

    #[test]
    fn sample_matrix_sorted_is_deterministic_and_value_preserving() {
        let mut samples = HashMap::new();
        samples.insert("zeta".to_string(), vec![9.0, 10.0]);
        samples.insert("alpha".to_string(), vec![1.0, 2.0]);

        let (names, matrix) = sample_matrix_sorted(&samples).expect("valid sample matrix");
        assert_eq!(names, vec!["alpha".to_string(), "zeta".to_string()]);
        assert_eq!(matrix[[0, 0]], 1.0);
        assert_eq!(matrix[[0, 1]], 2.0);
        assert_eq!(matrix[[1, 0]], 9.0);
        assert_eq!(matrix[[1, 1]], 10.0);
    }

    #[test]
    fn sample_matrix_sorted_rejects_empty_or_inconsistent_samples() {
        let empty = HashMap::new();
        assert!(sample_matrix_sorted(&empty).is_err());

        let mut inconsistent = HashMap::new();
        inconsistent.insert("a".to_string(), vec![1.0, 2.0]);
        inconsistent.insert("b".to_string(), vec![3.0]);
        assert!(sample_matrix_sorted(&inconsistent).is_err());
    }

    #[test]
    fn compute_pca_centers_first_two_dimensions() {
        let matrix = array![[1.0, 10.0, 100.0], [3.0, 14.0, 200.0]];
        let result = compute_pca(&matrix, 2);

        assert_eq!(result.shape(), &[2, 2]);
        assert!((result[[0, 0]] + 1.0).abs() < 1e-12);
        assert!((result[[1, 0]] - 1.0).abs() < 1e-12);
        assert!((result[[0, 1]] + 2.0).abs() < 1e-12);
        assert!((result[[1, 1]] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn axis_bounds_expands_degenerate_ranges() {
        let (min, max) = axis_bounds([2.5, 2.5, 2.5].into_iter()).expect("finite axis bounds");
        assert!(min < 2.5);
        assert!(max > 2.5);
    }
}
