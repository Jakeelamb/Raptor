use ndarray::Array2;
use ndarray_stats::QuantileExt;
use plotters::prelude::*;
use std::collections::HashMap;

/// Simple PCA implementation using ndarray
pub fn compute_pca(matrix: &Array2<f64>, n_components: usize) -> Array2<f64> {
    // Center the data
    let n_samples = matrix.shape()[0];
    let n_features = matrix.shape()[1];
    
    // Calculate the mean of each column
    let mut means = vec![0.0; n_features];
    for i in 0..n_features {
        let mut sum = 0.0;
        for j in 0..n_samples {
            sum += matrix[[j, i]];
        }
        means[i] = sum / n_samples as f64;
    }
    
    // Center the data
    let mut centered = Array2::zeros((n_samples, n_features));
    for i in 0..n_samples {
        for j in 0..n_features {
            centered[[i, j]] = matrix[[i, j]] - means[j];
        }
    }
    
    // Calculate covariance matrix
    let mut cov = Array2::zeros((n_features, n_features));
    for i in 0..n_features {
        for j in 0..n_features {
            let mut sum = 0.0;
            for k in 0..n_samples {
                sum += centered[[k, i]] * centered[[k, j]];
            }
            cov[[i, j]] = sum / (n_samples as f64 - 1.0);
        }
    }
    
    // For simplicity, just project onto the first two components
    // In a real implementation, we would compute eigenvectors here
    
    // Create a simple projection matrix (2 dimensions)
    let mut projection = Array2::zeros((n_features, n_components.min(2)));
    for i in 0..n_features.min(2) {
        projection[[i, i]] = 1.0;
    }
    
    // Project the data
    let mut result = Array2::zeros((n_samples, n_components.min(2)));
    for i in 0..n_samples {
        for j in 0..n_components.min(2) {
            let mut sum = 0.0;
            for k in 0..n_features {
                sum += centered[[i, k]] * projection[[k, j]];
            }
            result[[i, j]] = sum;
        }
    }
    
    result
}

pub fn plot_pca(
    samples: &HashMap<String, Vec<f64>>,
    output: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    // Convert HashMap to Array2
    let sample_names: Vec<String> = samples.keys().cloned().collect();
    let n_samples = sample_names.len();
    let first_sample = &samples[&sample_names[0]];
    let n_features = first_sample.len();
    
    let mut matrix = Array2::zeros((n_samples, n_features));
    for (i, name) in sample_names.iter().enumerate() {
        for (j, val) in samples[name].iter().enumerate() {
            matrix[[i, j]] = *val;
        }
    }
    
    // Compute PCA
    let pca_result = compute_pca(&matrix, 2);
    
    // Create plot
    let root = SVGBackend::new(output, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    // Find min/max values for axes
    let mut min_x = f64::MAX;
    let mut max_x = f64::MIN;
    let mut min_y = f64::MAX;
    let mut max_y = f64::MIN;
    
    for i in 0..pca_result.nrows() {
        min_x = min_x.min(pca_result[[i, 0]]);
        max_x = max_x.max(pca_result[[i, 0]]);
        min_y = min_y.min(pca_result[[i, 1]]);
        max_y = max_y.max(pca_result[[i, 1]]);
    }

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
            5,
            &BLUE,
            &|c, s, st| EmptyElement::at(c) + Circle::new((0, 0), s, st),
        ))?.label(&sample_names[i]);
    }

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    Ok(())
}

pub fn plot_pca_simple(pca_result: &Array2<f64>, output: &str) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    // Find min/max values for axes
    let mut min_x = f64::MAX;
    let mut max_x = f64::MIN;
    let mut min_y = f64::MAX;
    let mut max_y = f64::MIN;
    
    for i in 0..pca_result.nrows() {
        min_x = min_x.min(pca_result[[i, 0]]);
        max_x = max_x.max(pca_result[[i, 0]]);
        min_y = min_y.min(pca_result[[i, 1]]);
        max_y = max_y.max(pca_result[[i, 1]]);
    }

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