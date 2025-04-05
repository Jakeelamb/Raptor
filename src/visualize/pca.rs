use ndarray::{Array2, Axis};
use ndarray_stats::QuantileExt;
use linfa::prelude::*;
use linfa_reduction::Pca;
use plotters::prelude::*;
use std::collections::HashMap;

pub fn compute_pca(matrix: &Array2<f64>, n: usize) -> Array2<f64> {
    let pca = Pca::params(n).fit(matrix).expect("PCA failed");
    pca.transform(matrix)
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

    let (min_x, max_x) = pca_result.column(0).minmax().into_option().unwrap();
    let (min_y, max_y) = pca_result.column(1).minmax().into_option().unwrap();

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

    let (min_x, max_x) = pca_result.column(0).minmax().into_option().unwrap();
    let (min_y, max_y) = pca_result.column(1).minmax().into_option().unwrap();

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