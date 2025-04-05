use smartcore::linalg::basic::matrix::DenseMatrix;
use smartcore::linalg::basic::arrays::Array;
use smartcore::decomposition::pca::*;
use plotters::prelude::*;
use std::collections::HashMap;

pub fn plot_pca(
    samples: &HashMap<String, Vec<f64>>,
    output: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let sample_names: Vec<String> = samples.keys().cloned().collect();
    let data: Vec<Vec<f64>> = sample_names.iter().map(|s| samples[s].clone()).collect();
    let matrix = DenseMatrix::from_2d_vec(&data);

    let pca = PCA::fit(&matrix, PCAParameters::default().with_n_components(2))?;
    let proj = pca.transform(&matrix)?;

    let root = SVGBackend::new(output, (800, 600)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("PCA of TPM Matrix", ("sans-serif", 30))
        .margin(20)
        .x_label_area_size(40)
        .y_label_area_size(40)
        .build_cartesian_2d(-1.0..1.0, -1.0..1.0)?;

    chart.configure_mesh().draw()?;

    let (nrows, _) = proj.shape();
    for i in 0..nrows {
        let x = *proj.get((i, 0));
        let y = *proj.get((i, 1));
        chart.draw_series(std::iter::once(Circle::new(
            (x, y),
            5,
            BLUE.filled(),
        )))?.label(&sample_names[i]);
    }

    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .draw()?;

    Ok(())
} 