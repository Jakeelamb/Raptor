use plotters::prelude::*;
use ndarray::Array2;

/// Plot a heatmap visualization of the TPM matrix
pub fn plot_heatmap(matrix: &Array2<f64>, output: &str) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_val = matrix.iter().cloned().fold(f64::MIN, f64::max);
    let min_val = matrix.iter().cloned().fold(f64::MAX, f64::min);

    let chart = ChartBuilder::on(&root)
        .caption("TPM Heatmap", ("sans-serif", 30))
        .build_cartesian_2d(0..matrix.ncols(), 0..matrix.nrows())?;

    chart.configure_mesh().draw()?;

    for (i, row) in matrix.axis_iter(ndarray::Axis(0)).enumerate() {
        for (j, &val) in row.iter().enumerate() {
            let intensity = (val - min_val) / (max_val - min_val + 1e-6);
            chart.draw_series(std::iter::once(Rectangle::new(
                [(j, i), (j + 1, i + 1)],
                HSLColor(intensity, 0.9, 0.4).filled(),
            )))?;
        }
    }

    Ok(())
}

/// Enhanced heatmap with row and column labels
pub fn plot_heatmap_with_labels(
    matrix: &Array2<f64>, 
    row_labels: &[String], 
    col_labels: &[String], 
    output: &str
) -> Result<(), Box<dyn std::error::Error>> {
    let root = BitMapBackend::new(output, (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;

    let max_val = matrix.iter().cloned().fold(f64::MIN, f64::max);
    let min_val = matrix.iter().cloned().fold(f64::MAX, f64::min);

    let drawing_area = root.margin(40, 40, 60, 120);
    
    let chart = ChartBuilder::on(&drawing_area)
        .caption("TPM Heatmap", ("sans-serif", 30))
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(0..matrix.ncols(), 0..matrix.nrows())?;

    let mut chart = chart
        .configure_mesh()
        .disable_mesh()
        .x_labels(matrix.ncols())
        .y_labels(matrix.nrows())
        .x_label_formatter(&|x| {
            if *x < col_labels.len() {
                col_labels[*x].clone()
            } else {
                String::new()
            }
        })
        .y_label_formatter(&|y| {
            if *y < row_labels.len() {
                row_labels[*y].clone()
            } else {
                String::new()
            }
        })
        .x_label_style(("sans-serif", 8).into_font().transform(FontTransform::Rotate90))
        .draw()?;

    for (i, row) in matrix.axis_iter(ndarray::Axis(0)).enumerate() {
        for (j, &val) in row.iter().enumerate() {
            let intensity = (val - min_val) / (max_val - min_val + 1e-6);
            chart.draw_series(std::iter::once(Rectangle::new(
                [(j, i), (j + 1, i + 1)],
                HSLColor(intensity, 0.9, 0.4).filled(),
            )))?;
        }
    }

    // Add a color scale legend
    let legend_area = root.margin(40, 40, 40, 40);
    let mut legend_chart = ChartBuilder::on(&legend_area)
        .margin(5)
        .set_all_label_area_size(0)
        .build_cartesian_2d(0..100, 0..1)?;

    for i in 0..100 {
        let intensity = i as f64 / 100.0;
        legend_chart.draw_series(std::iter::once(Rectangle::new(
            [(i, 0), (i + 1, 1)],
            HSLColor(intensity, 0.9, 0.4).filled(),
        )))?;
    }

    // Add min/max labels
    root.draw(&Text::new(
        format!("{:.1}", min_val),
        (45, 30),
        ("sans-serif", 10).into_font(),
    ))?;
    
    root.draw(&Text::new(
        format!("{:.1}", max_val),
        (135, 30),
        ("sans-serif", 10).into_font(),
    ))?;

    Ok(())
} 