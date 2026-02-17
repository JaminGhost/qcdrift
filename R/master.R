#' Process Mass Spectrometry Data and Generate QC Plots
#'
#' This function performs a full analysis pipeline on mass spectrometry data,
#' including reading and cleaning the data, drift correction, normalization,
#' and generation of various quality control plots.
#'
#' @param filepath The path to the Excel file containing the raw mass spectrometry data.
#' @param sheet_name The name or index of the sheet to process within the Excel file. Default is 1.
#' @param correction_function The function to use for drift correction. Default is `correct_linear`, which performs linear interpolation between QC samples.
#' @param norm_method The method to use for normalization. Default is 'tsn' (Total Sum Normalization), but can also be set to 'auto' for autoscaling (Z-score normalization).
#'
#' @return A list containing two main elements:
#'   \item{data}{A list of data frames: `raw` (the initial cleaned data),
#'     `corrected` (data after drift correction), and `scaled` (data after normalization).}
#'   \item{plots}{A list of ggplot2 objects for quality control analysis,
#'     including PCA plots, CV plots, violin plots, and heatmaps.}
#'
#' @importFrom dplyr group_by mutate ungroup summarise left_join
#' @importFrom patchwork plot_annotation
#'
#' @export
process_runs <- function(filepath, sheet_name=1, correction_function = correct_linear,
                         norm_method = 'tsn') {

  # 1. READ & CLEAN
  long_data <- read_and_clean_data(filepath, sheet_name)

  # 2. DRIFT CORRECTION
  long_corrected <- qc_drift_correction(long_data)

  # 3. NORMALIZE (TSN)
  sample_sums <- long_corrected |>
    dplyr::group_by(sample) |>
    dplyr::summarise(total_sum = sum(corrected_abundance, na.rm=TRUE), .groups = "drop")

  long_scaled <- long_corrected |>
    dplyr::left_join(sample_sums, by="sample") |>
    dplyr::group_by(metabolites) |>
    dplyr::mutate(
      norm = corrected_abundance / total_sum,
      scaled_abundance = norm * mean(corrected_abundance, na.rm=TRUE)
    ) |>
    dplyr::ungroup()

  # 4. PLOT GENERATION
  pca_raw    <- run_pca(long_data, abundance)
  pca_scaled <- run_pca(long_scaled, scaled_abundance)

  # Heatmaps
  p_heat_raw  <- generate_qc_heatmap(long_data, abundance, title_prefix="Raw QC Heatmap")
  p_heat_corr <- generate_qc_heatmap(long_scaled, scaled_abundance, title_prefix="Corrected QC Heatmap")

  return(list(
    data = list(raw = long_data, corrected = long_corrected, scaled = long_scaled),
    plots = list(
      pca_raw = generate_pca_plot(pca_raw, title_prefix="Raw Data"),
      pca_corrected = generate_pca_plot(pca_scaled, title_prefix="Corrected & Scaled"),
      pca_pairs = generate_pca_pairs(pca_scaled, npcs = 4),
      cv_raw = generate_cv_plot(long_data, abundance, title_prefix="CV (Raw)"),
      cv_corrected = generate_cv_plot(long_scaled, scaled_abundance, title_prefix="CV (Corrected)"),
      violin_raw = generate_qc_violin(long_data, abundance, title_prefix="Raw QC Dist"),
      violin_corrected = generate_qc_violin(long_scaled, scaled_abundance, title_prefix="Corrected QC Dist"),
      heatmap_raw = p_heat_raw,
      heatmap_corr = p_heat_corr,
      # Combined Heatmap Page
      heatmap_comparison = p_heat_raw + p_heat_corr +
        patchwork::plot_annotation(title = "Heatmap: Metabolite Intensities (Z-Score)")
    )
  ))
}
