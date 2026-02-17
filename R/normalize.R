# Normalization functions

#' Normalize and Scale Data for Statistical Analysis
#'
#' Applies a scaling method (i.e. total sum normalization or autoscaling) to metabolite values.
#'
#' @param corrected_data A data frame of corrected metabolite values as returned by `qc_drift_correction()`.
#' @param method The scaling method: "tsn" (Total Sum Normalization, default) or "auto" (Autoscaling).
#' 
#' @return The scaled data frame.
#' @importFrom dplyr left_join select
#' @importFrom stats sd
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom tibble column_to_rownames rownames_to_column
#' @export
normalize_data <- function(corrected_data, method = "tsn") {

  # --- Reshape to Wide Format for Scaling ---
  corrected_wide <- corrected_data |>
    tidyr::pivot_wider(id_cols = metabolites, names_from = sample, values_from = corrected_abundance) |>
    tibble::column_to_rownames("metabolites") |>
    as.matrix()
  
  # --- Apply Scaling Method ---
  if (method == "auto") {
    # Autoscale (Z-score Scaling)
    data_t <- t(corrected_wide)
    scaled_data <- scale(data_t, center = TRUE, scale = TRUE)

  } else if (method == "tsn") {
    # Total Sum Normalization (TSN)

    # 1. Calculate Sample Sums (Samples are in columns)
    col_sums <- apply(corrected_wide, 2, sum, na.rm = TRUE)
    eps <- .Machine$double.eps # Guard against zero sums
    col_sums <- pmax(col_sums, eps)

    # 2. Divide by Sample Sums (Normalization for dilution)
    norm_by_sum <- sweep(corrected_wide, 2, col_sums, "/")

    # 3. Rescale by Metabolite Averages (Row Means)
    row_means <- rowMeans(corrected_wide, na.rm = TRUE)
    scaled_data <- sweep(norm_by_sum, 1, row_means, FUN = "*")

  } else {
    stop("Invalid scaling method specified. Use 'tsn' or 'auto'.")
  }

  # --- Reshape Back to Long Format and Join into `corrected_data` ---
  long_scaled <- scaled_data |>
    as.data.frame() |>
    tibble::rownames_to_column("metabolites") |>
    tidyr::pivot_longer(cols = -metabolites, names_to = "sample", values_to = "scaled_abundance") |>
    dplyr::right_join(corrected_data, by = c("metabolites", "sample"))

  return(long_scaled)
}
