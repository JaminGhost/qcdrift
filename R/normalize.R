# Normalization functions

#' Normalize and Scale Data for Statistical Analysis
#'
#' Applies a scaling method (e.g., Total Sum Normalization) to metabolite values.
#'
#' @param corrected_data matrix of QC-corrected mass spec values with metabolites on the rows and one column per sample.
#' @param method The scaling method: "tsn" (Total Sum Normalization, default) or "auto" (Autoscaling).
#' @return The scaled data frame.
#' @importFrom stats sd
#' @export
normalize_data <- function(corrected_data, method = "tsn") {

  # --- Apply Scaling Method ---
  if (method == "auto") {
    # Autoscale (Z-score Scaling)
    data_t <- t(corrected_data)
    scaled_data_t <- scale(data_t, center = TRUE, scale = TRUE)

  } else if (method == "tsn") {
    # Total Sum Normalization (TSN)

    # 1. Calculate Sample Sums (Samples are in columns)
    col_sums <- apply(corrected_data, 2, sum, na.rm = TRUE)
    eps <- .Machine$double.eps # Guard against zero sums
    col_sums <- pmax(col_sums, eps)

    # 2. Divide by Sample Sums (Normalization for dilution)
    norm_by_sum <- sweep(corrected_data, 2, col_sums, "/")

    # 3. Rescale by Metabolite Averages (Row Means)
    row_means <- rowMeans(corrected_data, na.rm = TRUE)
    scaled_data <- sweep(norm_by_sum, 1, row_means, FUN = "*")

  } else {
    stop("Invalid scaling method specified. Use 'tsn' or 'auto'.")
  }

  return(scaled_data)
}
