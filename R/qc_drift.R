# functions for correcting qc drift

#' Correct for QC Drift
#' @rdname qc_drift_correction
#'
#' @param raw_data A data frame containing the raw mass spectrometry data with columns for injection order (`io`), abundance (`abundance`), and a logical column indicating QC samples (`qc`).
#' @param correction_function A function to apply for drift correction. Currently, the only supported function is `correct_linear`, which performs linear interpolation between QC samples.
#' @param io A vector indicating injection order
#' @param abundance A vector of uncorrected abundances
#' @param qc A logical vector indicating QC samples
#' @param metabolites A vector identifying each metabolite
#'
#' @description These functions apply a correction for analytical drift in mass spectrometry data.
#' They assume that `raw_data` has appropriate columns for the chosen correction function.
#' 
#' `correct_linear` corrects for analytical drift in mass spectrometry data
#' by linearly interpolating between QC samples. It assumes that the first and
#' last samples are QC samples.
#' 
#' @return `qc_drift_correction` returns a data frame with an additional column, `corrected_abundance`, containing the QC-corrected abundances returned by `correction_function`.
#' 
#' @export
#' @importFrom dplyr group_by mutate ungroup
qc_drift_correction <- function(raw_data, correction_function = correct_linear, ...)
{
  corrected_data <- group_by(raw_data, metabolites) |>
    mutate(corrected_abundance = correction_function(io, abundance, qc)) |>
    ungroup()

  return(corrected_data)
}

#' @rdname qc_drift_correction
#'
#' @export
correct_linear <- function(io, abundance, qc, metabolites)
{
  # assume io, abundance, and qc are of equal length from a data.frame
  if(length(io) != length(abundance) | length(abundance) != length(qc))
    stop("io, abundance and qc must all be vectors of the same length.")

  # check that there are qc values to use
  if(sum(qc) < 2)
  {
    warning("Insufficient number of QC values - skipping correction")
    return(abundance)
  }

  # this code currently assumes the first injection is a QC sample
  if(!qc[1])
    stop("Code assumes the first injection is a QC sample.")
  if(!qc[length(qc)])
    stop("Code assumes the last injection is a QC sample.")

  # vector of qc-corrected values
  retval <- rep(NA, length(abundance))

  # Loop over pairs of QC samples
  qc_indeces <- (1:length(qc))[qc]

  # QC baseline
  retval[1] <- mean(abundance[qc_indeces], na.rm = TRUE)

  for(i in 2:length(qc_indeces))
  {
    a <- qc_indeces[i-1]
    b <- qc_indeces[i]

    # formula for this segment of the sequence
    correction_fct <- function(x)
    {
      denom <- (io[b] - io[a])
      if(denom == 0) return(rep(1, length(x))) # Safety check for zero-distance

      # slope
      m <- (abundance[b] - abundance[a]) / denom

      # correction factor = baseline / point slope
      # Safety check to prevent division by zero if drift hits 0
      pred_val <- (m * (x - io[a]) + abundance[a])
      pred_val[pred_val == 0] <- 1e-9

      retval[1] / pred_val
    }

    # correct values
    retval[(a+1):b] <- abundance[(a+1):b] * correction_fct(io[(a+1):b])
  }

  return(retval)
}
