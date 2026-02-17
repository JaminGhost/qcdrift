# functions for correcting qc drift

#' Linearly Correct for QC Drift
#'
#' @description This function corrects for analytical drift in mass spectrometry data
#' by linearly interpolating between QC samples. It assumes that the first and
#' last samples are QC samples.
#'
#' @param io Numeric vector indicating injection order.
#' @param abundance Numeric vector containing uncorrected abundances.
#' @param qc Logical vector identifying QC samples.
#'
#' @return A numeric vector of QC-corrected abundances.
#' @export
correct_linear <- function(io, abundance, qc)
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
