# utils.R

#' @importFrom utils globalVariables
NULL

# Declares all internal column names and dynamically created variables
# this resolves "no visible binding for global variable" warnings
utils::globalVariables(c(

  # variables used in plotting functions
  "qc", "metabolites", "sd", "ZScore", "log_abund",

  # variables used in process_runs
  "io", "abundance", "corrected_abundance", "total_sum",

  # Internal variables needed for dplyr/tidyr/ggplot2
  ".", ".data", "CV", "Mean", "SD",

  # Column names from the input data:
  "Sample", "scaled_abundance"
))