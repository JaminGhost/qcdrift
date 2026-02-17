#' Read and Clean Mass Spectrometry Data
#'
#' @description This function reads a sheet from an Excel file containing mass spectrometry data,
#' handles duplicate metabolite names, and transforms the data into a long format.
#' It also adds injection order and identifies QC samples.
#'
#' @param filepath Path to the Excel file.
#' @param sheet_name The name or index of the sheet to read. Defaults to 1.
#' @param qc_starts_with A character string prefix to identify QC samples. Defaults to 'QC'.
#'
#' @return A tibble in long format with columns: `metabolites`, `sample`, `abundance`, `io` (injection order), and `qc` (a boolean indicating if the sample is a QC sample).
#'
#' @importFrom dplyr select
#' @importFrom stats var
#' @importFrom tidyr pivot_longer
#' @export
read_and_clean_data <- function(filepath, sheet_name = 1, qc_starts_with = 'QC')
{

  # --- STEP 1: READ DATA and HANDLE DUPLICATES (FINAL FIXED ALIGNMENT) ---
  # Read the entire sheet without header assumption
  full_data <- readxl::read_excel(filepath, sheet = sheet_name)

  # 1. Extract Sample Names and Injection Orders
  # These should exclude the first column (Metabolite Name).
  injection_orders <- as.numeric(full_data[1, -1])
  names(injection_orders) <- names(full_data)[-1]

  # drop injection order row
  full_data <- full_data[-1,]

  # 2. Process metabolite name column
  names(full_data)[1] <- 'metabolites'

  # Create UNIQUE metabolite names (e.g. "Alanine_1", "Alanine_2")
  tmp <- table(full_data$metabolites)
  tmp <- tmp[tmp > 1]

  if(length(tmp) > 0)
  {
    for(i in 1:length(tmp))
    {
      full_data$metabolites[full_data$metabolites == names(tmp)[i]] <- paste(names(tmp)[i], 1:tmp, sep = '_')
    }
  }

  # 3. Convert to long format
  full_data <- pivot_longer(full_data, -1, names_to = 'sample', values_to = 'abundance')

  # 4. Add injection order and identify QC samples
  full_data <- full_data|>
    mutate(io = injection_orders[sample],
           qc = grepl(paste0('^', qc_starts_with), sample))
  
  return(full_data)
}