library(testthat)
library(readxl)

# Tests for qc_interpolate()
test_that("qc_interpolate calculates correct value", {
  expect_equal(qc_interpolate(3, 1, 5, 100, 200), 150)
  expect_equal(qc_interpolate(1, 1, 5, 100, 200), 100)
  expect_equal(qc_interpolate(5, 1, 5, 100, 200), 200)
})

# Tests for normalize_data()
test_that("normalize_data works correctly", {
  corrected_data <- data.frame(
    Metabolite = c("A", "B"),
    Sample1 = c(10, 20),
    Sample2 = c(30, 40)
  )
  
  # Test TSN normalization
  tsn_data <- normalize_data(corrected_data, method = "tsn")
  expect_equal(ncol(tsn_data), 3)
  expect_equal(nrow(tsn_data), 2)
  expect_true("Metabolite" %in% names(tsn_data))
  
  # Test auto scaling (z-score)
  auto_data <- normalize_data(corrected_data, method = "auto")
  expect_equal(ncol(auto_data), 3)
  expect_equal(nrow(auto_data), 2)
  expect_true("Metabolite" %in% names(auto_data))
})

# Tests for calculate_stats()
test_that("calculate_stats works correctly", {
  scaled_data <- data.frame(
    Metabolite = c("A", "B"),
    Sample1 = c(10, 20),
    Sample2 = c(30, 40)
  )
  
  stats <- calculate_stats(scaled_data)
  
  expect_true("stats_summary" %in% names(stats))
  expect_true("z_scores_table" %in% names(stats))
  expect_equal(nrow(stats$stats_summary), 2)
  expect_equal(nrow(stats$z_scores_table), 2)
})

# Integration test for process_runs()
test_that("process_runs executes without error", {
  # Path to the test data file included in the package
  filepath <- system.file("extdata", "RawMassSpec.xlsx", package = "msdata")
  
  # Check if the file exists
  if (file.exists(filepath)) {
    # Run the processing function
    result <- process_runs(filepath, sheet_name = "Sheet1")
    
    # Check the structure of the output
    expect_true("corrected_data" %in% names(result))
    expect_true("stats_summary" %in% names(result))
    expect_true("z_scores_table" %in% names(result))
    expect_true("raw_qc_plot" %in% names(result))
    expect_true("corrected_qc_plot" %in% names(result))
    expect_true("pca_biplot" %in% names(result))
    
    # Check that the output tables are not empty
    expect_gt(nrow(result$corrected_data), 0)
    expect_gt(nrow(result$stats_summary), 0)
    expect_gt(nrow(result$z_scores_table), 0)
    
    # Check that the function returns a list
    expect_true(is.list(result))
  } else {
    skip("Test data file not found.")
  }
})

test_that("process_runs returns a list", {
  filepath <- system.file("extdata", "RawMassSpec.xlsx", package = "msdata")
  if (file.exists(filepath)) {
    result <- process_runs(filepath, sheet_name = "Sheet1")
    expect_type(result, type="list")
  } else {
    skip("Test data file not found.")
  }
})
