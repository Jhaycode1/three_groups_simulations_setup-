# Check calibration of posterior probability of being non-null.
# A non-null probability is 1 - posterior_null.

# Compare predicted non-null probabilities with observed non-null frequencies.
calibration_metrics <- function(comparison, breaks = seq(0, 1, by = 0.1)) {
  required_columns <- c("method", "posterior_null", "true_group")
  missing_columns <- setdiff(required_columns, names(comparison))
  if (length(missing_columns) > 0) {
    stop("Comparison is missing: ", paste(missing_columns, collapse = ", "))
  }
  if (length(breaks) < 2 || any(diff(breaks) <= 0) ||
      breaks[1] < 0 || breaks[length(breaks)] > 1) {
    stop("breaks must be increasing values within 0 and 1.")
  }

  result_rows <- list()
  row_number <- 1L
  for (method_name in unique(comparison$method)) {
    method_data <- comparison[comparison$method == method_name, , drop = FALSE]
    # Calibration is for P(non-null), not P(null).
    predicted <- 1 - method_data$posterior_null
    observed <- as.numeric(method_data$true_group != "null")
    # Group predictions into fixed probability intervals.
    bin_id <- cut(
      predicted,
      breaks = breaks,
      include.lowest = TRUE,
      labels = FALSE
    )

    counts <- tabulate(bin_id, nbins = length(breaks) - 1L)
    for (bin_number in seq_len(length(breaks) - 1L)) {
      in_bin <- bin_id == bin_number
      number_in_bin <- counts[bin_number]
      # Empty bins are retained and reported as NA rather than invented values.
      result_rows[[row_number]] <- data.frame(
        method = method_name,
        bin_lower = breaks[bin_number],
        bin_upper = breaks[bin_number + 1L],
        number_in_bin = number_in_bin,
        mean_predicted_nonnull = if (number_in_bin == 0) NA_real_ else
          mean(predicted[in_bin]),
        observed_nonnull_rate = if (number_in_bin == 0) NA_real_ else
          mean(observed[in_bin]),
        absolute_calibration_error = if (number_in_bin == 0) NA_real_ else
          abs(mean(predicted[in_bin]) - mean(observed[in_bin]))
      )
      row_number <- row_number + 1L
    }
  }

  result <- do.call(rbind, result_rows)
  result$weighted_calibration_error <- NA_real_
  for (method_name in unique(result$method)) {
    method_rows <- result$method == method_name
    total <- sum(result$number_in_bin[method_rows])
    result$weighted_calibration_error[method_rows] <- if (total == 0) {
      NA_real_
    } else {
      sum(
        result$number_in_bin[method_rows] *
          ifelse(
            is.na(result$absolute_calibration_error[method_rows]),
            0,
            result$absolute_calibration_error[method_rows]
          )
      ) / total
    }
  }

  result
}

# Example command:
# source("truth_comparison.R")
# source("calibration_calc.R")
# posterior_summary <- readRDS(
#   "test_run_results/sample_test_1_posterior_summary.rds"
# )
# metadata <- readRDS("test_run_results/sample_test_1_metadata.rds")
# comparison <- compare_posterior_to_truth(
#   posterior_summary,
#   metadata = metadata
# )
# calibration <- calibration_metrics(comparison)
# View(calibration)
