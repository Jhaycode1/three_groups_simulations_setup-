# Calculate posterior-threshold power and false discovery rate.
# A gene is selected when posterior_null <= a cutoff.

# Select genes by posterior-null probability and calculate non-null power/FDR.
power_fdr_metrics <- function(comparison, cutoffs = c(0.05, 0.10, 0.20, 0.50)) {
  # This analysis treats beneficial and deleterious genes together as non-null.
  required_columns <- c("method", "gene", "posterior_null", "true_group")
  missing_columns <- setdiff(required_columns, names(comparison))
  if (length(missing_columns) > 0) {
    stop("Comparison is missing: ", paste(missing_columns, collapse = ", "))
  }
  if (!is.numeric(cutoffs) || any(cutoffs < 0 | cutoffs > 1)) {
    stop("cutoffs must be numeric values between 0 and 1.")
  }

  result_rows <- list()
  row_number <- 1L
  for (method_name in unique(comparison$method)) {
    method_data <- comparison[comparison$method == method_name, , drop = FALSE]
    true_nonnull <- method_data$true_group != "null"
    number_true_nonnull <- sum(true_nonnull)

    for (cutoff in cutoffs) {
      # Small posterior-null probability means evidence against the null group.
      selected <- method_data$posterior_null <= cutoff
      number_selected <- sum(selected)
      number_true_selected <- sum(selected & true_nonnull)
      number_false_selected <- sum(selected & !true_nonnull)

      result_rows[[row_number]] <- data.frame(
        method = method_name,
        cutoff = cutoff,
        number_selected = number_selected,
        number_true_selected = number_true_selected,
        number_false_selected = number_false_selected,
        power = if (number_true_nonnull == 0) NA_real_ else
          number_true_selected / number_true_nonnull,
        fdr = if (number_selected == 0) 0 else
          number_false_selected / number_selected
      )
      row_number <- row_number + 1L
    }
  }

  do.call(rbind, result_rows)
}

# Select beneficial and deleterious genes separately by their group probabilities.
directional_power_fdr_metrics <- function(
    comparison,
    cutoffs = c(0.50, 0.80, 0.95)
) {
  required_columns <- c(
    "method", "posterior_beneficial", "posterior_deleterious", "true_group"
  )
  missing_columns <- setdiff(required_columns, names(comparison))
  if (length(missing_columns) > 0) {
    stop("Comparison is missing: ", paste(missing_columns, collapse = ", "))
  }
  if (!is.numeric(cutoffs) || any(cutoffs < 0 | cutoffs > 1)) {
    stop("cutoffs must be numeric values between 0 and 1.")
  }

  result_rows <- list()
  row_number <- 1L
  for (method_name in unique(comparison$method)) {
    method_data <- comparison[comparison$method == method_name, , drop = FALSE]
    for (direction in c("beneficial", "deleterious")) {
      posterior_column <- paste0("posterior_", direction)
      true_direction <- method_data$true_group == direction
      number_true_direction <- sum(true_direction)

      for (cutoff in cutoffs) {
        # Directional selection requires high probability for this direction.
        selected <- method_data[[posterior_column]] >= cutoff
        number_selected <- sum(selected)
        number_true_selected <- sum(selected & true_direction)
        number_false_selected <- sum(selected & !true_direction)

        result_rows[[row_number]] <- data.frame(
          method = method_name,
          direction = direction,
          cutoff = cutoff,
          number_selected = number_selected,
          number_true_selected = number_true_selected,
          number_false_selected = number_false_selected,
          power = if (number_true_direction == 0) NA_real_ else
            number_true_selected / number_true_direction,
          fdr = if (number_selected == 0) 0 else
            number_false_selected / number_selected
        )
        row_number <- row_number + 1L
      }
    }
  }

  do.call(rbind, result_rows)
}

# Return non-null and directional results in one labeled data frame.
all_power_fdr_metrics <- function(
    comparison,
    null_cutoffs = c(0.05, 0.10, 0.20, 0.50),
    directional_cutoffs = c(0.50, 0.80, 0.95)
) {
  null_results <- power_fdr_metrics(comparison, cutoffs = null_cutoffs)
  null_results$type <- "non_null"
  null_results$direction <- "non_null"

  directional_results <- directional_power_fdr_metrics(
    comparison,
    cutoffs = directional_cutoffs
  )
  directional_results$type <- "directional"

  result_columns <- c(
    "method", "type", "direction", "cutoff", "number_selected",
    "number_true_selected", "number_false_selected", "power", "fdr"
  )
  rbind(
    null_results[, result_columns],
    directional_results[, result_columns]
  )
}

# Example command:
# source("truth_comparison.R")
# source("power_fdr_calc.R")
# posterior_summary <- readRDS(
#   "test_run_results/sample_test_1_posterior_summary.rds"
# )
# metadata <- readRDS("test_run_results/sample_test_1_metadata.rds")
# comparison <- compare_posterior_to_truth(
#   posterior_summary,
#   metadata = metadata
# )

# For non-directional power and FDR metrics:
# power_fdr <- power_fdr_metrics(comparison)
# View(power_fdr)

# Example command for directional power and FDR metrics:
# directional_power_fdr <- directional_power_fdr_metrics(comparison)
# View(directional_power_fdr)

# Recommended if you want to see all power and FDR metrics together in one table: 
# all_power_fdr <- all_power_fdr_metrics(comparison)
# View(all_power_fdr)
