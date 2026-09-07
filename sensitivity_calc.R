# Calculate one-vs-rest classification metrics and a three-group confusion
# matrix for every model/prior method in the comparison table.
classification_metrics <- function(comparison) {
  # These columns are produced by truth_comparison.R.
  required_columns <- c("method", "true_group", "most_probable_group")
  missing_columns <- setdiff(required_columns, names(comparison))
  if (length(missing_columns) > 0) {
    stop("Comparison is missing: ", paste(missing_columns, collapse = ", "))
  }

  # Keep a fixed group order so every confusion matrix is comparable.
  groups <- c("null", "beneficial", "deleterious")
  methods <- unique(comparison$method)
  metric_rows <- vector("list", length(methods) * length(groups))
  confusion_matrices <- setNames(vector("list", length(methods)), methods)
  row_number <- 1L

  for (method_name in methods) {
    method_data <- comparison[comparison$method == method_name, , drop = FALSE]
    # Rows are truth; columns are model predictions.
    confusion <- table(
      truth = factor(method_data$true_group, levels = groups),
      prediction = factor(method_data$most_probable_group, levels = groups)
    )
    confusion_matrices[[method_name]] <- confusion
    accuracy <- mean(method_data$true_group == method_data$most_probable_group)

    for (group_name in groups) {
      # Compute one-vs-rest counts for the current group.
      true_positive <- confusion[group_name, group_name]
      actual_positive <- sum(confusion[group_name, ])
      predicted_positive <- sum(confusion[, group_name])
      actual_negative <- sum(confusion) - actual_positive
      true_negative <- sum(confusion) - actual_positive -
        predicted_positive + true_positive

      metric_rows[[row_number]] <- data.frame(
        method = method_name,
        group = group_name,
        sensitivity = if (actual_positive == 0) NA_real_ else
          true_positive / actual_positive,
        specificity = if (actual_negative == 0) NA_real_ else
          true_negative / actual_negative,
        precision = if (predicted_positive == 0) NA_real_ else
          true_positive / predicted_positive,
        recall = if (actual_positive == 0) NA_real_ else
          true_positive / actual_positive,
        accuracy = accuracy
      )
      row_number <- row_number + 1L
    }
  }

  list(
    metrics = do.call(rbind, metric_rows),
    confusion_matrices = confusion_matrices
  )
}



##Use this command to calculate and get this code runnning  in the terminal


#posterior_summary <- readRDS(
#  "test_run_results/sample_test_2_posterior_summary.rds"
#)
#
#metadata <- readRDS(
#  "test_run_results/sample_test_2_metadata.rds"
#)

#comparison <- compare_posterior_to_truth(
#  posterior_summary,
#  metadata = metadata
#)

#classification_results <- classification_metrics(comparison)

#View(comparison)
#View(classification_results$metrics)

#To view the confusion matrices of each model, Use 
#View(classification_results$confusion_matrices)
