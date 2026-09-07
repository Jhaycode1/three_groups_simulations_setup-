# Calculate probability scores for posterior null probabilities.
# Lower log and Brier scores indicate better probabilistic predictions.

# Score the probability assigned to the observed null/non-null outcome.
log_score <- function(posterior_null, true_group, epsilon = 1e-10) {
  if (length(posterior_null) != length(true_group)) {
    stop("posterior_null and true_group must have the same length.")
  }
  if (any(!is.finite(posterior_null)) ||
      any(posterior_null < 0 | posterior_null > 1)) {
    stop("posterior_null must contain probabilities between 0 and 1.")
  }
  if (any(!true_group %in% c("null", "beneficial", "deleterious"))) {
    stop("true_group contains an unknown group label.")
  }

  # For a non-null truth, the relevant probability is 1 - posterior_null.
  true_probability <- ifelse(
    true_group == "null",
    posterior_null,
    1 - posterior_null
  )
  -mean(log(pmax(epsilon, true_probability)))
}

# Measure squared probability error for null versus non-null truth.
brier_score <- function(posterior_null, true_group) {
  if (length(posterior_null) != length(true_group)) {
    stop("posterior_null and true_group must have the same length.")
  }
  if (any(!is.finite(posterior_null)) ||
      any(posterior_null < 0 | posterior_null > 1)) {
    stop("posterior_null must contain probabilities between 0 and 1.")
  }
  if (any(!true_group %in% c("null", "beneficial", "deleterious"))) {
    stop("true_group contains an unknown group label.")
  }

  true_null <- as.numeric(true_group == "null")
  mean((posterior_null - true_null)^2)
}

# Apply both probability scores separately to every model/prior method.
probability_scores <- function(comparison) {
  required_columns <- c("method", "posterior_null", "true_group")
  missing_columns <- setdiff(required_columns, names(comparison))
  if (length(missing_columns) > 0) {
    stop("Comparison is missing: ", paste(missing_columns, collapse = ", "))
  }

  methods <- unique(comparison$method)
  result_rows <- lapply(methods, function(method_name) {
    method_data <- comparison[comparison$method == method_name, , drop = FALSE]
    data.frame(
      method = method_name,
      log_score = log_score(
        method_data$posterior_null,
        method_data$true_group
      ),
      brier_score = brier_score(
        method_data$posterior_null,
        method_data$true_group
      )
    )
  })

  do.call(rbind, result_rows)
}

# Example command:
# source("truth_comparison.R")
# source("scores_calc.R")
# posterior_summary <- readRDS(
#   "test_run_results/sample_test_1_posterior_summary.rds"
# )
# metadata <- readRDS("test_run_results/sample_test_1_metadata.rds")
# comparison <- compare_posterior_to_truth(
#   posterior_summary,
#   metadata = metadata
# )
# scores <- probability_scores(comparison)
# View(scores)
