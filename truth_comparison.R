# Compare posterior classifications with the known simulation truth.
# Truth encoding: -1 = beneficial, 0 = null, 1 = deleterious.
# The function joins by gene ID, so row order cannot silently change the result.
library(dplyr)

compare_posterior_to_truth <- function(
  posterior_summary,
  groups = NULL,
  metadata = NULL
) {
  # Require the probability and prediction columns created by posterior.R.
  required_columns <- c(
    "method", "gene", "posterior_null", "posterior_beneficial",
    "posterior_deleterious", "most_probable_group"
  )
  missing_columns <- setdiff(required_columns, names(posterior_summary))
  if (length(missing_columns) > 0) {
    stop("Posterior summary is missing: ",
         paste(missing_columns, collapse = ", "))
  }
  # Prefer truth saved with the run; direct groups are supported for testing.
  if (is.null(groups) && !is.null(metadata)) {
    groups <- metadata$groups
  }
  if (is.null(groups)) {
    stop("Provide groups or simulation metadata containing groups.")
  }
  if (!is.numeric(groups) || any(!groups %in% c(-1, 0, 1))) {
    stop("groups must contain only -1, 0, and 1.")
  }

  # Convert numeric simulation labels into readable truth labels.
  truth <- data.frame(
    gene = seq_along(groups),
    true_group = c("beneficial", "null", "deleterious")[groups + 2]
  )
  # Reject missing, duplicated, or mismatched gene indices before joining.
  expected_genes <- truth$gene
  method_genes <- split(posterior_summary$gene, posterior_summary$method)
  invalid_methods <- names(method_genes)[vapply(
    method_genes,
    function(gene_ids) {
      length(gene_ids) != length(expected_genes) ||
        anyDuplicated(gene_ids) > 0 ||
        !identical(sort(as.integer(gene_ids)), expected_genes)
    },
    logical(1)
  )]
  if (length(invalid_methods) > 0) {
    stop("Gene indices do not match the simulation truth for: ",
         paste(invalid_methods, collapse = ", "))
  }

  # Join by gene and add binary and three-group correctness indicators.
  posterior_summary %>%
    dplyr::left_join(truth, by = "gene") %>%
    dplyr::mutate(
      true_nonnull = true_group != "null",
      predicted_nonnull = most_probable_group != "null",
      group_correct = most_probable_group == true_group,
      nonnull_correct = predicted_nonnull == true_nonnull
    ) %>%
    dplyr::ungroup()
}

# Example command:
# posterior_summary <- readRDS(
#   "test_run_results/sample_test_1_posterior_summary.rds"
# )
# metadata <- readRDS("test_run_results/sample_test_1_metadata.rds")
# comparison <- compare_posterior_to_truth(
#   posterior_summary,
#   metadata = metadata
# )
# View(comparison)
