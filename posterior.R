# Posterior summary and visualization script.
# Set this to a saved posterior-sample RDS file before sourcing the script.
input_path <- "sample_test_1.rds"

if (!file.exists(input_path)) {
  stop("Results file does not exist: ", input_path)
}

# The raw object is a named list of MCMC matrices, one per method/prior fit.
raw_samples <- readRDS(input_path)
if (!is.list(raw_samples) || length(raw_samples) == 0) {
  stop("Results file must contain a non-empty list of MCMC sample matrices.")
}

# Convert one MCMC matrix into one row per gene of posterior group summaries.
summarize_samples <- function(samples, method_name) {
  if (!is.matrix(samples) && !is.data.frame(samples)) {
    stop("Samples for ", method_name, " must be a matrix or data frame.")
  }

  samples <- as.data.frame(samples)
  not_null_cols <- grep("^not_null\\[", names(samples))
  not_ben_cols <- grep("^not_ben\\[", names(samples))
  if (length(not_null_cols) == 0 || length(not_ben_cols) == 0) {
    stop("Samples for ", method_name,
         " must contain not_null and not_ben indicator columns.")
  }
  if (length(not_null_cols) != length(not_ben_cols)) {
    stop("Indicator column counts do not match for ", method_name, ".")
  }

  # Each indicator column is one gene; rows are retained MCMC draws.
  not_null <- as.matrix(samples[, not_null_cols, drop = FALSE])
  not_ben <- as.matrix(samples[, not_ben_cols, drop = FALSE])
  posterior_null <- 1 - colMeans(not_null)
  posterior_beneficial <- colMeans((1 - not_ben) * not_null)
  posterior_deleterious <- colMeans(not_ben * not_null)

  group_probabilities <- cbind(
    posterior_null,
    posterior_beneficial,
    posterior_deleterious
  )
  # Convert indicator frequencies into mutually exclusive group probabilities.
  summary_df <- data.frame(
    method = method_name,
    gene = seq_along(posterior_null),
    posterior_null = posterior_null,
    posterior_beneficial = posterior_beneficial,
    posterior_deleterious = posterior_deleterious,
    posterior_inclusion = 1 - posterior_null,
    most_probable_group = c("null", "beneficial", "deleterious")[
      max.col(group_probabilities, ties.method = "first")
    ]
  )

  effect_columns <- list(
    posterior_mean_log_fc = grep("^log_fc\\[", names(samples)),
    posterior_mean_beta_gwas = grep("^beta_GWAS\\[", names(samples))
  )
  for (column_name in names(effect_columns)) {
    columns <- effect_columns[[column_name]]
    summary_df[[column_name]] <- if (length(columns) == nrow(summary_df)) {
      colMeans(as.matrix(samples[, columns, drop = FALSE]))
    } else {
      rep(NA_real_, nrow(summary_df))
    }
  }

  if (any(abs(rowSums(group_probabilities) - 1) > 1e-8)) {
    stop("Posterior group probabilities do not sum to one for ", method_name,
         ".")
  }

  summary_df
}

# Summarize every saved model/prior combination.
posterior_summaries <- lapply(names(raw_samples), function(method_name) {
  summarize_samples(raw_samples[[method_name]], method_name)
})
names(posterior_summaries) <- names(raw_samples)

summary_df <- do.call(rbind, posterior_summaries)
summary_path <- file.path(
  dirname(input_path),
  paste0(sub("\\.rds$", "", basename(input_path)), "_posterior_summary.rds")
)
# Save one combined table so readRDS() can be sent directly to View().
saveRDS(summary_df, summary_path)
print(summary_df)

plot_method <- if ("combined_piMOM_samples" %in% names(posterior_summaries)) {
  "combined_piMOM_samples"
} else {
  names(posterior_summaries)[1]
}
plot_df <- posterior_summaries[[plot_method]]

# Plot only when the optional visualization packages are installed.
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("tidyr", quietly = TRUE)) {
  posterior_long <- tidyr::pivot_longer(
    plot_df,
    cols = c("posterior_null", "posterior_beneficial",
             "posterior_deleterious"),
    names_to = "group",
    values_to = "posterior"
  )

  posterior_plot <- ggplot2::ggplot(
    posterior_long,
    ggplot2::aes(x = gene, y = posterior, fill = group)
  ) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::labs(
      title = plot_method,
      x = "Gene",
      y = "Posterior probability",
      fill = "Group"
    ) +
    ggplot2::theme_minimal()

  output_path <- file.path(
    dirname(input_path),
    paste0(sub("\\.rds$", "", basename(input_path)), "_posterior.png")
  )
  grDevices::png(output_path, width = 1000, height = 700, res = 120)
  print(posterior_plot)
  grDevices::dev.off()
  message("Posterior plot saved to: ", output_path)
} else {
  message("Skipping plot: ggplot2 and tidyr are required.")
}

message("Posterior summary saved to: ", summary_path)


# Calcuating the posterior probabilities for each gene being beneficial or deleterious based on the posterior summaries.
# Set the value of the imput_path in line 1 of this file to the path of the posterior summary RDS file you want to use for this calculation.

# Then source this script from the R console.
# From the project root, use:
# source("posterior.R") and then call the function like this:

# From inside the result folder, use:
# source("../posterior.R")
