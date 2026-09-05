# Set this to the RDS file you want to inspect.
input_path <- "sample_test_1.rds"
if (!file.exists(input_path)) {
  stop("Results file does not exist: ", input_path)
}

raw_samples <- readRDS(input_path)

#Extract the table or data needed from the raw sample file
combined_pimom_samples <- raw_samples$combined_piMOM_samples

# locate columns for the group indicators
not_null_cols <- grep("not_null", colnames(combined_pimom_samples))
not_ben_cols  <- grep("not_ben", colnames(combined_pimom_samples))
if (interactive()) {
  View(not_null_cols)
}

# posterior probability gene is NULL
posterior_null <- 1 - colMeans(combined_pimom_samples[, not_null_cols])

# posterior probability gene is BENEFICIAL
posterior_beneficial <- colMeans(
    (1 - combined_pimom_samples[, not_ben_cols]) * combined_pimom_samples[, not_null_cols]
)

# posterior probability gene is DELETERIOUS
posterior_deleterious <- colMeans(
  combined_pimom_samples[, not_ben_cols] * combined_pimom_samples[, not_null_cols]
)

# assemble into a table
gene_index <- seq_along(posterior_null)

posterior_df <- data.frame(
  gene = gene_index,
  null_gene = posterior_null,
  beneficial_gene = posterior_beneficial,
  deleterious_gene = posterior_deleterious
)

posterior_df



library(ggplot2)
library(tidyr)

posterior_long <- pivot_longer(
  posterior_df,
  cols = c("null_gene", "beneficial_gene", "deleterious_gene"),
  names_to = "group",
  values_to = "posterior"
)

posterior_plot <- ggplot(posterior_long, aes(x = gene, y = posterior, fill = group)) +
  geom_col(position = "stack") +
  labs(
    x = "Gene",
    y = "Posterior probability",
    fill = "Group"
  ) +
  theme_minimal()

output_path <- file.path(
  dirname(input_path),
  paste0(sub("\\.rds$", "", basename(input_path)), "_posterior.png")
)
png(output_path, width = 1000, height = 700, res = 120)
print(posterior_plot)
dev.off()

message("Posterior plot saved to: ", output_path)