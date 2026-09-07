# Three-Groups Simulation Project Guide

This guide explains the repository for a new contributor. It describes what each file does, how data move through the project, what the saved outputs mean, and which parts are intentionally deferred.

## 1. Project purpose

The project simulates a three-groups model for combining GWAS and RNA-seq evidence. Every simulated gene belongs to exactly one truth group:

- `-1`: beneficial
- `0`: null
- `1`: deleterious

The current simulator creates the truth vector with:

```r
groups <- c(
  rep(c(-1, 1), each = num_beneficial),
  rep(0, num_genes - 2 * num_beneficial)
)
```

Therefore:

```text
beneficial genes  = num_beneficial
deleterious genes = num_beneficial
null genes        = num_genes - 2 * num_beneficial
```

For example, `num_genes = 200` and `num_beneficial = 50` gives 50 beneficial, 50 deleterious, and 100 null genes. The current `num_genes = 100`, `num_beneficial = 50` setting gives zero null genes and is useful for directional testing, but it is not a sparse gene-discovery scenario.

## 2. End-to-end data flow

```text
parent_script.R
  |
  |-- sets seed, sample sizes, gene counts, effects, MCMC settings
  |
  |-- source generate_data.R
  |     |-- creates groups, GWAS genotypes/outcomes, RNA counts, covariates
  |
  |-- source TG_nimble_model.R
  |     |-- builds the NIMBLE model
  |     |-- configures reversible-jump samplers
  |     |-- runs each requested model/prior combination
  |     |-- saves posterior samples as sample_test_<n>.rds
  |
  |-- saves simulation metadata as sample_test_<n>_metadata.rds

posterior.R
  |-- reads posterior samples
  |-- calculates per-gene group probabilities
  |-- saves sample_test_<n>_posterior_summary.rds
  |-- creates a posterior probability plot

truth_comparison.R
  |-- joins model predictions to the saved simulation truth

sensitivity_calc.R
  |-- calculates classification metrics and confusion matrices

power_fdr_calc.R
  |-- calculates non-null and directional power/FDR

scores_calc.R
  |-- calculates log and Brier probability scores

calibration_calc.R
  |-- checks whether predicted non-null probabilities are calibrated
```

## 3. Main executable scripts

### `parent_script.R`

This is the main driver. It:

1. Finds the script directory and sets it as the working directory.
2. Reads optional command-line arguments for seed and output name.
3. Creates the next available result name.
4. Loads required R packages.
5. Defines sample sizes, gene counts, truth composition, effect sizes, and MCMC settings.
6. Stops early if `2 * num_beneficial > num_genes`.
7. Sources `generate_data.R`.
8. Sources `TG_nimble_model.R`.
9. Saves metadata containing the seed, settings, and exact `groups` vector.
10. Leaves competitor methods and the author's legacy evaluation code commented out.

Run it from the repository root with:

```bash
Rscript parent_script.R
```

Optional arguments are:

```bash
Rscript parent_script.R <seed> <output_name>
```

For example:

```bash
Rscript parent_script.R 2 sample_test_seed2
```

The result name is used without an extension; the model adds `.rds`.

### `generate_data.R`

This script creates the simulated inputs in the parent script's environment.

- `groups` defines the known truth for each gene.
- `beta_GWAS` converts group labels and GWAS effect sizes into signed logistic-regression effects.
- `generate_GWAS_data()` creates binary genotype predictors, an intercept, disease probabilities, and binary outcomes.
- `gene_number_index` maps flattened RNA observations back to genes.
- `PD_indicator` creates the RNA disease/control design.
- `generate_RNA_pickrell()` samples real count structure from the Montgomery/Pickrell data and adds group-specific RNA signal through `thin_diff()`.
- `library_offset` and `sex_indicator` provide RNA covariates used by the model.

The source dataset is read from `Montgomery_and_Pickrell.rds` and removed from memory after RNA generation.

### `TG_nimble_model.R`

This is the model engine.

The file has five major parts:

1. Custom probability and random-generation functions for the positive-half piMOM prior and vectorized likelihoods.
2. `my_sampler_RJ_indicator`, a reversible-jump sampler that switches a gene between null and active states while proposing its effect coefficients.
3. `three_groups_code`, the NIMBLE model for RNA-seq, GWAS, group indicators, effects, and hyperparameters.
4. Constants, monitors, initial values, and model/MCMC containers.
5. Loops over `models_for_mcmc` and `priors_for_mcmc`, compiles each model, runs `runMCMC()`, and saves the list of posterior sample matrices.

The three-group logic is hierarchical:

```text
not_null = 0                  -> null
not_null = 1 and not_ben = 0  -> beneficial
not_null = 1 and not_ben = 1  -> deleterious
```

The `not_ben` name is historical and slightly confusing: `not_ben = 0` corresponds to beneficial, while `not_ben = 1` corresponds to deleterious.

The posterior probabilities are calculated from the monitored indicators:

```r
posterior_null <- 1 - colMeans(not_null)
posterior_beneficial <- colMeans((1 - not_ben) * not_null)
posterior_deleterious <- colMeans(not_ben * not_null)
```

These are joint probabilities, so the three probabilities sum to one for each gene.

## 4. Posterior outputs

### `sample_test_<n>.rds`

A list with one MCMC sample matrix per model/prior combination, such as:

```text
combined_piMOM_samples
combined_local_samples
RNA_only_piMOM_samples
RNA_only_local_samples
GWAS_only_piMOM_samples
GWAS_only_local_samples
```

Rows are retained MCMC draws. Columns are monitored parameters, including `not_null[gene]`, `not_ben[gene]`, and model-specific effects.

### `sample_test_<n>_metadata.rds`

A list containing the exact simulation settings and truth used to produce the matching posterior file. Important fields include:

```text
seed, num_genes, num_beneficial, GWAS_effect, RNA_effect,
groups, niter, nburnin, thin, models_for_mcmc, priors_for_mcmc
```

Always pair the posterior RDS and metadata RDS with the same base name.

### `sample_test_<n>_posterior_summary.rds`

A single data frame with one row per method and gene. Important columns are:

- `posterior_null`: probability of the null group
- `posterior_beneficial`: probability of the beneficial group
- `posterior_deleterious`: probability of the deleterious group
- `posterior_inclusion`: probability of either active group
- `most_probable_group`: largest of the three probabilities
- posterior mean effect columns when monitored by that model

## 5. Current evaluation modules

The current workflow is deliberately split into small files:

### `truth_comparison.R`

Validates the posterior summary, obtains `groups` from metadata, checks every method has exactly the expected gene indices, and joins truth to predictions by gene ID rather than row position.

### `sensitivity_calc.R`

`classification_metrics()` creates one-vs-rest sensitivity, specificity, precision, recall, overall accuracy, and a three-by-three confusion matrix for every method. Rows of a confusion matrix are truth; columns are predictions.

### `power_fdr_calc.R`

- `power_fdr_metrics()` evaluates non-null discovery using `posterior_null <= cutoff`.
- `directional_power_fdr_metrics()` evaluates beneficial and deleterious discovery using their own posterior probabilities.
- `all_power_fdr_metrics()` combines both result types into one labeled table.

FDR is calculated as:

```text
number_false_selected / number_selected
```

The functions return zero FDR when no genes are selected and `NA` power when there are no true genes of the target type.

### `scores_calc.R`

Calculates one log score and one Brier score per method for null versus non-null probability. Lower scores are better.

### `calibration_calc.R`

Bins predicted non-null probabilities and compares the mean predicted probability with the observed non-null rate. Empty bins have `NA` bin statistics because no genes are available for an average.

### `evaluation_metrics.R`

Preserves the original author's legacy functions for compatibility and reference. It is intentionally not sourced by the current workflow. Current analyses use the separated modules above.

## 6. Complete current evaluation command

```r
source("truth_comparison.R")
source("sensitivity_calc.R")
source("power_fdr_calc.R")
source("scores_calc.R")
source("calibration_calc.R")

posterior_summary <- readRDS(
  "test_run_results/sample_test_1_posterior_summary.rds"
)
metadata <- readRDS("test_run_results/sample_test_1_metadata.rds")

comparison <- compare_posterior_to_truth(
  posterior_summary,
  metadata = metadata
)

classification_results <- classification_metrics(comparison)
power_fdr_results <- all_power_fdr_metrics(comparison)
probability_scores_results <- probability_scores(comparison)
calibration_results <- calibration_metrics(comparison)
```

## 7. Files intentionally deferred

- `competitors_methods.R` is not part of the current evaluation workflow.
- Competitor integration is deferred by project decision.
- MCMC convergence and mixing analysis is also deferred until the project owner chooses to begin it.
- The current short or exploratory runs must not be presented as final scientific evidence until convergence and repeated-simulation work is done.

## 8. Interpretation rules

- A model result is a prediction; `metadata$groups` is the simulation truth.
- `num_beneficial` creates the same number of deleterious genes. Null genes are whatever remains.
- A posterior probability of exactly zero or one can result from very few retained MCMC draws.
- `NA` in a calibration bin means the bin had no observations.
- `NA` in a classification metric usually means its denominator was zero.
- Always check the matching metadata file before interpreting a result.
