# Three-groups Simulations
This repository houses simulation scripts for the three groups GWAS and RNA-seq model (Wixson et al. 2026). 

For a file-by-file explanation of the model, data flow, result files, and
current evaluation workflow, see [PROJECT_GUIDE.md](PROJECT_GUIDE.md).

The code for the real data analysis is the same as in the TG_nimble_model.R except for two noticable changes. First, the real data includes addition of other covariates. Second, the script used for the real data anlysis does not loop over different models, each model was run on its own script. Due to data privacy concerns we cannot share the data and thus sharing the script for that analysis is not useful. 


### **WARNING**: these simulations take a long time to run. It takes around 65 minutes to run on an M2pro macbook laptop with the `parent_script.R` options as they are. 


## To Run
The simulation can be run entirely from R. The `parent_script.R` file controls the simulation, calls the data-generation and model scripts, and creates the output filename automatically.

## To run a single simulation:
1. The following files ought to be in the same directory:
   a) `parent_script.R` - calls the relevant scripts and controls one simulation run
   b) `generate_data.R` - generates synthetic data for simulation studies
   c) `Montgomery_and_Pickrell.rds` - RNA-seq data from Pickrell et. al. (2010) and Montgomery et. al. (2010) [accessible here](https://bowtie-bio.sourceforge.net/recount/)
   d) `TG_nimble_model.R` - three-groups model in nimble, this will loop over models and priors and save results after running
   e) `competitors_methods.R` - code to run conventional GWAS, edgeR, limma-voom, and DESeq2 with both the Fisher and Cauchy p-value combinations 
   f) `truth_comparison.R` - compares posterior classifications with the known simulated truth
   g) `sensitivity_calc.R` - calculates sensitivity, specificity, precision, recall, accuracy, and confusion matrices
   h) `power_fdr_calc.R` - calculates overall and directional power and false discovery rate
   i) `scores_calc.R` - calculates logarithmic and Brier scores
   j) `calibration_calc.R` - checks calibration of posterior non-null probabilities
   k) `evaluation_metrics.R` - preserves the author's original evaluation helpers; not used by the current workflow
   l) `posterior.R` - reads a saved posterior-sample file and creates posterior summaries and a stacked plot

The matching metadata file, `sample_test_2_metadata.rds`, contains the exact
simulation truth and settings used to create the posterior samples. Always use
the posterior summary and metadata files with the same base name.
2. The simulation requires these packages: `edgeR`, `seqgendiff`, `nimble`, `tictoc`, `limma`, `qvalue`, `MASS`, `dplyr`, and `DESeq2`. The posterior plot script additionally requires `ggplot2` and `tidyr`.
3. Open a terminal in the repository directory and run:

   ```bash
   Rscript parent_script.R
   ```

4. To choose which models and priors are run, edit these vectors in `parent_script.R`:

   ```r
   models_for_mcmc <- c("combined")
   priors_for_mcmc <- c("piMOM", "local")
   ```

   This example runs only the `combined` model with both priors. To run only the combined model with the `piMOM` prior, use:

   ```r
   models_for_mcmc <- c("combined")
   priors_for_mcmc <- c("piMOM")
   ```

5. Adjust `niter`, `nburnin`, and `thin` in `parent_script.R` as needed. Running the script will:
   - Generate synthetic data using the `generate_data.R` script. GWAS data are generated from a logistic regression model. RNA-seq data are generated through randomly selecting a subset of a real RNA-seq dataset and adding signal to the specified genes (this is done using the binomial thinning of Gerard (2020)).
   - Run `TG_nimble_model.R` which performs MCMC sampling of the three groups method on the same generated data for each combination of model and gene-effect prior you have specified.
   - Run the competitor methods and compute evaluation metrics when those sections are enabled in `parent_script.R`.
   - Save posterior samples in `test_run_results/sample_test_<n>.rds`.

The script chooses the lowest available counter. For example, if `sample_test_1.rds` and `sample_test_3.rds` exist, the next output is `sample_test_2.rds`. If `sample_test_1.rds` is later deleted, the next run reuses `sample_test_1.rds`. No Python file is needed. The competitor-method and evaluation sections are currently commented out in `parent_script.R`; uncomment them if those additional analyses are needed.

The current example settings use `num_genes = 100` and `num_beneficial = 50`.
Because the simulator creates the same number of deleterious genes, this
scenario has 50 beneficial, 50 deleterious, and zero null genes. For FDR and
null-calibration studies, use more genes than `2 * num_beneficial`, for example
`num_genes = 200` and `num_beneficial = 50`.


## Simualation with missingness
The real data has 34 genes that are only measured in the GWAS experiment. To run simulations under this setting we simply delete the RNA-seq information (replace with `NA`'s) for a handful of genes. There is code in `parent_script.R` to do this that is commented out.

## Calculate posterior probabilities
To calculate posterior probabilities interactively, first open `posterior.R` and set
`input_path` to the result file you want to inspect:

```r
input_path <- "sample_test_2.rds"
```

The script validates and summarizes every model/prior sample matrix in the result
file. It saves per-method summaries as an RDS file and creates a stacked plot for
the combined piMOM fit when that fit is available. From the R console, set the
working directory to the directory containing the result file and source the script:

```r
setwd("test_run_results")
source("../posterior.R")
```

The script prints the posterior probability table, saves
`sample_test_2_posterior_summary.rds`, and saves the stacked bar plot as
`sample_test_2_posterior.png` in the same directory as the input file.

The evaluation files are run in this order:

```r
source("truth_comparison.R")
source("sensitivity_calc.R")
source("power_fdr_calc.R")
source("scores_calc.R")
source("calibration_calc.R")
```

`evaluation_metrics.R` is intentionally not sourced here. It preserves the
author's original functions for reference and compatibility; the separate files
above are the current workflow and do not depend on the legacy implementation.

Then create the truth comparison and classification metrics:

```r
posterior_summary <- readRDS(
   "test_run_results/sample_test_2_posterior_summary.rds"
)
metadata <- readRDS("test_run_results/sample_test_2_metadata.rds")
comparison <- compare_posterior_to_truth(
   posterior_summary,
   metadata = metadata
)
classification_results <- classification_metrics(comparison)
View(comparison)
View(classification_results$metrics)

power_fdr_results <- all_power_fdr_metrics(comparison)
View(power_fdr_results)

probability_scores_results <- probability_scores(comparison)
View(probability_scores_results)

calibration_results <- calibration_metrics(comparison)
View(calibration_results)
```

`power_fdr_results` is one table. The `type` column identifies whether a row
evaluates all non-null genes or one direction, and `direction` identifies
`non_null`, `beneficial`, or `deleterious`.

Citation: Troy P Wixson, Benjamin A Shaby, Daisy L Philtron, Leandro A Lima, Stacia K Wyman, Julia A Kaye, Steven Finkbeiner, International Parkinson Disease Genomics Consortium (IPDGC), A three-groups non-local model for combining heterogeneous data sources to identify genes associated with Parkinson’s disease, Biometrics, Volume 82, Issue 2, June 2026, ujag090, https://doi.org/10.1093/biomtc/ujag090
