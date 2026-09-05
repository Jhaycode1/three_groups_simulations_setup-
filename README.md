# Three-groups Simulations
This repository houses simulation scripts for the three groups GWAS and RNA-seq model (Wixson et al. 2026). 

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
   f) `evaluation_metrics.R` - code to compute the logarithmic score and brier score from the output of `TG_nimble_model.R` and `competitors_methods.R`
   g) `posterior.R` - reads a saved posterior-sample file and creates a stacked posterior-probability plot
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


## Simualation with missingness
The real data has 34 genes that are only measured in the GWAS experiment. To run simulations under this setting we simply delete the RNA-seq information (replace with `NA`'s) for a handful of genes. There is code in `parent_script.R` to do this that is commented out.

## Calculate posterior probabilities
To calculate posterior probabilities interactively, first open `posterior.R` and set
`input_path` to the result file you want to inspect. The script currently expects
the saved file to contain `combined_piMOM_samples`, so run the `combined` model
with the `piMOM` prior before using it:

```r
input_path <- "sample_test_2.rds"
```

Then, from the R console, set the working directory to the directory containing
the result file and source the script:

```r
setwd("test_run_results")
source("../posterior.R")
```

The script displays the posterior probability table, opens `not_null_cols` with
`View()` in an interactive session, and saves the stacked bar plot as
`sample_test_2_posterior.png` in the same directory as the input file.


Citation: Troy P Wixson, Benjamin A Shaby, Daisy L Philtron, Leandro A Lima, Stacia K Wyman, Julia A Kaye, Steven Finkbeiner, International Parkinson Disease Genomics Consortium (IPDGC), A three-groups non-local model for combining heterogeneous data sources to identify genes associated with Parkinson’s disease, Biometrics, Volume 82, Issue 2, June 2026, ujag090, https://doi.org/10.1093/biomtc/ujag090
