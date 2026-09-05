#############################
# Parent script
### calls:
###### generate_data.R, three_groups_nimble.R, pickrell_data.rds, 
######     competitors_methods.R, evaluation_metrics.R
#############################
# # save arguments from bash 
# bash_args = commandArgs(trailingOnly = T)
# # test if there is at least one argument: if not, return an error
# if (length(bash_args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# } else if (length(bash_args)==1) {
#   # default output file
#   bash_args[2] = "out.txt"
# }
# #############################
# if(!exists("bash_args")){           # for debugging
#   bash_args <- c(127, ".....TEMP....")
# }
# #############################
# my_seed1 <- as.numeric(bash_args[1])/127
script_args <- commandArgs(trailingOnly = TRUE)

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
script_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else {
  getwd()
}
setwd(script_dir)

next_save_name <- function(results_dir, prefix = "sample_test_") {
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }

  existing_files <- list.files(results_dir, full.names = FALSE)
  matching_files <- grep(
    paste0("^", prefix, "[0-9]+(?:\\..*)?$"),
    existing_files,
    value = TRUE,
    perl = TRUE
  )

  used_numbers <- if (length(matching_files) == 0) {
    integer(0)
  } else {
    as.integer(sub(paste0("^", prefix, "([0-9]+).*$"), "\\1", matching_files))
  }

  counter <- 1L
  while (counter %in% used_numbers) {
    counter <- counter + 1L
  }

  paste0(prefix, counter)
}


library(edgeR)
library(seqgendiff)
library(nimble)
library(tictoc)
library(limma)
library(qvalue)
library(MASS)
library(dplyr)
library(DESeq2)
tic("One simulation at preset values")

my_seed1 <- if (length(script_args) >= 1) as.numeric(script_args[1]) else 1
set.seed(my_seed1*127) # Used to set the number of times you want the simulation to run as the same value of seed will always produce the same data simulations. Best for comparisons and confirmation of data generation.


# !!! SAVE NAMES !!! 
save_location <- file.path(script_dir, "test_run_results")
save_name <- if (length(script_args) >= 2) {
  script_args[2]
} else {
  next_save_name(save_location)
} # This is used to create a unique output filename.



# data shape
num_individuals_RNA  <- 100 # The total number of individual or people's RNA-seq data needed that we want to simulate. The number is configurable. 
num_individuals_GWAS <- 1000 # The total number of individual or people's GWAS genetic data needed that we want to simulate. The number is configurable.Notice that the number of individuals for RNA-seq and GWAS can be different. This is because the data is collected from different sources and the number of individuals may not be the same. That is because the RNA-seq data is always lower in noise and efficient and takes place in small number of people comapred to the GWAS data that the differences between individual is low, so much number is needed for high statistical power.

num_genes            <- 10 #This set the total number of genes that we want to simulate in total. 
num_beneficial       <- 2 # set the number of genes that are beneficial or have a positive effect on the outcome. This is used to simulate the data and to evaluate the performance of the methods.
GWAS_effect          <- c(rep(0.5, 2*num_beneficial), 
                          rep(0, num_genes - 2*num_beneficial)) # set the effect size of the GWAS data. The first 2*num_beneficial genes have a positive effect of 0.5 (so as to ensure consistency and they all detected to avoid noise) and the rest have no effect. This is used to simulate the data and to evaluate the performance of the methods.
RNA_effect           <- c(rep(1.4, 2*num_beneficial), 
                          rep(0, num_genes - 2*num_beneficial)) # set the effect size of the RNA-seq data. The first 2*num_beneficial genes have a positive effect of 1.4 (so as to ensure consistency and they all detected to avoid noise) and the rest have no effect. This is used to simulate the data and to evaluate the performance of the methods.

# data generation:
source("generate_data.R") 
  # This script loads Montgomery_and_Pickrell.rds which are published data 
  #   Pickrell et. al. (2010) and Montgomery et. al. (2010) and can be accessed
  #   here:   https://bowtie-bio.sourceforge.net/recount/

# # To induce missingness uncomment the following code: 
# Y_RNA[seq(2,9902, length.out = 100)] <- NA
# Y_RNA[seq(9,9909, length.out = 100)] <- NA # nolint
# Y_RNA[seq(19,9919, length.out = 100)] <- NA
# Y_RNA[seq(29,9929, length.out = 100)] <- NA
# Y_RNA[seq(39,9939, length.out = 100)] <- NA

#############################
# Three_groups
models_for_mcmc <- c("combined", "RNA_only", "GWAS_only") #This is used to specify the models that we want to run in the MCMC simulation. The models are "combined", "RNA_only", and "GWAS_only". The combined model uses both RNA-seq and GWAS data for joint inference analysis, while the RNA_only model uses only RNA-seq data, and the GWAS_only model uses only GWAS data. This is used to evaluate the performance of the methods. # nolint
priors_for_mcmc <- c("piMOM", "local")
niter           <- 5 # This set the number of iteration that want you want the MCMC sampling to run amd simnulate the data. The number of iterations is configurable and can be set to any number. The higher the number of iterations, the more accurate the results will be, but it will also take longer to run. The lower the number of iterations, the less accurate the results will be, but it will also take less time to run. # nolint # nolint
nburnin         <- 2 # This set the number of burn-in iterations that you want the MCMC sampling to run and simulate the data. The number of burn-in iterations is configurable and can be set to any number. The higher the number of burn-in iterations, the more accurate the results will be, but it will also take longer to run. The lower the number of burn-in iterations, the less accurate the results will be, but it will also take less time to run. # nolint
thin            <- 1 # This set the thinning interval that you want the MCMC sampling to run and simulate the data. The thinning interval is configurable and can be set to any number. The higher the thinning interval, the more accurate the results will be, but it will also take longer to run. The lower the thinning interval, the less accurate the results will be, but it will also take less time to run, which is the number of steps or the result that want the model simulation to keep. Eg, here on 100 simulation, the model will keep 2 results, one at 51 and the other at 101, and discard the rest. This is used to reduce the autocorrelation in the MCMC samples and to improve the mixing of the chains. The thinning interval is set to 50, which means that every 50th sample will be kept and the rest will be discarded. This is used to reduce the autocorrelation in the MCMC samples and to improve the mixing of the chains. # nolint
num_samples     <- (niter-nburnin)/thin #
3
mcmc_save_name <- file.path(save_location, save_name)
 
source("TG_nimble_model.R")
print("Saved")
# results in post_probs_null_TG


toc()
# #############################
# # competitors
# source("competitors_methods.R")
# print("competitors_methods.R completed")

# #############################
# # Evaluate
# source("evaluation_metrics.R")
# three_groups_combined_piMOM <- 
#   evaluation_function(post_probs_null_TG$combined_piMOM_samples, 
#                       groups^2, method = "TG_combined_piMOM", my_seed=my_seed1)
# three_groups_RNA_piMOM      <- 
#   evaluation_function(post_probs_null_TG$RNA_only_piMOM_samples,
#                       groups^2, method = "TG_RNA_piMOM", my_seed=my_seed1)
# three_groups_GWAS_piMOM      <- 
#   evaluation_function(post_probs_null_TG$GWAS_only_piMOM_samples, 
#                       groups^2, method = "TG_GWAS_piMOM", my_seed=my_seed1)
# three_groups_combined_local <- 
#   evaluation_function(post_probs_null_TG$combined_local_samples,
#                       groups^2, method = "TG_combined_local", my_seed=my_seed1)
# three_groups_RNA_local      <- 
#   evaluation_function(post_probs_null_TG$RNA_only_local_samples, 
#                       groups^2, method = "TG_RNA_local", my_seed=my_seed1)
# three_groups_GWAS_local      <- 
#   evaluation_function(post_probs_null_TG$GWAS_only_local_samples, 
#                       groups^2, method = "TG_GWAS_local", my_seed=my_seed1)
# GWAS_only                   <- 
#   evaluation_function(post_probs_null_GWAS, groups^2, 
#                       method = "GWAS", my_seed=my_seed1)
# edgeR_only                  <- 
#   evaluation_function(post_probs_null_edgeR, groups^2, 
#                       method = "edgeR", my_seed=my_seed1)
# voom_only                   <- 
#   evaluation_function(post_probs_null_voom, groups^2, 
#                       method = "voom", my_seed=my_seed1)
# DEseq2_only                 <- 
#   evaluation_function(post_probs_null_DEseq2, groups^2,
#                       method = "DEseq2", my_seed=my_seed1)
# GWAS_edgeR                  <- 
#   evaluation_function(post_probs_null_GWAS_edgeR, groups^2,
#                       method = "GWAS_edgeR", my_seed=my_seed1)
# GWAS_voom                   <- 
#   evaluation_function(post_probs_null_GWAS_voom, groups^2,
#                       method = "GWAS_voom", my_seed=my_seed1)
# GWAS_DEseq2                 <- 
#   evaluation_function(post_probs_null_GWAS_DEseq2, groups^2,
#                       method = "GWAS_DEseq2", my_seed=my_seed1)

# myresults <- rbind(three_groups_combined_piMOM, three_groups_RNA_piMOM, 
#                    three_groups_GWAS_piMOM, 
#                    three_groups_combined_local, three_groups_RNA_local, 
#                    three_groups_GWAS_local,
#                    GWAS_only, edgeR_only, voom_only, DEseq2_only, 
#                    GWAS_edgeR, GWAS_voom, GWAS_DEseq2)

# saveRDS(myresults, 
#         file=paste0(save_location, "results_", save_name, ".rds"))

# myresults_null <- 
#   cbind(TG_comb_piMOM = post_probs_null_TG$combined_piMOM_samples, 
#         TG_comb_local = post_probs_null_TG$combined_local_samples, 
#         TG_RNA_poMOM = post_probs_null_TG$RNA_only_piMOM_samples,
#         TG_RNA_local = post_probs_null_TG$RNA_only_local_samples,
#         TG_GWAS_piMOM = post_probs_null_TG$GWAS_only_piMOM_samples,
#         TG_GWAS_local = post_probs_null_TG$GWAS_only_local_samples,
#         GWAS = post_probs_null_GWAS, 
#         edgeR = post_probs_null_edgeR,
#         voom = post_probs_null_voom,
#         DEseq2 = post_probs_null_DEseq2,
#         GWAS_edgeR = post_probs_null_GWAS_edgeR,
#         GWAS_voom = post_probs_null_GWAS_voom,
#         GWAS_DEseq2 = post_probs_null_GWAS_DEseq2)
# myresults_ben <- 
#   cbind(TG_comb_piMOM = post_probs_ben_TG$combined_piMOM_samples,
#         TG_comb_local = post_probs_ben_TG$combined_local_samples,
#         TG_RNA_poMOM = post_probs_ben_TG$RNA_only_piMOM_samples,
#         TG_RNA_local = post_probs_ben_TG$RNA_only_local_samples,
#         TG_GWAS_piMOM = post_probs_ben_TG$GWAS_only_piMOM_samples,
#         TG_GWAS_local = post_probs_ben_TG$GWAS_only_local_samples)

# saveRDS(myresults_null, 
#         file=paste0(save_location, "results_null_", save_name, ".rds"))
# saveRDS(myresults_ben, 
#         file=paste0(save_location, "results_ben_", save_name, ".rds"))