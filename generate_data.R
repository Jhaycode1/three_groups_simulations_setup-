#############################
# Synthetic data generation.
#
# This file is sourced by parent_script.R and intentionally writes the data
# objects used by TG_nimble_model.R into the parent environment.
#############################
library(seqgendiff)
library(edgeR)
######################
# Generate data
### needs in environment: 
##### num_individuals_RNA, num_individuals_GWAS, num_genes, num_beneficial, 
##### GWAS_effect, RNA_effect 
######################

# RNA uses pickrell dataset
# This line loads Montgomery_and_Pickrell.rds which are published data 
#   Pickrell et. al. (2010) and Montgomery et. al. (2010) and can be accessed
#   here:   https://bowtie-bio.sourceforge.net/recount/
pickrell_data     <- readRDS("Montgomery_and_Pickrell.rds") 

#############################
# Define truth labels and generate GWAS data first. The GWAS outcome is also
# used to establish the disease setting for the RNA-seq simulation.
#############################
# three groups: -1 is beneficial, 0 is null, 1 is deleterious
groups           <- c(rep(c(-1, 1), each = num_beneficial), ## E.g -1,-1, 1, 1, 0, 0, 0, 0, 0, 0)
                      rep(0, num_genes - 2*num_beneficial))
effect_size_GWAS <- GWAS_effect
beta_GWAS        <- groups*effect_size_GWAS

# Generate binary genotype predictors and a binary disease outcome.
generate_GWAS_data <- function(num_ind = num_individuals_GWAS, 
                               num_gen = num_genes, 
                               beta = beta_GWAS){
  X_GWAS <<- matrix(NA, num_ind, num_gen)
  
  # allow each gene to have a different proportion of minor alleles
  for(i in 1:num_gen){
    temp_proportion_alleles <- rbeta(1, 20, 35)
    X_GWAS[,i] <<- rbinom(num_ind, 1, temp_proportion_alleles)
  }
  
  intercept_GWAS <<- rnorm(1, 0, sd = 1)
  
  p_GWAS <<- exp(intercept_GWAS + X_GWAS%*%beta) / 
    (1+exp(intercept_GWAS + X_GWAS%*%beta))
  
  Y_GWAS       <<- rbinom(num_ind, 1, p_GWAS)
}
generate_GWAS_data()



#############################
# RNA-seq setup and simulation.
#############################
gene_number_index <- rep(1:num_genes, num_individuals_RNA)
PD_indicator <- rbinom(num_individuals_RNA, 1, 0.5)

generate_RNA_pickrell <- function(data = pickrell_data, 
                                  num_ind = num_individuals_RNA, 
                                  num_gen = num_genes, 
                                  PD_ind = PD_indicator, 
                                  fc_effect_size = RNA_effect, 
                                  gene_groups = groups){
  pickrell_data_subsampled <- select_counts(mat = data, 
                                            nsamp = num_ind, 
                                            ngene = num_gen, 
                                            filter_first = T)
  
  # Build the disease design and gene-specific log fold-change matrix.
  design_mat <- matrix(PD_ind, ncol = 1)
  coef_mat   <- matrix(log2(fc_effect_size^gene_groups), ncol = 1)
  
  # Thin real count data to add the requested group-specific signal.
  thinned_output <- thin_diff(mat = pickrell_data_subsampled, 
                              design_fixed = design_mat, 
                              coef_fixed = coef_mat)
  
  Y_RNA          <<- as.vector(thinned_output$mat)
  
  library_offset <- colSums(thinned_output$mat)
  library_offset <<- rep(library_offset, each = num_gen)
  
  # Generate individual-level covariates used by the RNA likelihood.
  sex_indicator   <<- rbinom(num_ind, 1, 0.5)  # covariate (sex)
  sex_indicator   <<- rep(sex_indicator, each=num_genes)
}
generate_RNA_pickrell()

rm(pickrell_data)
