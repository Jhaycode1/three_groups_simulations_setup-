#############################
# Legacy evaluation metrics from the original project.
#
# These functions are preserved for compatibility with the author's original
# evaluation workflow. The current workflow uses truth_comparison.R,
# sensitivity_calc.R, power_fdr_calc.R, and scores_calc.R instead. Do not
# source this file for the new workflow unless the legacy functions are needed.
#
# The original evaluation_function() has its own historical FDR implementation;
# the current, tested FDR implementation is in power_fdr_calc.R.
#############################
library(dplyr)

# Original author functions. See the header above before using this file.
log.score <- function(post_probs_null, true_null_status) {
  n <- length(true_null_status)
  eps <- 1e-10
  S <- 0
  for(i in 1:n){
    S <- S + sum(-log(max(eps, post_probs_null[i])) *
                   as.numeric(true_null_status[i]==0) -
                   log(max(eps,1-post_probs_null[i])) *
                   as.numeric(true_null_status[i]==1))
  }
  return(S/n)
}
brier.score <- function(post_probs_null, true_null_status) {
  n <- length(true_null_status)  # post_probs_null is the probability forecast P(G>1), so 1-lfdr or 1-P(null), C is 0 if in null and 1 if non-null (as.numeric(G>1))
  S <- sum((post_probs_null-as.numeric(true_null_status==0))^2)
  return(S/n)
}

evaluation_function <- function(post_probs_null, true_null_status, method, my_seed){
  num_true   <- sum(true_null_status)

  num_true_5               <- length(intersect(which(true_null_status != 0), which(post_probs_null < .5)))
  num_true_2               <- length(intersect(which(true_null_status != 0), which(post_probs_null < .2)))
  num_true_1               <- length(intersect(which(true_null_status != 0), which(post_probs_null < .1)))
  num_true_05              <- length(intersect(which(true_null_status != 0), which(post_probs_null < .05)))
  num_false_5              <- length(intersect(which(true_null_status == 0), which(post_probs_null < .5)))
  num_false_2              <- length(intersect(which(true_null_status == 0), which(post_probs_null < .2)))
  num_false_1              <- length(intersect(which(true_null_status == 0), which(post_probs_null < .1)))
  num_false_05             <- length(intersect(which(true_null_status == 0), which(post_probs_null < .05)))

  temp_data_frame              <- data.frame(my_seed=rep(my_seed,4))
  temp_data_frame$cut_off      <- c(.05,.10,.20,.50)
  temp_data_frame$power        <- c(num_true_05,num_true_1,num_true_2,num_true_5)/num_true
  temp_data_frame$fdr          <- c(num_false_05/length(post_probs_null<.05),
                                    num_false_1/length(post_probs_null<.1),
                                    num_false_2/length(post_probs_null<.2),
                                    num_false_5/length(post_probs_null<.5))
  temp_data_frame$logscore     <- log.score(post_probs_null,true_null_status)
  temp_data_frame$brierscore   <-  brier.score(post_probs_null,true_null_status)
  temp_data_frame$GWAS_effect  <- max(GWAS_effect)
  temp_data_frame$RNA_effect   <- max(RNA_effect)
  temp_data_frame$method       <- method
  temp_data_frame$prop_nonnull <- num_true/num_genes
  temp_data_frame$max_post_prob_null <- max(post_probs_null)
  temp_data_frame$min_post_prob_null <- min(post_probs_null)

  temp_results <- temp_data_frame %>%
    dplyr::select(method,my_seed,max_post_prob_null,min_post_prob_null,cut_off,power,fdr,logscore,brierscore,GWAS_effect,RNA_effect,prop_nonnull)
  return(temp_results)
}

