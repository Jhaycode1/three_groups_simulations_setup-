#############################
# Optional competitor-method analysis.
#
# Run after generate_data.R has created Y_GWAS, X_GWAS, Y_RNA, and the design
# variables. This file is intentionally deferred from the current workflow.
# It compares conventional GWAS, edgeR, limma-voom, DESeq2, and p-value
# combination methods against the three-groups model.
#############################
library(edgeR)
library(limma)
library(qvalue)
library(MASS)
library(dplyr)
library(DESeq2)
#############################

######
# GWAS 
######
tic("GWAS")
logit_pval <- gwas_beta <- gwas_beta_sd <- rep(NA,num_genes)
for(i in 1:num_genes){
  temp            <- glm(Y_GWAS ~ X_GWAS[,i], family = binomial())
  gwas_beta[i]    <- summary(temp)$coefficients[2]
  gwas_beta_sd[i] <- summary(temp)$coefficients[4]
  logit_pval[i]   <- summary(temp)$coefficients[8]}
post_probs_null_GWAS <- qvalue(logit_pval)$lfdr
t1 <- toc()
gwas_time <- unname(t1$toc - t1$tic)


######
# edgeR
######
tic("edgeR")
groups_factor <- factor(PD_indicator)
  # store simulated count data as x
x <- matrix(Y_RNA, nrow=num_genes)
  # convert x, design into DGEList object for edgeR
y <- DGEList(counts=x,group=groups_factor) 
  # calculate normalization factors to scale raw library sizes
y <- calcNormFactors(y)             
  # create design matrix
design <- model.matrix(~ PD_indicator)   
  # estimate dispersions for each 'tag' (gene)
y <- estimateDisp(y,design)                

fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
edger_pval <- lrt$table$PValue
post_probs_null_edgeR <- qvalue(lrt$table$PValue)$lfdr
t1 <- toc()
edgeR_time <- unname(t1$toc - t1$tic)


######
# voom
######
tic("voom")
v <- voom(x, design, plot=F)
vfit <- lmFit(v, design)
efit <- eBayes(vfit)
coef_limma <- efit$coefficients[,2]
se_limma <- (sqrt(efit$s2.post) * efit$stdev.unscaled)[,2]
voom_pval <- efit$p.value[,2]
post_probs_null_voom <- qvalue(efit$p.value[,2])$lfdr  
t1 <- toc()
voom_time <- unname(t1$toc - t1$tic)

###### 
# DEseq2
###### 
tic("DEseq2")
des_x <- x    # matrix of counts
# !!!columns of the count matrix MUST match, in order,  rows of the column data
colnames_vector <- rep(NA, dim(des_x)[2])
j_control <- j_PD <- 0
for(i in 1:length(colnames_vector)){
  if(PD_indicator[i] == 0){
    j_control <- j_control + 1
    colnames_vector[i] <- paste0("control_", j_control)
  } else {
    j_PD <- j_PD + 1
    colnames_vector[i] <- paste0("PD_", j_PD)
  }
}
colnames(des_x) <- colnames_vector
des_coldata <- data.frame(condition=factor(PD_indicator))
rownames(des_coldata) <- colnames(des_x)
des_data <- DESeqDataSetFromMatrix(countData = des_x,
                                   colData = des_coldata,
                                   design = ~ condition)
des_fit <- DESeq(des_data)
#resultsNames(des_fit) # lists the coefficients
des_res <- results(des_fit, name="condition_1_vs_0")
des_res$pvalue[is.na(des_res$pvalue)] <- 1 
coef_deseq2 <- des_res$log2FoldChange
se_deseq2   <- des_res$lfcSE
post_probs_null_DEseq2 <- qvalue(des_res$pvalue)$lfdr
t1 <- toc()
DEseq2_time <- unname(t1$toc - t1$tic)

###### 
# combine results from GWAS and edgeR 
######
gwas_edger1 <- -2*(log(logit_pval)+log(edger_pval))
gwas_edger2 <- pchisq(gwas_edger1,lower.tail = F, df = 4)
gwas_edger_cauchy_stat <- 1/2 * (tan((0.5-logit_pval)*pi) + 
                                   tan((0.5-edger_pval)*pi))
gwas_edger_cauchy_pval <- pcauchy(gwas_edger_cauchy_stat, lower.tail = F)
post_probs_null_GWAS_edgeR <- qvalue(gwas_edger2)$lfdr 
post_probs_null_GWAS_edgeR_cauchy <- qvalue(gwas_edger_cauchy_pval)$lfdr

###### 
# combine results from gwas and voom
###### 
gwas_voom1 <- -2*(log(logit_pval)+log(voom_pval))
gwas_voom2 <- pchisq(gwas_voom1,lower.tail = F, df = 4)
gwas_voom_cauchy_stat <- 1/2 * (tan((0.5-logit_pval)*pi) + 
                                  tan((0.5-voom_pval)*pi))
gwas_voom_cauchy_pval <- pcauchy(gwas_voom_cauchy_stat, lower.tail = F)
post_probs_null_GWAS_voom <- qvalue(gwas_voom2)$lfdr 
post_probs_null_GWAS_voom_cauchy <- qvalue(gwas_voom_cauchy_pval)$lfdr

###### 
# combine results from gwas and deseq2
###### 
gwas_deseq1 <- -2*(log(logit_pval)+log(des_res$pvalue))
gwas_deseq2 <- pchisq(gwas_deseq1,lower.tail = F, df = 4)
gwas_deseq2_cauchy_stat <- 1/2 * (tan((0.5-logit_pval)*pi) + 
                                    tan((0.5-des_res$pvalue)*pi))
gwas_deseq2_cauchy_pval <- pcauchy(gwas_deseq2_cauchy_stat, lower.tail = F)
post_probs_null_GWAS_DEseq2 <- qvalue(gwas_deseq2)$lfdr 
post_probs_null_GWAS_DEseq2_cauchy <- qvalue(gwas_deseq2_cauchy_pval)$lfdr



compet_results <- list("logit_pval" = logit_pval, 
                       "edger_pval" = edger_pval, 
                       "voom_pval" = voom_pval, 
                       "deseq2_pval" = des_res$pval,
                       "gwas_edger_pval" = gwas_edger2, 
                       "gwas_edger_cauchy_pval" = gwas_edger_cauchy_pval, 
                       "gwas_voom_pval" = gwas_voom2,
                       "gwas_voom_cauchy_pval" = gwas_voom_cauchy_pval, 
                       "gwas_deseq_pval" = gwas_deseq2, 
                       "gwas_deseq2_cauchy_pval" = gwas_deseq2_cauchy_pval, 
                       "post_probs_null_GWAS" = post_probs_null_GWAS,
                       "post_probs_null_edgeR" = post_probs_null_edgeR, 
                       "post_probs_null_voom" = post_probs_null_voom, 
                       "post_probs_null_DEseq2" = post_probs_null_DEseq2, 
                       "post_probs_null_GWAS_edgeR" = post_probs_null_GWAS_edgeR, 
                       "post_probs_null_GWAS_edgeR_cauchy" = 
                         post_probs_null_GWAS_edgeR_cauchy, 
                       "post_probs_null_GWAS_voom" = 
                         post_probs_null_GWAS_voom, 
                       "post_probs_null_GWAS_voom_cauchy" = 
                         post_probs_null_GWAS_voom_cauchy, 
                       "post_probs_null_GWAS_DEseq2" = 
                         post_probs_null_GWAS_DEseq2, 
                       "post_probs_null_GWAS_DEseq2_cauchy" = 
                         post_probs_null_GWAS_DEseq2_cauchy, 
                       "gwas_time" = gwas_time, 
                       "edgeR_time" = edgeR_time, 
                       "voom_time" = voom_time, 
                       "DEseq2_time" = DEseq2_time)
saveRDS(compet_results, 
        paste0(save_location, "competitors_results_", save_name,"_seed", my_seed1, ".rds"))


