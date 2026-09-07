#############################
# Three-groups NIMBLE model engine.
#
# Expected objects are created by parent_script.R and generate_data.R.
# This file defines the model and runs the requested model/prior combinations.
#############################
### needs in environment:
##### models_for_mcmc, priors_for_mcmc, niter, nburnin, thin, num_samples, 
#####    dirichlet_prior, generate_GWAS_data(), generate_RNA_pickrell()
##############

library(nimble)
library(tictoc)
tic()
#############################
# Positive-half piMOM prior and helper functions.
#############################
dhalfpiMOM <- nimbleFunction(
  run = function(x=double(0), t=double(0), r=double(0), 
                 log=integer(0, default=0)){
    returnType(double(0))
    if(x == 0){
      logProb = 0
    } else {
      logProb <- (r/2)*log(2*t)-log(gamma(r/2))-(r+1)*log(x)-t/(x^2)
    }
    
    if(log) return(logProb)
    else return(exp(logProb))
  })

# helper functions to get the piMOM random number generation working
piMOM_cdf <- function(bound, t, r){
  return(integrate(function(x, t.=t, r.=r){
    return(2*t.^(r./2)/gamma(1)*x^(-(r.+1))*exp(-t./x^2))}, 
    lower=0, upper=bound)$value) }
phalfpiMOM <- nimbleRcall(function(bound=double(0), t=double(0), r=double(0)){}, 
                          Rfun = "piMOM_cdf", returnType=double(0))
piMOM_inverse <- function(a, t, r){
  uniroot(function(x) phalfpiMOM(x, t=t, r=r) - a, 
          c(0.0001, 100), extendInt = "yes")$root} 

nimble_piMOM_inverse <- nimbleRcall(function(a=double(0), 
                                             t=double(0), r=double(0)){}, 
                                    Rfun = "piMOM_inverse", 
                                    returnType=double(0))

# random draw function: needs helpers above, needed for NA's in RNA-seq data 
rhalfpiMOM <- nimbleFunction(
  run = function(n=integer(0), t = double(0), r = double(0)){
    returnType(double(0))
    if(n != 1) print("rhalfpiMOM only allows n = 1; using n = 1.")
    dev <<- runif(1, 0.01, 0.99)
    return(nimble_piMOM_inverse(dev, t=t, r=r))
  })

#############################
# Vectorized Bernoulli and negative-binomial distributions.
#############################
dBernoulliVector <- nimbleFunction(
  run = function(x    = double(1),
                 prob = double(1),
                 log  = integer(0, default = 0)) { ## default should be zero
    
    returnType(double(0))
    logProb <- sum(dbinom(x, size = 1, prob = prob, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  }
)


# for a vector of negative binomial variables
dNegBinVector <- nimbleFunction(
  run = function(x    = double(1),
                 size = double(0),
                 prob = double(1),
                 log  = integer(0, default = 0)) { ## default should be zero
    
    returnType(double(0))
    logProb <- sum(dnbinom(x, size = size, prob = prob, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  }
)

rNegBinVector <- nimbleFunction(
  run = function(n    = integer(0, default=1),
                 size = double(0),
                 prob = double(1)){
    
    returnType(double(1))
    temp_vec <- nimNumeric(length = length(prob))
    for(i in 1:length(prob)) {
      temp_vec[i] <- rnbinom(1, size=size, prob=prob[i])}
    return(temp_vec)
  })

registerDistributions(list(
  dhalfpiMOM = list(
    BUGSdist = "dhalfpiMOM(t, r)",
    range = c(0, Inf)), 
  dNegBinVector = list(
    BUGSdist = "dNegBinVector(size, prob)", 
    discrete = FALSE,
    types    = c("value=double(1)", "prob=double(1)"),
    pqAvail  = FALSE)
))

# for a vector of bernoulli variables
dBernoulliVector <- nimbleFunction(
  run = function(x    = double(1),
                 prob = double(1),
                 log  = integer(0, default = 0)) { ## default should be zero
    
    returnType(double(0))
    logProb <- sum(dbinom(x, size = 1, prob = prob, log = TRUE))
    if(log) return(logProb) else return(exp(logProb))
  }
)

#############################
# Sparse multiplication used to include only active GWAS gene effects.
#############################
sparseMult <- nimbleFunction(
  run = function(nonzeros_ind = double(1),
                 dataMat = double(2), 
                 effSizes = double(1)){
    returnType(double(1))
    nonzeros <- which(nonzeros_ind == 1)
    out <- dataMat[ , nonzeros] %*% effSizes[nonzeros]
    return(out[,1])
  })

#############################
# Reversible-jump sampler for switching a gene between null and active states.
#############################
###################
my_sampler_RJ_indicator <- nimbleFunction(
  name = 'my_sampler_RJ_indicator',
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    ## note: target is the indicator variable,
    ## control$targetNode is the variable conditionally in the model
    ## control list extraction
    coefNode      <- model$expandNodeNames(control$targetNode, 
                                           returnScalarComponents = TRUE)
    # coefNodes     <- control$targetNode
    
    proposalScale <- control$scale
    proposalMean  <- control$mean
    len_coefNode <- length(coefNode) 
                      # It is better to do this in setup code and use below
    
    
    ## node list generation
    calcNodes <- model$getDependencies(c(coefNode, target))
    calcNodesReduced <- model$getDependencies(target)
  },
  run = function() {
    currentIndicator <- model[[target]]
    if(currentIndicator == 0) {   ## propose addition of coefNode
      currentLogProb <- model$getLogProb(calcNodesReduced)
      proposalCoef <- numeric(len_coefNode)
      logProbForwardProposal <- 0
      for(l in 1:len_coefNode) {
        
        proposalCoef[l] <- rnorm(1, proposalMean, proposalScale)
        logProbForwardProposal <- logProbForwardProposal + 
          dnorm(proposalCoef[l], proposalMean, proposalScale, log = TRUE)
      }
      values(model, coefNode) <<- proposalCoef
      
      model[[target]] <<- 1
      proposalLogProb <- model$calculate(calcNodes)
      logAcceptanceProb <- proposalLogProb - currentLogProb - 
        logProbForwardProposal
    } else {                      ## propose removal of coefNode
      currentLogProb <- model$getLogProb(calcNodes)
      currentCoef <-  values(model, coefNode)
      logProbReverseProposal<- 0
      for(l in 1:len_coefNode) {
        
        logProbReverseProposal <- logProbReverseProposal + 
          dnorm(currentCoef[l], proposalMean, proposalScale, log = TRUE)
      }
      values(model, coefNode) <<- rep(0, len_coefNode) 
      # these zeros are handled in the "if(x == 0)" part of dhalfpiMOM()
      
      model[[target]] <<- 0
      model$calculate(calcNodes)
      logAcceptanceProb <- model$getLogProb(calcNodesReduced) - 
        currentLogProb + logProbReverseProposal
    }
    accept <- decide(logAcceptanceProb)
    if(accept) { copy(from = model, to = mvSaved, row = 1, 
                      nodes = calcNodes, logProb = TRUE)
    } else     { copy(from = mvSaved, to = model, row = 1, 
                      nodes = calcNodes, logProb = TRUE) }
  },
  methods = list(
    reset = function() { }
  )
)




#############################
# NIMBLE model definition.
#############################
###########
three_groups_code <- nimbleCode({
  
  # RNA-seq - iterates through genes - vectorized on individuals 
  if(which_model != "GWAS_only"){
    for(i in 1:num_genes){
      Y_RNA[1:num_individuals_RNA, i]       ~ 
        dNegBinVector(size = r_nb[gene_number_index[i]],
                      prob = p_nb[1:num_individuals_RNA, i])
      p_nb[1:num_individuals_RNA, i]        <- 
        1 / (1 + gene_effect[1:num_individuals_RNA, i] * 
               dispersion[gene_number_index[i]])
      gene_effect[1:num_individuals_RNA, i] <- 
        exp(sex_indicator[1:num_individuals_RNA] * 
              sex_effect + 
              library_offset[1:num_individuals_RNA] +
              gene_effect_baseline[gene_number_index[i]] +
              log_fc[gene_number_index[i]] *
            PD_indicator[1:num_individuals_RNA])
    }
  }
  
  # GWAS - iterates through individuals - vectorized on genes  
  if(which_model != "RNA_only"){
    # GWAS - patient level  
    Y_GWAS[1:num_individuals_GWAS] ~ 
      dBernoulliVector(p_bern[1:num_individuals_GWAS])
    p_bern[1:num_individuals_GWAS] <-
      expit(intercept_GWAS +
              sparseMult(nonzeros_ind = not_null[1:num_genes],
                         dataMat = X_GWAS[1:num_individuals_GWAS,1:num_genes],
                         effSizes = beta_GWAS[1:num_genes]))
  }
  
  # both RNA-seq and GWAS - iterates through genes
  for(j in 1:num_genes){
    # calculate the group-wise effect
    not_null[j]  ~ dbern(lambda_null)
    not_ben[j]   ~ dbern(lambda_ben)
    
  # RNA-seq
    if(which_model != "GWAS_only"){
    # gene-wise dispersion
    dispersion[j]  ~ dlnorm(meanlog = overall_dispersion, sdlog = sigma_disp)
    r_nb[j]        <- 1/dispersion[j]
    # gene effect
    gene_effect_baseline[j]  ~ dnorm(0,10^(-3))
    # group effect size
    if(which_prior == "piMOM"){
      delta_deleterious_RNA[j] ~ dhalfpiMOM(t=t_deleterious_RNA, r=2)
      delta_beneficial_RNA[j]  ~ dhalfpiMOM(t=t_beneficial_RNA, r=2)
    } else if(which_prior == "local"){
      delta_deleterious_RNA[j] ~ T(dnorm(0, sd = 0.5), 0, )
      delta_beneficial_RNA[j]  ~ T(dnorm(0, sd = 0.5), 0, )
    }
    log_fc[j]    <- not_null[j]*((1-not_ben[j])*(-1)*delta_beneficial_RNA[j] + 
                                   not_ben[j]*delta_deleterious_RNA[j]) 
    }
    
  # GWAS
    if(which_model != "RNA_only"){
    # group effect size
      if(which_prior == "piMOM"){
        delta_deleterious_GWAS[j] ~ dhalfpiMOM(t=t_deleterious_GWAS, r=2) 
        delta_beneficial_GWAS[j]  ~ dhalfpiMOM(t=t_beneficial_GWAS, r=2)
      }else if(which_prior == "local"){
        delta_deleterious_GWAS[j] ~ T(dnorm(0, sd=1), 0, ) 
        delta_beneficial_GWAS[j]  ~ T(dnorm(0, sd=1), 0, )
      }
    beta_GWAS[j] <- not_null[j]*((1-not_ben[j])*(-1)*delta_beneficial_GWAS[j] + 
                                   not_ben[j]*delta_deleterious_GWAS[j])
    }
  }
  
  # hyperparameters
  lambda_null  ~ dbeta(1, 2)
  lambda_ben   ~ dbeta(1, 1)
  
  # RNA-seq hyper-parameters
  if(which_model != "GWAS_only"){
    sex_effect    ~ dnorm(0, 10^(-3))
    overall_dispersion ~ dnorm(mean = 0, 10^(-2))
    sigma_disp         ~ T(dt(0, 1, df = 4), 0, )
    
    t_deleterious_RNA   ~ dinvgamma(shape = 2, rate = 4)
    t_beneficial_RNA    ~ dinvgamma(shape = 2, rate = 4)
  }
  
  # GWAS hyper-parameters
  if(which_model != "RNA_only"){
    intercept_GWAS ~ dnorm(0, 10^(-3))
    
    t_deleterious_GWAS   ~ dinvgamma(shape = 2, rate = 4)
    t_beneficial_GWAS    ~ dinvgamma(shape = 2, rate = 4)
  }
})  


#############################
# Build constants, monitors, data, and initial values.
#############################
consts  <- list(num_individuals_RNA=num_individuals_RNA, 
                num_individuals_GWAS=num_individuals_GWAS,
                library_offset=
                  log(library_offset[((1:num_individuals_RNA)*num_genes - 1)]), 
                num_genes=num_genes, X_GWAS=X_GWAS,
                gene_number_index=gene_number_index[1:num_genes],  
                PD_indicator=PD_indicator,
                sex_indicator=
                  sex_indicator[((1:num_individuals_RNA)*num_genes - 1)])

# set monitors 
monitors_GWAS <- c("not_null", "not_ben", "beta_GWAS", "intercept_GWAS",
                   "delta_deleterious_GWAS", "delta_beneficial_GWAS", 
                   "lambda_null", "lambda_ben", 
                   "t_deleterious_GWAS", "t_beneficial_GWAS")
monitors_RNA  <- c("not_null", "not_ben", "log_fc", "sex_effect", "dispersion", 
                   "gene_effect_baseline", 
                   "delta_deleterious_RNA", "delta_beneficial_RNA", 
                   "lambda_null", "lambda_ben", "overall_dispersion", 
                   "sigma_disp", "t_deleterious_RNA", "t_beneficial_RNA")
monitors_both <- union(monitors_GWAS, monitors_RNA)
  
# intitialize lists
models_list          <- list()
mcmc_list            <- list()
compiled_models_list <- list()
compiled_mcmc_list   <- list()
mcmc_samples_list    <- list()

# set data and inits
data    <- list(Y_RNA=matrix(Y_RNA, nrow = num_individuals_RNA, byrow = T), 
                Y_GWAS=Y_GWAS)
inits   <- list(dispersion=rgamma(num_genes,2,2), 
                gene_effect_baseline=rgamma(num_genes,2,2), 
                delta_deleterious_RNA=rep(0.5, num_genes), 
                delta_beneficial_RNA=rep(0.5, num_genes), 
                delta_beneficial_GWAS=rep(1, num_genes), 
                delta_deleterious_GWAS=rep(1, num_genes), 
                intercept_GWAS=rnorm(1, 0, 1),
                lambda_null=0.01, 
                lambda_ben=0.5,
                sex_effect=rnorm(1), 
                not_null=rbinom(num_genes, 1, 0.5), 
                not_ben=rbinom(num_genes, 1, 0.5), 
                overall_dispersion=rnorm(1, 0, 1),
                sigma_disp = rgamma(1, 2, 10), 
                t_deleterious_RNA = 3,
                t_beneficial_RNA = 3, 
                t_deleterious_GWAS = 3,
                t_beneficial_GWAS = 3)

#############################
# Run every requested model/prior combination and save posterior samples.
#############################
#############
for(model in models_for_mcmc){
  
  which_model <- model # for nimble code if statements

  if(model == "GWAS_only"){
    monitors <- monitors_GWAS
  } else if(model == "RNA_only"){
    monitors <- monitors_RNA
  } else {
    monitors <- monitors_both
  }
  
  for(prior in priors_for_mcmc){
    tic()
    which_prior         <- prior # for nimble code if statements
    model_name          <- paste(model, prior, "model", sep = "_")
    mcmc_name           <- paste(model, prior, "mcmc", sep = "_")
    conf_name           <- paste(model, prior, "conf", sep = "-")
    compiled_model_name <- paste("c", model_name, sep = "")
    compiled_mcmc_name  <- paste("c", mcmc_name, sep = "")
    
    # nimble mcmc setup
    models_list[[model_name]] <- 
      nimbleModel(code=three_groups_code, 
                  constant = consts, 
                  data = data, inits = inits)
    mcmc_list[[conf_name]]    <- 
      configureMCMC(models_list[[model_name]], 
                    monitors = monitors)
    
    ###############################
    # configure RJ samplers
    if(which_model == "combined"){
      nodeControl  = list(mean = 0, scale = 1)
      for(b in 1:length(models_list[[model_name]]$not_null)){
        nodeControl$targetNode <- 
          c(paste0("delta_beneficial_GWAS[", b, "]"),
            paste0("delta_deleterious_GWAS[", b, "]"),
            paste0("delta_beneficial_RNA[", b, "]"),
            paste0("delta_deleterious_RNA[", b, "]"))
        mcmc_list[[conf_name]]$removeSamplers(paste0("not_null[", b,"]"))
        mcmc_list[[conf_name]]$addSampler(type = my_sampler_RJ_indicator,
                                          target = paste0("not_null[", b,"]"), 
                                          control= nodeControl)
        
        # add toggled sampler for when the gene is not_null
        for(d in 1:length(nodeControl$targetNode)){
          currentConf <- 
            mcmc_list[[conf_name]]$getSamplers(nodeControl$targetNode[d])
          mcmc_list[[conf_name]]$removeSamplers(nodeControl$targetNode[d])
          mcmc_list[[conf_name]]$addSampler(type = sampler_RJ_toggled,
                                            target = nodeControl$targetNode[d],
                                            control = list(samplerType = 
                                                             currentConf[[1]]))  
        }
      }
    } else if(which_model == "RNA_only") {
      nodeControl  = list(mean = 0, scale = 1)
      for(b in 1:length(models_list[[model_name]]$not_null)){
        nodeControl$targetNode <- c(paste0("delta_beneficial_RNA[", b, "]"),
                                    paste0("delta_deleterious_RNA[", b, "]"))
        mcmc_list[[conf_name]]$removeSamplers(paste0("not_null[", b,"]"))
        mcmc_list[[conf_name]]$addSampler(type = my_sampler_RJ_indicator,
                                          target = paste0("not_null[", b,"]"), 
                                          control= nodeControl)
        
        # add toggled sampler for when the gene is not_null
        for(d in 1:length(nodeControl$targetNode)){
          currentConf <- 
            mcmc_list[[conf_name]]$getSamplers(nodeControl$targetNode[d])
          mcmc_list[[conf_name]]$removeSamplers(nodeControl$targetNode[d])
          mcmc_list[[conf_name]]$addSampler(type = sampler_RJ_toggled,
                                            target = nodeControl$targetNode[d],
                                            control = list(samplerType = 
                                                             currentConf[[1]]))  
        }
      }
    } else if(which_model == "GWAS_only"){
      nodeControl  = list(mean = 0, scale = 1)
      for(b in 1:length(models_list[[model_name]]$not_null)){
        nodeControl$targetNode <- c(paste0("delta_beneficial_GWAS[", b, "]"),
                                    paste0("delta_deleterious_GWAS[", b, "]"))
        mcmc_list[[conf_name]]$removeSamplers(paste0("not_null[", b,"]"))
        mcmc_list[[conf_name]]$addSampler(type = my_sampler_RJ_indicator,
                                          target = paste0("not_null[", b,"]"), 
                                          control= nodeControl)
        
        # add toggled sampler for when the gene is not_null
        for(d in 1:length(nodeControl$targetNode)){
          currentConf <- 
            mcmc_list[[conf_name]]$getSamplers(nodeControl$targetNode[d])
          mcmc_list[[conf_name]]$removeSamplers(nodeControl$targetNode[d])
          mcmc_list[[conf_name]]$addSampler(type = sampler_RJ_toggled,
                                            target = nodeControl$targetNode[d],
                                            control = list(samplerType = 
                                                             currentConf[[1]]))  
        }
      }
    }
    ###################################
    
    
    
    mcmc_list[[mcmc_name]]    <- buildMCMC(mcmc_list[[conf_name]], 
                                                  monitors = monitors)
    compiled_models_list[[compiled_model_name]] <- 
      compileNimble(models_list[[model_name]])
    compiled_mcmc_list[[compiled_mcmc_name]] <- 
      compileNimble(mcmc_list[[mcmc_name]], 
                    project = models_list[[model_name]])
    
    # sample from posterior
    mcmc_samples_list[[paste0(model, "_", prior, "_samples")]] <- 
      runMCMC(compiled_mcmc_list[[compiled_mcmc_name]], niter = niter, 
              nburnin = nburnin, thin = thin)
    print(paste0("sampling for ", model, " ", prior, " took: "))
    toc()
  }
}
print(paste0("Whole script took: "))
toc()

# save all samples
saveRDS(mcmc_samples_list, paste0(mcmc_save_name, ".rds"))

#############################
# Legacy in-memory posterior calculations retained for the original workflow.
# The current reusable summaries are produced by posterior.R.
#############################
#############
post_probs_null_TG <- 
  as.data.frame(lapply(mcmc_samples_list, function(x){
    incl_ind <- grep("not_null", colnames(x)) # lIne used to extract the columns of the not_null indicators from the MCMC samples
    #The line for the calculation of posterior probabbility of null group 
    return(1-colMeans(x[,incl_ind]))
    }))

post_probs_ben_TG <- 
  as.data.frame(lapply(mcmc_samples_list, function(x){
    incl_ind_null <- grep("not_null", colnames(x)) #Line used to extract the columns of the not_null indicators from the MCMC samples
    incl_ind_ben  <- grep("not_ben", colnames(x)) #Line used to extract the columns of the not_ben indicators from the MCMC samples
    # create a matrix to store the product of (1-not_ben) and not_null
    prod_mat <- matrix(NA, nrow = dim(x)[1], ncol = length(incl_ind_ben))
    for(i in 1:length(incl_ind_ben)){
      #the real code line for the calculation of the posterrior probability of the beneficial group is below
      prod_mat[,i] <- (1-x[,incl_ind_ben[i]]) * x[,incl_ind_null[i]]
    }
    # Posterior probability of the beneficial group is the mean of the product matrix across all MCMC samples, in it simple forrm
    return(colMeans(prod_mat))
    }))



