# Calculating the 5-fold cross validation log posterior likelihood
# for each of the models we want to compare. 
# _V2 to reflect the post-R update and adding additional models. 

library(foreach)
library(doParallel)
library(PanCanVarSel)
library(spBayesSurv)
library(survival)
library(Matrix)

# Saving important working directories
currentwd = "~/PanTCGA/SpikeAndSlabResults/SpikeAndSlabGibbsSampler/ModelWithInterceptOutsideSS/ResultsAfterFiltering_V3/ModelComparison/"
PanTCGAwd = "~/PanTCGA/"
modelwd = "~/PanTCGA/SpikeAndSlabResults/SpikeAndSlabGibbsSampler/ModelWithInterceptOutsideSS/"
datawd = "~/PanTCGA/SpikeAndSlabResults/SpikeAndSlabGibbsSampler/ModelWithInterceptOutsideSS/ResultsAfterFiltering_V3/"

# load in the Covariates, Survival, Censored data 
load(paste(datawd, "XYC_V2_WithAge_StandardizedPredictors.rda", sep = ""))

# load in helper functions
source(paste(PanTCGAwd, "HelperFunctions_V2.R", sep = ""))

# load in model with spike and slab and intercept outside s/s framework
source(paste(modelwd, "ExtendedHierarchicalModelLogNormalSpikeAndSlabInterceptOutsideSS_V3.R", sep = ""))

# Important variables
covariates_by_cancer_full = lapply(Covariates, function(cancer) as.numeric(colnames(cancer))) # the covariates each cancer type has (each should have a 1 for the intercept)
covariates_table = table(unlist(covariates_by_cancer_full)) # gives the number of cancer types who have data on each PC
covariates_in_model_full = as.numeric(names(covariates_table))

covariates_in_model = covariates_in_model_full
p = length(covariates_in_model)
covariates_by_cancer = covariates_by_cancer_full # the name of this variable in the HelperFunctions.R file
iters = 100000
n_cancer = length(Covariates)

cancer_types = names(Covariates)

# Generating starting values
beta_tilde_start = rnorm(p, 0, 1)
sigma2_start = 1/rgamma(1, 1, 1)
lambda2_start = 1/rgamma(p, 1, 1)
pi_start = rep(0.5, p)

# Initialize inclusion indicators for every cancer type as in the slab
gamma_start = lapply(1:n_cancer, function(type) {
  indicators = rbinom(length(covariates_by_cancer[[type]]), size = 1, prob = 1)
  indicators[1] = 1
  names(indicators) = covariates_by_cancer[[type]]
  list(indicators)
})

# For hierarchical spike/slab model
starting_values = list(sigma2_start = sigma2_start,
                       beta_tilde_start = beta_tilde_start,
                       lambda2_start = lambda2_start,
                       gamma_start = gamma_start, # this will get passed on to be the list that stores posterior draws
                       pi_start = pi_start)

# For fixed at 0.5 model
pi_start_fixed_at_0.5 = rep(0.5, 1)

starting_values_fixed_at_0.5 = list(sigma2_start = sigma2_start,
                                    beta_tilde_start = beta_tilde_start,
                                    lambda2_start = lambda2_start,
                                    gamma_start = gamma_start, # this will get passed on to be the list that stores posterior draws
                                    pi_start = pi_start_fixed_at_0.5)

# For fixed at 1.0 model
pi_start_fixed_at_1.0 = rep(1, 1)

starting_values_fixed_at_1.0 = list(sigma2_start = sigma2_start,
                                    beta_tilde_start = beta_tilde_start,
                                    lambda2_start = lambda2_start,
                                    gamma_start = gamma_start, # this will get passed on to be the list that stores posterior draws
                                    pi_start = pi_start_fixed_at_1.0)

# For shared across betas model
pi_start_shared_across_betas = rep(0.5, 1)

starting_values_shared_across_betas = list(sigma2_start = sigma2_start,
                                           beta_tilde_start = beta_tilde_start,
                                           lambda2_start = lambda2_start,
                                           gamma_start = gamma_start, # this will get passed on to be the list that stores posterior draws
                                           pi_start = pi_start_shared_across_betas)

# For null (intercept only) model
Covariates_null = lapply(Covariates, function(type) type[, 1, drop = FALSE])

covariates_in_model_null = covariates_in_model_full[1]
p_null = length(covariates_in_model_null)
covariates_by_cancer_null = lapply(covariates_by_cancer, function(type) type[1]) # just storing the intercept for every cancer

# Generating starting values
beta_tilde_start_null = rnorm(p_null, 0, 1)
sigma2_start_null = 1/rgamma(1, 1, 1)
lambda2_start_null = 1/rgamma(p_null, 1, 1)
pi_start_null = rep(0.5, p_null)

# Initialize inclusion indicators for every cancer type as in the slab
gamma_start_null = lapply(1:n_cancer, function(type) {
  indicators = rbinom(length(covariates_by_cancer_null[[type]]), size = 1, prob = 1)
  indicators[1] = 1
  names(indicators) = covariates_by_cancer_null[[type]]
  list(indicators)
})

starting_values_null = list(sigma2_start = sigma2_start_null,
                            beta_tilde_start = beta_tilde_start_null,
                            lambda2_start = lambda2_start_null,
                            gamma_start = gamma_start_null, # this will get passed on to be the list that stores posterior draws
                            pi_start = pi_start_null)


# Creating the priors argument
priors = list(betatilde_priorvar_intercept = 10^2,    # \tilde\beta_0 ~ N(0,100). both beta tildes are centered at 0
              betatilde_priorvar_coefficient = 1,     # -> \tilde\beta_p ~ N(0,1), p > 0
              lambda2_priorshape_intercept = 1,       # -> \lambda^2_0 ~ IG(1,1)
              lambda2_priorrate_intercept = 1,         
              lambda2_priorshape_coefficient = 5,     # -> \lambda^2_p ~ IG(5,1) p > 0
              lambda2_priorrate_coefficient = 1,
              sigma2_priorshape = .01,                # -> \sigma^2 ~ IG(0.01, 0.01)
              sigma2_priorrate = .01,
              spike_priorvar = 1/10000)

# Generating indices for the training data
TrainingSet = Generate5FoldTrainingSet(Covariates, Survival, cancer_types, n_vec)

###################################################################################
### Combining all models runs into a single parallel computation ##################
###################################################################################

# Attempting to print progress results during run
writeLines(c(""), paste(currentwd, "PosteriorLikelihoodForAllModelsLog", ".txt", sep = ""))

start <- Sys.time()
cl <- makeCluster(5, outfile=paste0(currentwd, "LogAllModels.txt"))
registerDoParallel(cl)
PL = foreach(cv_iter = 1:5, .packages = c("MASS", "truncnorm", "EnvStats", "Matrix", "PanCanVarSel")) %dopar% {
  current_train_ind = sapply(TrainingSet[[1]], '[', cv_iter)
  Covariates_Training = lapply(1:length(Covariates), function(k) Covariates[[k]][current_train_ind[[k]], ])
  Survival_Training = lapply(1:length(Survival), function(k) Survival[[k]][current_train_ind[[k]]])
  Censored_Training = lapply(1:length(Censored), function(k) Censored[[k]][current_train_ind[[k]]])
  
  sink(paste(currentwd, "PosteriorLikelihoodForAllModelsLog", ".txt", sep = ""), append=TRUE)
  cat(paste("iteration",cv_iter, "has started at time", Sys.time(), "\n"))
  
  # ---------------------------------------------------------------------------
  # Hierarchical Spike/Slab model
  posteriors_hierarchical_ss = HierarchicalLogNormalSpikeSlab(Covariates_Training, Survival_Training, Censored_Training, 
                                        starting_values, iters, covariates_in_model, 
                                        covariates_by_cancer, priors, pi_generation = "shared_across_cancers",
                                        progress = "model_comparison")
  
  cat(paste("hierarchical s/s model in parallel process",cv_iter, "has finished at time", Sys.time(), "\n"))
  
  # ---------------------------------------------------------------------------
  # Model with prior inclusion probability fixed at 0.5
  posteriors_fixed_0.5 = HierarchicalLogNormalSpikeSlab(Covariates_Training, Survival_Training, Censored_Training, 
                                        starting_values_fixed_at_0.5, iters, covariates_in_model, 
                                        covariates_by_cancer, priors, pi_generation = "fixed_at_0.5",
                                        progress = "model_comparison")
  
  cat(paste("fixed at 0.5 model in parallel process",cv_iter, "has finished at time", Sys.time(), "\n"))
  
  # ---------------------------------------------------------------------------
  # Model with prior inclusion probability fixed at 1.0 (the full model)
  posteriors_fixed_1.0 = HierarchicalLogNormalSpikeSlab(Covariates_Training, Survival_Training, Censored_Training, 
                                        starting_values_fixed_at_1.0, iters, covariates_in_model, 
                                        covariates_by_cancer, priors, pi_generation = "fixed_at_1.0",
                                        progress = "model_comparison")
  
  cat(paste("fixed at 1.0 model in parallel process",cv_iter, "has finished at time", Sys.time(), "\n"))
  
  # ---------------------------------------------------------------------------
  # Model where prior inclusion probability is the same for all betas
  posteriors_shared_across_betas = HierarchicalLogNormalSpikeSlab(Covariates_Training, Survival_Training, Censored_Training, 
                                        starting_values_shared_across_betas, iters, covariates_in_model, 
                                        covariates_by_cancer, priors, pi_generation = "shared_across_betas",
                                        progress = "model_comparison")
  
  cat(paste("shared across betas model in parallel process",cv_iter, "has finished at time", Sys.time(), "\n"))
  
  # ---------------------------------------------------------------------------
  # Null model
  Covariates_Training_Null = lapply(1:length(Covariates_null), function(k) Covariates_null[[k]][current_train_ind[[k]], , drop = FALSE])
  
  posteriors_null = HierarchicalLogNormalSpikeSlab(Covariates_Training_Null, Survival_Training, Censored_Training, 
                                              starting_values_null, iters, covariates_in_model_null, 
                                              covariates_by_cancer_null, priors, pi_generation = "fixed_at_1.0",
                                              progress = "model_comparison")
  
  cat(paste("null (intercept only) model in parallel process",cv_iter, "has finished at time", Sys.time(), "\n"))
  
  # ---------------------------------------------------------------------------
  # The joint approach (concatenating all data matrices together, applying model as though 1 cancer type). 
  CovariatesCombined <- list(do.call(rbind, CombineCancerDataTogether(Covariates, covariates_in_model_full)))
  SurvivalCombined <- list(unlist(Survival))
  CensoredCombined <- list(unlist(Censored))
  
  # Combining the training data
  CovariatesTrainingCombined <- list(do.call(rbind, CombineCancerDataTogether(Covariates_Training, covariates_in_model_full)))
  SurvivalTrainingCombined <- list(unlist(Survival_Training))
  CensoredTrainingCombined <- list(unlist(Censored_Training))
  
  # Adjusting some of the standard inputs for the hierarchical model:
  covariates_by_cancer_joint <- list(covariates_in_model_full) # since all the cancers are combined and have 0s filled in, I am just creating a new covariates_by_cancer specific to this scenario.
  gamma_start_joint <- list(list(rep(1, length(covariates_in_model_full)))) # since there is only "1" cancer type here, we only need one list for the inclusion indicators. We thus have a list of one list. 
  starting_values_joint <- starting_values; starting_values_joint$gamma_start <- gamma_start_joint
  
  posteriors_joint = HierarchicalLogNormalSpikeSlab(CovariatesTrainingCombined, SurvivalTrainingCombined, CensoredTrainingCombined,
                                                    starting_values_joint, iters, covariates_in_model_full, covariates_by_cancer_joint,
                                                    priors, pi_generation = "fixed_at_0.5",
                                                    progress = "model_comparison")
  
  cat(paste("the joint model in parallel process",cv_iter, "has finished at time", Sys.time(), "\n"))
  
  # ---------------------------------------------------------------------------
  # The separate approach (run the model separately for each cancer type and combine posterior likelihoods at the end)
  PosteriorsByCancerSeparately <- lapply(1:n_cancer, function(i) list()) # to store the posteriors for each cancer type
  names(PosteriorsByCancerSeparately) <- cancer_types
  
  for (type in  1:n_cancer) { # for each cancer type
    # Subsetting just the data for a single cancer type
    CovariatesTraining_type <- Covariates_Training[type]
    SurvivalTraining_type <- Survival_Training[type]
    CensoredTraining_type <- Censored_Training[type]
    
    # Adjusting these parameters to reflect just a single cancer type
    covariates_by_cancer_type <- covariates_by_cancer_full[type]
    covariates_in_model_type <- unlist(covariates_by_cancer_type)
    
    # Selecting just the starting values that correspond to this cancer type
    avail <- covariates_in_model_full %in% covariates_by_cancer_type[[1]]
    starting_values_type <- list(sigma2_start = starting_values$sigma2_start,
                                 beta_tilde_start = starting_values$beta_tilde_start[avail],
                                 lambda2_start = starting_values$lambda2_start[avail],
                                 gamma_start = list(list(starting_values$gamma_start[[type]][[1]])),
                                 pi_start = starting_values$pi_start[avail])
    
    # Running the model
    posteriors_type <- HierarchicalLogNormalSpikeSlab(CovariatesTraining_type, SurvivalTraining_type, CensoredTraining_type,
                                                      starting_values_type, iters, covariates_in_model_type, covariates_by_cancer_type,
                                                      priors, pi_generation = "fixed_at_0.5",
                                                      progress = "model_comparison")
    # Storing the results
    PosteriorsByCancerSeparately[[type]] <- posteriors_type
  }
  
  cat(paste("the separated model in parallel process",cv_iter, "has finished at time", Sys.time(), "\n"))
  
  # ---------------------------------------------------------------------------
  # Computing the posterior likelihoods
  Results = c(PosteriorLikelihood(posteriors_hierarchical_ss, Covariates, Survival, Censored, current_train_ind, iters/2),
  PosteriorLikelihood(posteriors_fixed_0.5, Covariates, Survival, Censored, current_train_ind, iters/2),
  PosteriorLikelihood(posteriors_fixed_1.0, Covariates, Survival, Censored, current_train_ind, burnin = iters/2),
  PosteriorLikelihood(posteriors_shared_across_betas, Covariates, Survival, Censored, current_train_ind, iters/2),
  PosteriorLikelihood(posteriors_null, Covariates_null, Survival, Censored, current_train_ind, burnin = iters/2),
  PosteriorLikelihood(posteriors_joint, Covariates, Survival, Censored, current_train_ind, burnin = iters/2, joint = TRUE))
  
  # Computing the posterior likelihood for the separated model
  joint_ll_for_each_cancer <- lapply(1:n_cancer, function(type) {
    JointLogLikelihoodForSeparateModel(posteriors = PosteriorsByCancerSeparately[[type]],
                                       X = Covariates[type],
                                       Y = Survival[type],
                                       Censored = Censored[type], 
                                       current_train_ind[type], 
                                       burnin = iters/2)
  })
  
  # Return the posterior likelihood
  iters_used_in_ll <- (iters/2/10) + 1
  PL_separated <- -log(iters_used_in_ll) + logSum(Reduce('+', joint_ll_for_each_cancer))
  
  # Combining the results
  Results <- c(Results, PL_separated)
  
  names(Results) = c("Hierarchical S/S", "Fixed at 0.5", "Fixed at 1.0", "Shared Across Betas", "Null Model", 
                     "Joint Model", "Separated Model")
  
  Results
}
stopCluster(cl)
end <- Sys.time()
end-start

# Computing the results
PL.Matrix = do.call(rbind, PL)
PL.Means = colMeans(PL.Matrix)


# -----------------------------------------------------------------------------
# Computing the out-of-sample posterior likelihood for the Maity et al. model
# separately. 
# -----------------------------------------------------------------------------

# Printing the progress while it runs
writeLines(c(""), paste(currentwd, "PosteriorLikelihoodForMaityModelLog", ".txt", sep = ""))

iters = 50000
start <- Sys.time()
cl <- makeCluster(5, outfile = paste0(currentwd, "Log.txt"))
registerDoParallel(cl)
PL.Maity = foreach(cv_iter = 1:5, .packages = c("MASS", "truncnorm", "EnvStats", "Matrix", "PanCanVarSel")) %dopar% {
  current_train_ind = sapply(TrainingSet[[1]], '[', cv_iter)
  Covariates_Training = lapply(1:length(Covariates), function(k) Covariates[[k]][current_train_ind[[k]], ])
  Survival_Training = lapply(1:length(Survival), function(k) Survival[[k]][current_train_ind[[k]]])
  Censored_Training = lapply(1:length(Censored), function(k) Censored[[k]][current_train_ind[[k]]])
  
  sink(paste(currentwd, "PosteriorLikelihoodForMaityModelLog", ".txt", sep = ""), append=TRUE)
  cat(paste("iteration",cv_iter, "has started at time", Sys.time(), "\n"))
  
  # ---------------------------------------------------------------------------
  # Implementing the Maity et al. (2019) approach. 
  # Combining the training data
  CovariatesTrainingCombinedPCVS <- PrepareCancerDataForPanCanVarSel(Covariates_Training, covariates_in_model)
  SurvivalTrainingCombined <- unlist(Survival_Training) 
  CensoredTrainingCombined <- unlist(Censored_Training)
  n_vec_training <- sapply(Covariates_Training, nrow)
  
  # Putting the survival and censor times together with a right-censor indicator
  SurvivalCensoredTraining <- SurvivalTrainingCombined
  SurvivalCensoredTraining[is.na(SurvivalCensoredTraining)] <- CensoredTrainingCombined[is.na(SurvivalTrainingCombined)] # filling in the censored times with the right-censor times
  ct <- as.matrix(cbind(time = SurvivalCensoredTraining, 
                        censored = !is.na(unlist(SurvivalTrainingCombined)))) # 1 = event time, 0 = censored (according to package documentation but double-check)
  
  start <- Sys.time()
  posteriors_horseshoe <- hsaftgroupcorr(ct = ct, X = CovariatesTrainingCombinedPCVS, method.tau = "truncatedCauchy",
                                         method.sigma = "Jeffreys", r = n_cancer, alpha = 0.05, burn = iters,
                                         nmc = iters, n.seq = n_vec_training, pk = p)
  end <- Sys.time()
  end-start
  
  cat(paste("Maity et al. Model in parallel process", cv_iter, "has finished at time", Sys.time(), "\n"))
  
  # Computing the posterior likelihood for the PanCanVarSel model
  # Putting the posterior samples in the format expected by the PosteriorLikelihood function
  ind_for_betas_by_cancer <- seq(1, p*n_cancer + p, by = p)
  betas_horseshoe <- lapply(1:n_cancer, function(type) { # for each cancer type
    betas_type <- posteriors_horseshoe$BetaSamples[ind_for_betas_by_cancer[type]:(ind_for_betas_by_cancer[(type+1)]-1),]
    betas_type_list <- as.list(data.frame(betas_type)) # converting each column to a list element
    
    # Checking that converting the matrix to a list worked
    # match <- c(); for (i in 1:50) match[i] <- all.equal(betas_type_list[[i]], betas_type[,i]); all(match)
    
    names(betas_type_list) <- NULL
    betas_type_list
  })
  PL_horseshoe <- PosteriorLikelihood(posteriors = list(betas = betas_horseshoe, sigma2 = posteriors_horseshoe$Sigma2Samples),
                                      X = CombineCancerDataTogether(Covariates, covariates_in_model), Y = Survival, Censored = Censored, current_train_ind, burnin = 1)
  
  # Returning the posterior likelihood for the model
  PL_horseshoe
}
stopCluster(cl)
end <- Sys.time()
end-start

# Taking the average
PL.Maity.Mean <- mean(unlist(PL.Maity))


# -----------------------------------------------------------------------------
# Computing the out-of-sample posterior likelihood for the hierarchical 
# spike-and-slab model by itself
# -----------------------------------------------------------------------------

writeLines(c(""), paste(currentwd, "PosteriorLikelihoodForHierarchicalSS", ".txt", sep = ""))

iters = 100000
start <- Sys.time()
cl <- makeCluster(5)
registerDoParallel(cl)
PL = foreach(cv_iter = 1:5, .packages = c("MASS", "truncnorm", "EnvStats", "Matrix")) %dopar% {
  current_train_ind = sapply(TrainingSet[[1]], '[', cv_iter)
  Covariates_Training = lapply(1:length(Covariates), function(k) Covariates[[k]][current_train_ind[[k]], ])
  Survival_Training = lapply(1:length(Survival), function(k) Survival[[k]][current_train_ind[[k]]])
  Censored_Training = lapply(1:length(Censored), function(k) Censored[[k]][current_train_ind[[k]]])
  
  sink(paste(currentwd, "PosteriorLikelihoodForHierarchicalSS", ".txt", sep = ""), append=TRUE)
  cat(paste("iteration",cv_iter, "has started at time", Sys.time(), "\n"))
  
  # ---------------------------------------------------------------------------
  # Hierarchical Spike/Slab model
  posteriors_hierarchical_ss = HierarchicalLogNormalSpikeSlab(Covariates_Training, Survival_Training, Censored_Training, 
                                                              starting_values, iters, covariates_in_model, 
                                                              covariates_by_cancer, priors, pi_generation = "shared_across_cancers",
                                                              progress = "model_comparison")
  
  cat(paste("hierarchical s/s model in iteration",cv_iter, "has finished at time", Sys.time(), "\n"))
  
  PosteriorLikelihood(posteriors_hierarchical_ss, Covariates, Survival, Censored, current_train_ind, iters/2)
}
stopCluster(cl)
end <- Sys.time()
end-start

# -----------------------------------------------------------------------------
# Computing the out-of-sample posterior likelihood for the separated model
# -----------------------------------------------------------------------------

writeLines(c(""), paste(currentwd, "PosteriorLikelihoodForSeparatedModel", ".txt", sep = ""))

iters = 1e5
start <- Sys.time()
cl <- makeCluster(5)
registerDoParallel(cl)
PL.Separated = foreach(cv_iter = 1:5, .packages = c("MASS", "truncnorm", "EnvStats", "Matrix")) %dopar% {
  current_train_ind = sapply(TrainingSet[[1]], '[', cv_iter)
  Covariates_Training = lapply(1:length(Covariates), function(k) Covariates[[k]][current_train_ind[[k]], ])
  Survival_Training = lapply(1:length(Survival), function(k) Survival[[k]][current_train_ind[[k]]])
  Censored_Training = lapply(1:length(Censored), function(k) Censored[[k]][current_train_ind[[k]]])
  
  sink(paste(currentwd, "PosteriorLikelihoodForSeparatedModel", ".txt", sep = ""), append=TRUE)
  cat(paste("iteration",cv_iter, "has started at time", Sys.time(), "\n"))
  
  # ---------------------------------------------------------------------------
  # Hierarchical Spike/Slab model
  PosteriorsByCancerSeparately <- lapply(1:n_cancer, function(i) list()) # to store the posteriors for each cancer type
  names(PosteriorsByCancerSeparately) <- cancer_types
  
  for (type in  1:n_cancer) { # for each cancer type
    # Subsetting just the data for a single cancer type
    CovariatesTraining_type <- Covariates_Training[type]
    SurvivalTraining_type <- Survival_Training[type]
    CensoredTraining_type <- Censored_Training[type]
    
    # Adjusting these parameters to reflect just a single cancer type
    covariates_by_cancer_type <- covariates_by_cancer_full[type]
    covariates_in_model_type <- unlist(covariates_by_cancer_type)
    
    # Selecting just the starting values that correspond to this cancer type
    avail <- covariates_in_model_full %in% covariates_by_cancer_type[[1]]
    starting_values_type <- list(sigma2_start = starting_values$sigma2_start,
                                 beta_tilde_start = starting_values$beta_tilde_start[avail],
                                 lambda2_start = starting_values$lambda2_start[avail],
                                 gamma_start = list(list(starting_values$gamma_start[[type]][[1]])),
                                 pi_start = starting_values$pi_start[avail])
    
    # Running the model
    posteriors_type <- HierarchicalLogNormalSpikeSlab(CovariatesTraining_type, SurvivalTraining_type, CensoredTraining_type,
                                                      starting_values_type, iters, covariates_in_model_type, covariates_by_cancer_type,
                                                      priors, pi_generation = "fixed_at_0.5",
                                                      progress = "model_comparison")
    # Storing the results
    PosteriorsByCancerSeparately[[type]] <- posteriors_type
  }
  
  cat(paste("the separated model in parallel process",cv_iter, "has finished at time", Sys.time(), "\n"))
  
  # Computing the posterior likelihood for the separated model
  joint_ll_for_each_cancer <- lapply(1:n_cancer, function(type) {
    JointLogLikelihoodForSeparateModel(posteriors = PosteriorsByCancerSeparately[[type]],
                                       X = Covariates[type],
                                       Y = Survival[type],
                                       Censored = Censored[type], 
                                       current_train_ind[type], 
                                       burnin = iters/2)
  })
  
  # Return the posterior likelihood
  iters_used_in_ll <- (iters/2/10) + 1
  -log(iters_used_in_ll) + logSum(Reduce('+', joint_ll_for_each_cancer))
}
stopCluster(cl)
end <- Sys.time()
end-start

