# Calculating the 5-fold cross validation log posterior likelihood
# for each of the models we want to compare. 
# Authors: Sarah Samorodnitsky and Eric Lock (2020)
# University of Minnesota

library(foreach)
library(doParallel)

# Saving important working directories
currentwd = "~/PanTCGA/SpikeAndSlabResults/SpikeAndSlabGibbsSampler/ModelWithInterceptOutsideSS/ResultsAfterFiltering_V3/ModelComparison/" # change this to your directory

# load in the Covariates, Survival, Censored data
load("XYC_V2_WithAge_StandardizedPredictors.rda")

# load in helper functions
source("HelperFunctions.R")

# load in model with spike and slab and intercept outside s/s framework
source("ExtendedHierarchicalModelLogNormalSpikeAndSlabInterceptOutsideSS.R")

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
writeLines(c(""), paste(currentwd, "PosteriorLikelihoodForAllModelsLog2", ".txt", sep = ""))

start <- Sys.time()
cl <- makeCluster(5)
registerDoParallel(cl)
PL = foreach(cv_iter = 1:5, .packages = c("MASS", "truncnorm", "EnvStats")) %dopar% {
  current_train_ind = sapply(TrainingSet[[1]], '[', cv_iter)
  Covariates_Training = lapply(1:length(Covariates), function(k) Covariates[[k]][current_train_ind[[k]], ])
  Survival_Training = lapply(1:length(Survival), function(k) Survival[[k]][current_train_ind[[k]]])
  Censored_Training = lapply(1:length(Censored), function(k) Censored[[k]][current_train_ind[[k]]])
  
  sink(paste(currentwd, "PosteriorLikelihoodForAllModelsLog2", ".txt", sep = ""), append=TRUE)
  cat(paste("iteration",cv_iter, "has started at time", Sys.time(), "\n"))
  
  # Hierarchical Spike/Slab model
  set.seed(2) # setting a seed so I get the same result each time
  posteriors_hierarchical_ss = HierarchicalLogNormalSpikeSlab(Covariates_Training, Survival_Training, Censored_Training, 
                                                              starting_values, iters, covariates_in_model, 
                                                              covariates_by_cancer, priors, pi_generation = "shared_across_cancers")
  
  cat(paste("hierarchical s/s model in parallel process",cv_iter, "has finished at time", Sys.time(), "\n"))
  
  # Model with prior inclusion probability fixed at 0.5
  set.seed(2) # setting a seed so I get the same result each time
  posteriors_fixed_0.5 = HierarchicalLogNormalSpikeSlab(Covariates_Training, Survival_Training, Censored_Training, 
                                                        starting_values_fixed_at_0.5, iters, covariates_in_model, 
                                                        covariates_by_cancer, priors, pi_generation = "fixed_at_0.5")
  
  cat(paste("fixed at 0.5 model in parallel process",cv_iter, "has finished at time", Sys.time(), "\n"))
  
  # Model with prior inclusion probability fixed at 1.0
  set.seed(2) # setting a seed so I get the same result each time
  posteriors_fixed_1.0 = HierarchicalLogNormalSpikeSlab(Covariates_Training, Survival_Training, Censored_Training, 
                                                        starting_values_fixed_at_1.0, iters, covariates_in_model, 
                                                        covariates_by_cancer, priors, pi_generation = "fixed_at_1.0")
  
  cat(paste("fixed at 1.0 model in parallel process",cv_iter, "has finished at time", Sys.time(), "\n"))
  
  # Model where prior inclusion probability is the same for all betas
  set.seed(2) # setting a seed so I get the same result each time
  posteriors_shared_across_betas = HierarchicalLogNormalSpikeSlab(Covariates_Training, Survival_Training, Censored_Training, 
                                                                  starting_values_shared_across_betas, iters, covariates_in_model, 
                                                                  covariates_by_cancer, priors, pi_generation = "shared_across_betas")
  
  cat(paste("shared across betas model in parallel process",cv_iter, "has finished at time", Sys.time(), "\n"))
  
  # Null model
  set.seed(2) # setting a seed so I get the same result each time
  Covariates_Training_Null = lapply(1:length(Covariates_null), function(k) Covariates_null[[k]][current_train_ind[[k]], , drop = FALSE])
  
  posteriors_null = HierarchicalLogNormalSpikeSlab(Covariates_Training_Null, Survival_Training, Censored_Training, 
                                                   starting_values_null, iters, covariates_in_model_null, 
                                                   covariates_by_cancer_null, priors, pi_generation = "fixed_at_1.0")
  
  cat(paste("null (intercept only) model in parallel process",cv_iter, "has finished at time", Sys.time(), "\n"))
  
  # Computing the posterior likelihoods
  Results = c(PosteriorLikelihood(posteriors_hierarchical_ss, Covariates, Survival, Censored, current_train_ind, iters/2),
              PosteriorLikelihood(posteriors_fixed_0.5, Covariates, Survival, Censored, current_train_ind, iters/2),
              PosteriorLikelihood(posteriors_fixed_1.0, Covariates, Survival, Censored, current_train_ind, burnin = iters/2),
              PosteriorLikelihood(posteriors_shared_across_betas, Covariates, Survival, Censored, current_train_ind, iters/2),
              PosteriorLikelihood(posteriors_null, Covariates_null, Survival, Censored, current_train_ind, burnin = iters/2))
  
  names(Results) = c("Hierarchical S/S", "Fixed at 0.5", "Fixed at 1.0", "Shared Across Betas", "Null Model")
  
  Results
}
stopCluster(cl)
end <- Sys.time()
end-start

# Computing the results
PL.Matrix = do.call(rbind, PL)
PL.Means = colMeans(PL.Matrix)
