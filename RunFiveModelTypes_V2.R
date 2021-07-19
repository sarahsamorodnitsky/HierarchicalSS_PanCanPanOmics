# This function runs each of the five model types in serial to be used for the 
# model simulations. 

library(PanCanVarSel)

modelwd = "~/PanTCGA/SpikeAndSlabResults/SpikeAndSlabGibbsSampler/ModelWithInterceptOutsideSS/"

# load in model with spike and slab and intercept outside s/s framework
source(paste(modelwd, "ExtendedHierarchicalModelLogNormalSpikeAndSlabInterceptOutsideSS_V3.R", sep = ""))

RunFiveModelTypes = function(Covariates, Survival, Censored, starting_values_for_every_model, 
                             iters, covariates_in_model,
                             covariates_by_cancer, priors) {
  # starting_values_for_every_model (list): provides list of starting value lists for each model type
  # iters (int): number of iterations to run every model for
  # priors (list): gives prior values for every parameter
  
  # Creating the starting_values variables
  starting_values = starting_values_for_every_model$starting_values # for the hierarchical s/s model
  starting_values_fixed_at_0.5 = starting_values_for_every_model$starting_values_fixed_at_0.5
  starting_values_fixed_at_1.0 = starting_values_for_every_model$starting_values_fixed_at_1.0
  starting_values_shared_across_betas = starting_values_for_every_model$starting_values_shared_across_betas
  starting_values_null = starting_values_for_every_model$starting_values_null
  starting_values_joint = starting_values_for_every_model$starting_values_joint
  starting_values_separate = starting_values_for_every_model$starting_values_separate
  
  # Running each of the models ------------------------------------------------
  # Hierarchical s/s ----------------------------------------------------------
  posteriors_hierarchical_ss = HierarchicalLogNormalSpikeSlab(Covariates, Survival, Censored, 
                                                              starting_values, iters, covariates_in_model, 
                                                              covariates_by_cancer, priors, 
                                                              pi_generation = "shared_across_cancers",
                                                              progress = "simulation")
  
  # Fixed at 0.5 --------------------------------------------------------------
  posteriors_fixed_0.5 = HierarchicalLogNormalSpikeSlab(Covariates, Survival, Censored, 
                                                        starting_values_fixed_at_0.5, iters, covariates_in_model, 
                                                        covariates_by_cancer, priors, 
                                                        pi_generation = "fixed_at_0.5",
                                                        progress = "simulation")
  
  # Fixed at 1.0 --------------------------------------------------------------
  posteriors_fixed_1.0 = HierarchicalLogNormalSpikeSlab(Covariates, Survival, Censored, 
                                                        starting_values_fixed_at_1.0, iters, covariates_in_model, 
                                                        covariates_by_cancer, priors, 
                                                        pi_generation = "fixed_at_1.0",
                                                        progress = "simulation")
  
  # Shared across betas -------------------------------------------------------
  posteriors_shared_across_betas = HierarchicalLogNormalSpikeSlab(Covariates, Survival, Censored, 
                                                                  starting_values_shared_across_betas, iters, covariates_in_model, 
                                                                  covariates_by_cancer, priors, 
                                                                  pi_generation = "shared_across_betas",
                                                                  progress = "simulation")
  
  # Null model ----------------------------------------------------------------
  Covariates_Null = lapply(Covariates, function(type) type[, 1, drop = FALSE])
  covariates_by_cancer_null = lapply(covariates_by_cancer, function(type) type[1]) # just storing the intercept for every cancer
  covariates_in_model_null = covariates_in_model[1]
  
  posteriors_null = HierarchicalLogNormalSpikeSlab(Covariates_Null, Survival, Censored, 
                                                   starting_values_null, iters, covariates_in_model_null, 
                                                   covariates_by_cancer_null, priors, 
                                                   pi_generation = "fixed_at_1.0",
                                                   progress = "simulation")
  
  # Joint model ---------------------------------------------------------------
  CovariatesCombined <- list(do.call(rbind, CombineCancerDataTogether(Covariates, covariates_in_model_full)))
  SurvivalCombined <- list(unlist(Survival))
  CensoredCombined <- list(unlist(Censored))
  covariates_by_cancer_joint <- list(covariates_in_model) # since all the cancers are combined and have 0s filled in, I am just creating a new covariates_by_cancer specific to this scenario.
  
  posteriors_joint = HierarchicalLogNormalSpikeSlab(CovariatesCombined, SurvivalCombined, CensoredCombined,
                                                    starting_values_joint, iters, covariates_in_model, covariates_by_cancer_joint,
                                                    priors, pi_generation = "fixed_at_0.5",
                                                    progress = "simulation")
  
  
  # Separated model -----------------------------------------------------------
  posteriors_by_each_cancer_separately <- lapply(1:n_cancer, function(i) list()) # to store the posteriors for each cancer type
  names(posteriors_by_each_cancer_separately) <- cancer_types
  
  for (type in  1:n_cancer) { # for each cancer type
    # Subsetting just the data for a single cancer type
    Covariates_type <- Covariates[type]
    Survival_type <- Survival[type]
    Censored_type <- Censored[type]
    
    # Adjusting these parameters to reflect just a single cancer type
    covariates_by_cancer_type <- covariates_by_cancer[type]
    covariates_in_model_type <- unlist(covariates_by_cancer_type)
    
    # Selecting just the starting values that correspond to this cancer type
    starting_values_type <- starting_values_separate[[type]]
    
    # Running the model
    posteriors_type <- HierarchicalLogNormalSpikeSlab(Covariates_type, Survival_type, Censored_type,
                                                      starting_values_type, iters, covariates_in_model_type, covariates_by_cancer_type,
                                                      priors, pi_generation = "fixed_at_0.5",
                                                      progress = "simulation")
    # Storing the results
    posteriors_by_each_cancer_separately[[type]] <- posteriors_type
  }
  
  
  # Maity et al. model --------------------------------------------------------
  CovariatesCombinedPCVS <- PrepareCancerDataForPanCanVarSel(Covariates, covariates_in_model)
  SurvivalCombined <- unlist(Survival) 
  CensoredCombined <- unlist(Censored)
  n_vec_training <- sapply(Covariates, nrow)
  
  # Putting the survival and censor times together with a right-censor indicator
  SurvivalCensored <- SurvivalCombined
  SurvivalCensored[is.na(SurvivalCensored)] <- CensoredCombined[is.na(SurvivalCombined)] # filling in the censored times with the right-censor times
  ct <- as.matrix(cbind(time = SurvivalCensored, 
                        censored = !is.na(unlist(SurvivalCombined)))) # 1 = had event, 0 = censored (according to package documentation)
  
  posteriors_horseshoe <- hsaftgroupcorr(ct = ct, X = CovariatesCombinedPCVS, method.tau = "truncatedCauchy",
                                         method.sigma = "Jeffreys", r = n_cancer, alpha = 0.05, burn = iters/2,
                                         nmc = iters/2, n.seq = n_vec_training, pk = p)
  
  # Combining all results into one list
  PosteriorsForEveryModel = list(hierarchical_ss = posteriors_hierarchical_ss,
                                 fixed_0.5 = posteriors_fixed_0.5,
                                 fixed_1.0 = posteriors_fixed_1.0,
                                 shared_across_betas = posteriors_shared_across_betas,
                                 null_model = posteriors_null,
                                 joint_model = posteriors_joint,
                                 separate_model = posteriors_by_each_cancer_separately,
                                 horseshoe_model = posteriors_horseshoe)
  
  return(PosteriorsForEveryModel)
}
