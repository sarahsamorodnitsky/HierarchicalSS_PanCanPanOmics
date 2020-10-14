# This function runs each of the five model types in serial to be used for the 
# model simulations. 

modelwd = "~/PanTCGA/SpikeAndSlabResults/SpikeAndSlabGibbsSampler/ModelWithInterceptOutsideSS/"

# load in model with spike and slab and intercept outside s/s framework
source(paste(modelwd, "ExtendedHierarchicalModelLogNormalSpikeAndSlabInterceptOutsideSS.R", sep = ""))

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
  
  # Running each of the models
  # Hierarchical s/s
  posteriors_hierarchical_ss = HierarchicalLogNormalSpikeSlab(Covariates, Survival, Censored, 
                                                              starting_values, iters, covariates_in_model, 
                                                              covariates_by_cancer, priors, 
                                                              pi_generation = "shared_across_cancers",
                                                              progress = "simulation")
  
  # Fixed at 0.5
  posteriors_fixed_0.5 = HierarchicalLogNormalSpikeSlab(Covariates, Survival, Censored, 
                                                        starting_values_fixed_at_0.5, iters, covariates_in_model, 
                                                        covariates_by_cancer, priors, 
                                                        pi_generation = "fixed_at_0.5",
                                                        progress = "simulation")
  
  # Fixed at 1.0
  posteriors_fixed_1.0 = HierarchicalLogNormalSpikeSlab(Covariates, Survival, Censored, 
                                                        starting_values_fixed_at_1.0, iters, covariates_in_model, 
                                                        covariates_by_cancer, priors, 
                                                        pi_generation = "fixed_at_1.0",
                                                        progress = "simulation")
  
  # Shared across betas
  posteriors_shared_across_betas = HierarchicalLogNormalSpikeSlab(Covariates, Survival, Censored, 
                                                                  starting_values_shared_across_betas, iters, covariates_in_model, 
                                                                  covariates_by_cancer, priors, 
                                                                  pi_generation = "shared_across_betas",
                                                                  progress = "simulation")
  
  # Null model
  Covariates_null = lapply(Covariates, function(type) type[, 1, drop = FALSE])
  covariates_by_cancer_null = lapply(covariates_by_cancer, function(type) type[1]) # just storing the intercept for every cancer
  covariates_in_model_null = covariates_in_model[1]
  
  posteriors_null = HierarchicalLogNormalSpikeSlab(Covariates_Null, Survival, Censored, 
                                                   starting_values_null, iters, covariates_in_model_null, 
                                                   covariates_by_cancer_null, priors, 
                                                   pi_generation = "fixed_at_1.0",
                                                   progress = "simulation")
  
  # Combining all results into one list
  PosteriorsForEveryModel = list(hierarchical_ss = posteriors_hierarchical_ss,
                                 fixed_0.5 = posteriors_fixed_0.5,
                                 fixed_1.0 = posteriors_fixed_1.0,
                                 shared_across_betas = posteriors_shared_across_betas,
                                 null_model = posteriors_null)
  
  return(PosteriorsForEveryModel)
}
