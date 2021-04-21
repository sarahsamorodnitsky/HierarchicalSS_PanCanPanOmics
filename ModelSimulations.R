# Simulations to compare how each of our 5 proposed models perform under
# different covariate conditions. Our 5 models are the hierarchical spike/slab,
# fixed at 0.5, shared across all predictors, full model, and the null model.

# Our 5 conditions under which we will assess each method are as follows:
# (1) Half of predictors are present across *all* cancer types that have it, 
#     half of predictors not present at all. eg. if 5 cancers have covariate 2,
#     then if covariate is either included or excluded for all 5 cancers.
# 
# (2) 0.1 of predictors are present across *all* cancer types that have it, 
#     0.9 of predictors not present at all. eg. if 5 cancers have covariate 2,
#     then if covariate is either included or excluded for all 5 cancers with 
#     probability 0.1
#
# (3) Each cancer type+predictor pair independently has a 0.5 chance of being 
#     included. eg. if 5 cancers have covariate 2, then covariate 2 is independently
#     included for each cancer, so it might be for cancer A but not cancer B.
#
# (4) Each cancer type+predictor pair has a 0.1 chance independently of being included.
#     Same thing as before but lower probability of inclusion
#
# (5) "Full model" -- all covariates are truly included 
#
# (6) "Null model" -- all covariates are truly not included

# For each scenario, run each of the 5 models with 10 replications. Assess each in terms of 
# squared deviation as defined in our notes and in terms of the posterior likelihood. 

# When assessing the squared deviation: run each model for 2000 iterations (may need to increase this), 
# 1000 iteration burn-in. Then, compute the squared deviation. Average across all 10 replications.

# When assessing the log-posterior likelihood, run the model on training data for 2000 iterations
# with 1000 iteration burn-in. Then, generate a test data set of the same size and using the results
# from the training data set run, compute the log-posterior likelihood. Average across 10 replications. 

library(foreach)
library(doParallel)

currentwd = "~/PanTCGA/SpikeAndSlabResults/SpikeAndSlabGibbsSampler/ModelWithInterceptOutsideSS/ResultsAfterFiltering_V3/ModelSimulations/"
PanTCGAwd = "~/PanTCGA/"
datawd = "~/PanTCGA/SpikeAndSlabResults/SpikeAndSlabGibbsSampler/ModelWithInterceptOutsideSS/ResultsAfterFiltering_V3/"

# load in the Covariates, Survival, Censored data to generate fake data with the same structure
load(paste(datawd, "XYC_V2_WithAge_StandardizedPredictors.rda", sep = ""))

# load in helper functions
source(paste(PanTCGAwd, "HelperFunctions.R", sep = ""))

# load in function that runs the model for each of the five types
source(paste(currentwd, "RunFiveModelTypes.R", sep = ""))


# Important variables
# Using these to generate fake data and also run the model
covariates_by_cancer_full = lapply(Covariates, function(cancer) as.numeric(colnames(cancer))) # the covariates each cancer type has (each should have a 1 for the intercept)
covariates_table = table(unlist(covariates_by_cancer_full)) # gives the number of cancer types who have data on each PC
covariates_in_model_full = as.numeric(names(covariates_table))

covariates_in_model = covariates_in_model_full
p = length(covariates_in_model)
covariates_by_cancer = covariates_by_cancer_full # the name of this variable in the HelperFunctions.R file
iters = 10000
burnin = iters/2
n_cancer = length(Covariates)

cancer_types = names(Covariates)

# Variables for generating fake data and for running the model
# Creating the priors argument
priors = list(betatilde_priorvar_intercept = 10^2,    # \tilde\beta_0 ~ N(0,100). both beta tildes are centered at 0
              betatilde_priorvar_coefficient = 1,     # -> \tilde\beta_p ~ N(0,1), p > 0
              lambda2_priorshape_intercept = 1,       # -> \lambda^2_0 ~ IG(1,1)
              lambda2_priorrate_intercept = 1,         
              lambda2_priorshape_coefficient = 5,     # -> \lambda^2_p ~ IG(5,1) p > 0
              lambda2_priorrate_coefficient = 1,
              sigma2_priorshape = 1,                # -> \sigma^2 ~ IG(1, 1) changed this to hopefully help with the simulation
              sigma2_priorrate = 1,
              spike_priorvar = 1/10000)

# Generating starting values - these will be the same for every replication and model type
starting_values_for_every_model = GenerateStartingValues(p, covariates_in_model, covariates_by_cancer)

# Running the 10 replications in parallel
writeLines(c(""), paste(currentwd, "ModelSimulationsLog", ".txt", sep = ""))

start <- Sys.time()
cl <- makeCluster(10)
registerDoParallel(cl)
# Changed from FiveModelsFiveDatasetsTenReps -> FiveModelsFiveDatasetsThirtyReps for 2nd revision. 
FiveModelsFiveDatasetsThirtyReps = foreach(repli = 1:30, .packages = c("MASS", "truncnorm", "EnvStats")) %dopar% {
  set.seed(repli) # set a seed so for each replication the same data and true values are generated
  
  sink(paste(currentwd, "ModelSimulationsLog", ".txt", sep = ""), append=TRUE)
  cat(paste("iteration", repli, "has started at time", Sys.time(), "\n"))
  
  # (1) Half of predictors are present across *all* cancer types that have it, 
  #     half of predictors not present at all. eg. if 5 cancers have covariate 2,
  #     then if covariate is either included or excluded for all 5 cancers.
  
  # Generating the true values under the current condition
  TrueValues_C1 = GenerateTrueValues(p, covariates_by_cancer, covariates_in_model, 
                                     n_cancer, priors, condition = "half_present_across_all")
  true_betas_C1 = TrueValues_C1$true_betas
  true_sigma2_C1 = TrueValues_C1$true_sigma2
  true_beta_tilde_C1 = TrueValues_C1$true_beta_tilde
  true_lambda2_C1 = TrueValues_C1$true_lambda2
  true_gammas_C1 = TrueValues_C1$true_gammas
  
  # Generating fake data
  SimulatedData_C1 = DataBasedOnTrueValues(true_betas_C1, true_sigma2_C1, n_vec, covariates_by_cancer)
  Design_C1 = SimulatedData_C1$Design
  Response_C1 = SimulatedData_C1$Response
  Censored_C1 = SimulatedData_C1$Censored
  
  # for testing
  # Covariates = Design_C1
  # Survival = Response_C1
  # Censored = Censored_C1
  
  PosteriorsForEveryModel_C1 = RunFiveModelTypes(Design_C1, Response_C1, Censored_C1, 
                                                 starting_values_for_every_model,
                                                 iters, covariates_in_model, 
                                                 covariates_by_cancer, priors)
  
  # Computing the sum of squared deviations
  model_types = names(PosteriorsForEveryModel_C1)
  SumOfSquaredDevs_C1 = ComputeSumOfSquaredDeviationsForEveryModel(true_gammas_C1, PosteriorsForEveryModel_C1, iters, model_types)
  
  # Computing the posterior likelihood based on test data
  SimulatedTestData_C1 = DataBasedOnTrueValues(true_betas_C1, true_sigma2_C1, n_vec, covariates_by_cancer)
  Design_C1_Test = SimulatedTestData_C1$Design
  Response_C1_Test = SimulatedTestData_C1$Response
  Censored_C1_Test = SimulatedTestData_C1$Censored
  
  # Calculating the posterior likelihood using the test data
  PLs_C1 = ComputePosteriorLikelihoodForEveryModel(PosteriorsForEveryModel_C1, Design_C1_Test, Response_C1_Test,
                                          Censored_C1_Test, iters, burnin, model_types)
  
  cat(paste("model runs for condition 1 in parallel process", repli, "has finished at time", Sys.time(), "\n"))
  
  # (2) 0.1 of predictors are present across *all* cancer types that have it, 
  #     0.9 of predictors not present at all. eg. if 5 cancers have covariate 2,
  #     then if covariate is either included or excluded for all 5 cancers with 
  #     probability 0.1
  
  # Generating the true values under the current condition
  TrueValues_C2 = GenerateTrueValues(p, covariates_by_cancer, covariates_in_model, 
                                     n_cancer, priors, condition = "0.1_present_across_all")
  true_betas_C2 = TrueValues_C2$true_betas
  true_sigma2_C2 = TrueValues_C2$true_sigma2
  true_beta_tilde_C2 = TrueValues_C2$true_beta_tilde
  true_lambda2_C2 = TrueValues_C2$true_lambda2
  true_gammas_C2 = TrueValues_C2$true_gammas
  
  # Generating fake data
  SimulatedData_C2 = DataBasedOnTrueValues(true_betas_C2, true_sigma2_C2, n_vec, covariates_by_cancer)
  Design_C2 = SimulatedData_C2$Design
  Response_C2 = SimulatedData_C2$Response
  Censored_C2 = SimulatedData_C2$Censored
  
  # for testing
  # Covariates = Design_C2
  # Survival = Response_C2
  # Censored = Censored_C2
  
  PosteriorsForEveryModel_C2 = RunFiveModelTypes(Design_C2, Response_C2, Censored_C2, 
                                                 starting_values_for_every_model,
                                                 iters, covariates_in_model, 
                                                 covariates_by_cancer, priors)
  
  # Computing the sum of squared deviations
  model_types = names(PosteriorsForEveryModel_C2)
  SumOfSquaredDevs_C2 = ComputeSumOfSquaredDeviationsForEveryModel(true_gammas_C2, PosteriorsForEveryModel_C2, iters, model_types)
  
  # Computing the posterior likelihood based on test data
  SimulatedTestData_C2 = DataBasedOnTrueValues(true_betas_C2, true_sigma2_C2, n_vec, covariates_by_cancer)
  Design_C2_Test = SimulatedTestData_C2$Design
  Response_C2_Test = SimulatedTestData_C2$Response
  Censored_C2_Test = SimulatedTestData_C2$Censored
  
  # Calculating the posterior likelihood using the test data
  PLs_C2 = ComputePosteriorLikelihoodForEveryModel(PosteriorsForEveryModel_C2, Design_C2_Test, Response_C2_Test,
                                                   Censored_C2_Test, iters, burnin, model_types)
  
  cat(paste("model runs for condition 2 in parallel process", repli, "has finished at time", Sys.time(), "\n"))
  
  # (3) Each cancer type+predictor pair independently has a 0.5 chance of being 
  #     included. eg. if 5 cancers have covariate 2, then covariate 2 is independently
  #     included for each cancer, so it might be for cancer A but not cancer B.
  
  TrueValues_C3 = GenerateTrueValues(p, covariates_by_cancer, covariates_in_model, 
                                     n_cancer, priors, condition = "independent_0.5")
  true_betas_C3 = TrueValues_C3$true_betas
  true_sigma2_C3 = TrueValues_C3$true_sigma2
  true_beta_tilde_C3 = TrueValues_C3$true_beta_tilde
  true_lambda2_C3 = TrueValues_C3$true_lambda2
  true_gammas_C3 = TrueValues_C3$true_gammas
  
  # Generating fake data
  SimulatedData_C3 = DataBasedOnTrueValues(true_betas_C3, true_sigma2_C3, n_vec, covariates_by_cancer)
  Design_C3 = SimulatedData_C3$Design
  Response_C3 = SimulatedData_C3$Response
  Censored_C3 = SimulatedData_C3$Censored
  
  # for testing
  # Covariates = Design_C1
  # Survival = Response_C1
  # Censored = Censored_C1
  
  PosteriorsForEveryModel_C3 = RunFiveModelTypes(Design_C3, Response_C3, Censored_C3, starting_values_for_every_model,
                                                 iters, covariates_in_model, covariates_by_cancer, priors)
  
  # Computing the sum of squared deviations
  model_types = names(PosteriorsForEveryModel_C3)
  SumOfSquaredDevs_C3 = ComputeSumOfSquaredDeviationsForEveryModel(true_gammas_C3, PosteriorsForEveryModel_C3, iters, model_types)
  
  # Computing the posterior likelihood based on test data
  SimulatedTestData_C3 = DataBasedOnTrueValues(true_betas_C3, true_sigma2_C3, n_vec, covariates_by_cancer)
  Design_C3_Test = SimulatedTestData_C3$Design
  Response_C3_Test = SimulatedTestData_C3$Response
  Censored_C3_Test = SimulatedTestData_C3$Censored
  
  # Calculating the posterior likelihood on the test data
  PLs_C3 = ComputePosteriorLikelihoodForEveryModel(PosteriorsForEveryModel_C3, Design_C3_Test, Response_C3_Test,
                                                   Censored_C3_Test, iters, burnin, model_types)
  
  cat(paste("model runs for condition 3 in parallel process", repli, "has finished at time", Sys.time(), "\n"))
  
  # (4) Each cancer type+predictor pair has a 0.1 chance independently of being included.
  #     Same thing as before but lower probability of inclusion
  
  TrueValues_C4 = GenerateTrueValues(p, covariates_by_cancer, covariates_in_model, 
                                     n_cancer, priors, condition = "independent_0.1")
  true_betas_C4 = TrueValues_C4$true_betas
  true_sigma2_C4 = TrueValues_C4$true_sigma2
  true_beta_tilde_C4 = TrueValues_C4$true_beta_tilde
  true_lambda2_C4 = TrueValues_C4$true_lambda2
  true_gammas_C4 = TrueValues_C4$true_gammas
  
  # Generating fake data
  SimulatedData_C4 = DataBasedOnTrueValues(true_betas_C4, true_sigma2_C4, n_vec, covariates_by_cancer)
  Design_C4 = SimulatedData_C4$Design
  Response_C4 = SimulatedData_C4$Response
  Censored_C4 = SimulatedData_C4$Censored
  
  # for testing
  # Covariates = Design_C4
  # Survival = Response_C4
  # Censored = Censored_C4
  
  PosteriorsForEveryModel_C4 = RunFiveModelTypes(Design_C4, Response_C4, Censored_C4, starting_values_for_every_model,
                                                 iters, covariates_in_model, covariates_by_cancer, priors)
  
  # Computing the sum of squared deviations
  model_types = names(PosteriorsForEveryModel_C4)
  SumOfSquaredDevs_C4 = ComputeSumOfSquaredDeviationsForEveryModel(true_gammas_C4, PosteriorsForEveryModel_C4, iters, model_types)
  
  # Computing the posterior likelihood based on test data
  SimulatedTestData_C4 = DataBasedOnTrueValues(true_betas_C4, true_sigma2_C4, n_vec, covariates_by_cancer)
  Design_C4_Test = SimulatedTestData_C4$Design
  Response_C4_Test = SimulatedTestData_C4$Response
  Censored_C4_Test = SimulatedTestData_C4$Censored
  
  # Calculating the posterior likelihood based on the test data
  PLs_C4 = ComputePosteriorLikelihoodForEveryModel(PosteriorsForEveryModel_C4, Design_C4_Test, Response_C4_Test,
                                                   Censored_C4_Test, iters, burnin, model_types)
  
  cat(paste("model runs for condition 4 in parallel process", repli, "has finished at time", Sys.time(), "\n"))
  
  # (5) "Full model" -- all covariates are truly included 
  
  TrueValues_C5 = GenerateTrueValues(p, covariates_by_cancer, covariates_in_model, 
                                     n_cancer, priors, condition = "full_all_included")
  true_betas_C5 = TrueValues_C5$true_betas
  true_sigma2_C5 = TrueValues_C5$true_sigma2
  true_beta_tilde_C5 = TrueValues_C5$true_beta_tilde
  true_lambda2_C5 = TrueValues_C5$true_lambda2
  true_gammas_C5 = TrueValues_C5$true_gammas
  
  # Generating fake data
  SimulatedData_C5 = DataBasedOnTrueValues(true_betas_C5, true_sigma2_C5, n_vec, covariates_by_cancer)
  Design_C5 = SimulatedData_C5$Design
  Response_C5 = SimulatedData_C5$Response
  Censored_C5 = SimulatedData_C5$Censored
  
  # for testing
  # Covariates = Design_C5
  # Survival = Response_C5
  # Censored = Censored_C5
  
  PosteriorsForEveryModel_C5 = RunFiveModelTypes(Design_C5, Response_C5, Censored_C5, starting_values_for_every_model,
                                                 iters, covariates_in_model, covariates_by_cancer, priors)
  
  # Computing the sum of squared deviations
  model_types = names(PosteriorsForEveryModel_C5)
  SumOfSquaredDevs_C5 = ComputeSumOfSquaredDeviationsForEveryModel(true_gammas_C5, PosteriorsForEveryModel_C5, iters, model_types)
  
  # Computing the posterior likelihood based on test data
  SimulatedTestData_C5 = DataBasedOnTrueValues(true_betas_C5, true_sigma2_C5, n_vec, covariates_by_cancer)
  Design_C5_Test = SimulatedTestData_C5$Design
  Response_C5_Test = SimulatedTestData_C5$Response
  Censored_C5_Test = SimulatedTestData_C5$Censored
  
  # Generating posterior based on the test data
  PLs_C5 = ComputePosteriorLikelihoodForEveryModel(PosteriorsForEveryModel_C5, Design_C5_Test, Response_C5_Test,
                                                   Censored_C5_Test, iters, burnin, model_types)
  
  cat(paste("model runs for condition 5 in parallel process", repli, "has finished at time", Sys.time(), "\n"))
  
  # (6) "Null model" -- all covariates are truly not included
  TrueValues_C6 = GenerateTrueValues(p, covariates_by_cancer, covariates_in_model, 
                                     n_cancer, priors, condition = "null_none_included")
  true_betas_C6 = TrueValues_C6$true_betas
  true_sigma2_C6 = TrueValues_C6$true_sigma2
  true_beta_tilde_C6 = TrueValues_C6$true_beta_tilde
  true_lambda2_C6 = TrueValues_C6$true_lambda2
  true_gammas_C6 = TrueValues_C6$true_gammas
  
  # Generating fake data
  SimulatedData_C6 = DataBasedOnTrueValues(true_betas_C6, true_sigma2_C6, n_vec, covariates_by_cancer)
  Design_C6 = SimulatedData_C6$Design
  Response_C6 = SimulatedData_C6$Response
  Censored_C6 = SimulatedData_C6$Censored
  
  # for testing
  # Covariates = Design_C6
  # Survival = Response_C6
  # Censored = Censored_C6
  
  PosteriorsForEveryModel_C6 = RunFiveModelTypes(Design_C6, Response_C6, Censored_C6, starting_values_for_every_model,
                                                 iters, covariates_in_model, covariates_by_cancer, priors)
  
  # Computing the sum of squared deviations
  model_types = names(PosteriorsForEveryModel_C6)
  SumOfSquaredDevs_C6 = ComputeSumOfSquaredDeviationsForEveryModel(true_gammas_C6, PosteriorsForEveryModel_C6, iters, model_types)
  
  # Computing the posterior likelihood based on test data
  SimulatedTestData_C6 = DataBasedOnTrueValues(true_betas_C6, true_sigma2_C6, n_vec, covariates_by_cancer)
  Design_C6_Test = SimulatedTestData_C6$Design
  Response_C6_Test = SimulatedTestData_C6$Response
  Censored_C6_Test = SimulatedTestData_C6$Censored
  
  # Calculating the posterior likelihood based on the test data
  PLs_C6 = ComputePosteriorLikelihoodForEveryModel(PosteriorsForEveryModel_C6, Design_C6_Test, Response_C6_Test,
                                                   Censored_C6_Test, iters, burnin, model_types)
  
  cat(paste("model runs for condition 6 in parallel process", repli, "has finished at time", Sys.time(), "\n"))
  
  # Appending the results from each of the 5 scenarios together
  SumOfSquaredDevsForEveryModel = rbind(SumOfSquaredDevs_C1, SumOfSquaredDevs_C2, SumOfSquaredDevs_C3, SumOfSquaredDevs_C4, SumOfSquaredDevs_C5, SumOfSquaredDevs_C6)
  rownames(SumOfSquaredDevsForEveryModel) = c("AllIncludedWithProb0.5", "AllIncludedWithProb0.1", "Independent0.5", 
                                              "Independent0.1", "FullModel", "NullModel")
  
  PLsForEveryModel = rbind(PLs_C1, PLs_C2, PLs_C3, PLs_C4, PLs_C5, PLs_C6)
  rownames(PLsForEveryModel) = c("AllIncludedWithProb0.5", "AllIncludedWithProb0.1", "Independent0.5", 
                                 "Independent0.1", "FullModel", "NullModel")
  
  # Return the results in a list
  list(SumOfSquaredDevsForEveryModel = SumOfSquaredDevsForEveryModel, PLsForEveryModel = PLsForEveryModel)
  
}
stopCluster(cl)
end <- Sys.time()
end - start

################################################################################
### Combining the results from each replication  ###############################
################################################################################

save(FiveModelsFiveDatasetsThirtyReps, file = paste(currentwd, "ModelSimulations30RepsResults_V2.rda", sep = ""))
save(FiveModelsFiveDatasetsTenReps, file = paste(currentwd, "ModelSimulations20RepsResults_V2.rda", sep = ""))

# Loading in the saved results for both 20 replicates and 30 replications. 
# Going to combine them together to get 50 total replications. 
load(paste(currentwd, "ModelSimulations30RepsResults_V2.rda", sep = ""))
load(paste(currentwd, "ModelSimulations20RepsResults_V2.rda", sep = ""))

# Combining the two sets of replications together
FiveModels50Reps <- c(FiveModelsFiveDatasetsTenReps,
                      FiveModelsFiveDatasetsThirtyReps)

# Checking the length looks correct (should be 50). 
length(FiveModels50Reps)

# Selecting just the sums of squared deviations and the posterior likelihoods. 
FiveModelsSumOfSquaredDevs = sapply(FiveModels50Reps, '[', 1)
FiveModelsPLs = sapply(FiveModels50Reps, '[', 2)

# Taking the mean of each list of matrices above
FiveModelsSumOfSquaredDevsMean = Reduce('+', FiveModelsSumOfSquaredDevs)/length(FiveModelsSumOfSquaredDevs)
FiveModelsPLsMean = Reduce('+', FiveModelsPLs)/length(FiveModelsPLs)

# Doing pairwise t-tests between the results within each condition for each model
# Storing the number of conditions and models
n_conditions = nrow(FiveModelsSumOfSquaredDevsMean)
n_models = ncol(FiveModelsSumOfSquaredDevsMean)

# Recall that the rows are the conditions and the column names are the model types
conditions = rownames(FiveModelsSumOfSquaredDevsMean)
models = colnames(FiveModelsSumOfSquaredDevsMean)

# Initializing a list of matrices to contain p-values
ListOfConditionsMatrixOfPValuesSumOfSquaredDevs = lapply(1:n_conditions, function(cond) list()) 
names(ListOfConditionsMatrixOfPValuesSumOfSquaredDevs) = rownames(FiveModelsSumOfSquaredDevsMean)

ListOfConditionsMatrixOfPValuesPLs = lapply(1:n_conditions, function(cond) list()) 
names(ListOfConditionsMatrixOfPValuesPLs) = rownames(FiveModelsPLsMean)


for (cond in 1:n_conditions) { # for each condition
  # Creating matrices to store the p-values in
  MatrixOfPValuesSSDForCond.i = matrix(nrow = n_models, ncol = n_models)
  rownames(MatrixOfPValuesSSDForCond.i) = colnames(FiveModelsSumOfSquaredDevsMean)
  colnames(MatrixOfPValuesSSDForCond.i) = colnames(FiveModelsSumOfSquaredDevsMean)
  
  MatrixOfPValuesPLsForCond.i = matrix(nrow = n_models, ncol = n_models)
  rownames(MatrixOfPValuesPLsForCond.i) = colnames(FiveModelsPLsMean)
  colnames(MatrixOfPValuesPLsForCond.i) = colnames(FiveModelsPLsMean)
  
  # Iterating through all pairs of models
  for (mod1 in 1:n_models) {
    for (mod2 in 1:n_models) {
      # Doing a t-test to compare the results from each model run for each condition
      # For sum of squared deviations
      # Need to add conditions because for full and null model, their ssds were entirely 1 or 0 
      # which throws an error. 
      if ((cond == 5 & mod1 == 3 & mod2 == 5) | (cond == 5 & mod1 == 5 & mod2 == 3) |
          (cond == 6 & mod1 == 3 & mod2 == 5) | (cond == 6 & mod1 == 5 & mod2 == 3) | 
          (cond == 5 & mod1 == 5 & mod2 == 5) | (cond == 6 & mod1 == 3 & mod2 == 3)) {
        
        MatrixOfPValuesSSDForCond.i[mod1, mod2] = 0
   
      } else {
        t.test.results.ssd = t.test(sapply(FiveModelsSumOfSquaredDevs, function(repl) repl[cond,mod1]), 
                                    sapply(FiveModelsSumOfSquaredDevs, function(repl) repl[cond,mod2]),
                                    paired = TRUE)
        
        # Storing the results in their respective matrices
        # For sum of squared deviations
        MatrixOfPValuesSSDForCond.i[mod1, mod2] = t.test.results.ssd$p.value
      }
      
      # Storing the results in their respective matrices
      # For posterior likelihoods
      t.test.results.pls = t.test(sapply(FiveModelsPLs, function(repl) repl[cond,mod1]), 
                                  sapply(FiveModelsPLs, function(repl) repl[cond,mod2]),
                                  paired = TRUE)
      
      # For posterior likelihood
      MatrixOfPValuesPLsForCond.i[mod1, mod2] = t.test.results.pls$p.value
    }
  }
  
  # Storing the matrix in the list
  # For sum of squared deviations
  ListOfConditionsMatrixOfPValuesSumOfSquaredDevs[[cond]] = MatrixOfPValuesSSDForCond.i
  
  # For posterior likelihood
  ListOfConditionsMatrixOfPValuesPLs[[cond]] = MatrixOfPValuesPLsForCond.i
}

# Creating tables
library(xtable)

# Adding more informative titles
colnames(FiveModelsSumOfSquaredDevsMean) = c("Hierarchical Spike & Slab",
                                                    "Inclusion Probability Fixed at 0.5",
                                                    "Full Model",
                                                    "Shared Inclusion Probability for all Predictors",
                                                    "Null Model")
rownames(FiveModelsSumOfSquaredDevsMean) = c("All Included w.p. 0.5",
                                                    "All Included w.p. 0.1",
                                                    "Indep. Inclusion w.p. 0.5",
                                                    "Indep. Inclusion w.p. 0.1",
                                                    "All Covariates Included",
                                                    "No Covariates Included")
  
xtable(FiveModelsSumOfSquaredDevsMean, digits = 4)
xtable(FiveModelsPLsMean)

# Checking certain p-values to see if they are significant
round(ListOfConditionsMatrixOfPValuesSumOfSquaredDevs$Independent0.5, 3)
