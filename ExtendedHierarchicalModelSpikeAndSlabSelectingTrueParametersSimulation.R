### Purpose: Testing the selection accuracy of spike-and-slab self-written model. This means 
### Checking how well the self-written sampler correctly identifies which covariates
### should be in the model for each cancer type and which should not be. 
# Author: Sarah Samorodnitsky and Eric Lock (2020)
# University of Minnesota

library(MASS)

source("ExtendedHierarchicalModelNormalSpikeAndSlab.R")

SimulateSelectionAccuracy = function(nsim, giter) {
  # Simulate the accuracy of my self-written model by selecting variables that are relevant 
  # to the response. Using a log-normal survival model here.
  # (have non-zero \beta coefficient). 
  # nsim (int): the number of iterations of the simulation to run
  # giter (int): the number of iterations to burn-in for JAGS and the number of additional posterior samples to draw
  
  DataBasedOnTrueValuesAndInclusion = function(true_betas, true_sigma2, n_vec, covariates_by_cancer) {
    # Generate data in the same format as the TCGA data so as to use the function written above
    # Survival outcomes include 50% censored observations.
    # true_betas (list of vectors of floats): beta values for each covariate, will have unnecessary betas
    # true_sigma2 (float): true variance of survival times across cancer types
    # n_vec (int vector): sample size of each cancer type, should have same length of number of cancers
    
    # Initialize lists to store generated data
    Design = list()
    Response = list()
    Censored = list()
    
    I = length(true_betas) # number of cancer types
    
    for (i in 1:I) { # for each cancer type
      ps_i = covariates_by_cancer[[i]] # the covariates available for cancer type i
      n_i = n_vec[i] # sample size for the ith cancer type
      beta_coefs = as.matrix(true_betas[[i]]) # extracting just the betas that exist for the ith cancer type
      num_covariates = length(ps_i)
      
      X = matrix(nrow = n_i, ncol = num_covariates)
      colnames(X) = ps_i
      
      # Fill in the feature matrix
      for (j in 1:num_covariates) { # iterate through the covariates
        X[,j] = rnorm(n_i, 0, 1)
      }
      
      # Multiply the vector of beta coefficients by the inclusion indicator, will set betas equal to 0 if they are not relevant to the response
      # Thus, the dimension of beta_coefs does not change after multiplication
      Y = matrix(rnorm(n_i, mean = X %*% beta_coefs, sd = sqrt(true_sigma2)), nrow = n_i, ncol = 1) # generating survival times
      C = matrix(rnorm(n_i, mean = X %*% beta_coefs, sd = sqrt(true_sigma2)), nrow = n_i, ncol = 1) # generating censor times
      
      # Decide which observations to censor: 
      # if censor time comes before survival time, then subject is censored and we don't know how long they survived until
      num_censored = sum(Y > C)
      num_not_censored = length(Y) - num_censored
      which_to_censor = Y > C # the observations which are censored
      which_to_not_censor = Y <= C # the observations which are not censored
      
      Y[which_to_censor,] = rep(NA, num_censored) # censor all observations whose survival times exceeded last contact times
      C[which_to_not_censor,] = rep(NA, num_not_censored) # ignore all censored times for observations whose survival times were less than their censor times
      
      
      Design[[i]] = X
      Response[[i]] = Y
      Censored[[i]] = C
    }
    return(list(Design = Design, Response = Response, Censored = Censored))
  }
  
  # Characteristics of the data I selected for the simulation
  # n_betas = n_cancer*4 # total number of parameters, for iterating through intervals (4 coefficients for each cancer)
  # covariates_by_cancer = sapply(1:n_cancer, function(i) list(1:4))
  covariates_by_cancer = list(c(1, 2, 4), c(1, 2, 3, 4), c(1, 3), c(1, 2, 3), c(1, 4), c(1, 2, 4),
                              c(1, 3), c(1, 2), c(1, 2, 3), c(1, 2, 3, 4), c(1, 2), c(1, 3)) 
  n_cancer = length(covariates_by_cancer) # number of cancer types
  n_vec = sample(200:500, size = n_cancer) # sample size for each of 3 cancer types, randomly generated.
  
  # Useful values
  iters = giter
  covariates_in_model = 1:4
  n_betas = length(unlist(covariates_by_cancer))
  p = length(covariates_in_model)

  SelPerc = rep(0, n_betas) # initialize vector to store the selection percentage of each beta
  Differences = lapply(1:nsim, function(i) list()) # to store any differences that occur between the true gammas and the inclusion probabilities
  
  for (i in 1:nsim) {
    svMisc::progress(i/(nsim/100)) 
    
    priors = list(betatilde_priorvar_intercept = 10^2,    # \tilde\beta_0 ~ N(0,100). both beta tildes are centered at 0
                  betatilde_priorvar_coefficient = 1,     # -> \tilde\beta_p ~ N(0,1), p > 0
                  lambda2_priorshape_intercept = 1,       # -> \lambda^2_0 ~ IG(1,1)
                  lambda2_priorrate_intercept = 1,         
                  lambda2_priorshape_coefficient = 5,     # -> \lambda^2_p ~ IG(5,1) p > 0
                  lambda2_priorrate_coefficient = 1,
                  sigma2_priorshape = 1,                # -> \sigma^2 ~ IG(1, 1)
                  sigma2_priorrate = 1,
                  spike_priorvar = 1/10000)
    
    betatilde_priorvar = c(priors$betatilde_priorvar_intercept, rep(priors$betatilde_priorvar_coefficient, p-1))
    true_beta_tilde = rnorm(p, 0, sqrt(betatilde_priorvar)) # generate from the prior
    
    true_sigma2 = 1/rgamma(1, priors$sigma2_priorshape, priors$sigma2_priorrate) # generate from the prior
    
    lambda2_priorshape = c(priors$lambda2_priorshape_intercept, rep(priors$lambda2_priorshape_coefficient, p-1))
    lambda2_priorrate = c(priors$lambda2_priorrate_intercept, rep(priors$lambda2_priorrate_coefficient, p-1))
    true_lambda2 = 1/rgamma(p, lambda2_priorshape, lambda2_priorrate) # generate from the prior
    
    true_pis = rbeta(p, 1, 1)
    true_pis[1] = 1
    true_gammas = lapply(1:n_cancer, function(type) {
      gammas = rbinom(length(covariates_by_cancer[[type]]), 1, true_pis)
      gammas
    })
    names(true_gammas) = 1:n_cancer
    true_betas = lapply(1:n_cancer, function(type) {
      true_gammas_i = true_gammas[[type]]
      ps_i = covariates_by_cancer[[type]]
      avail = covariates_in_model %in% ps_i
      var_vec = true_lambda2[avail]
      var_vec[true_gammas_i == 0] = 1/10000
      mvrnorm(1, (true_beta_tilde[avail]) * true_gammas_i, diag(var_vec))
    })
    
    # Generate the data based on the generated true values. 
    SimulatedData = DataBasedOnTrueValuesAndInclusion(true_betas, true_sigma2, n_vec, covariates_by_cancer)
    Covariates = SimulatedData$Design # The covariate values for each cancer type, each cancer may have a varying number of columns. The columns have a name indicating which predictor they correspond to. 
    Survival = SimulatedData$Response # the survival time for each observation
    Censored = SimulatedData$Censored # the censor times for each observation
    
    # Starting values
    beta_tilde_start = rnorm(p, 0, 1)
    sigma2_start = 1/rgamma(1,1,1)
    lambda2_start = 1/rgamma(p,1,1)
    pi_start = rbeta(p,1,1)
    
    # Initialize inclusion indicators for every cancer type as in the slab
    gamma_start = lapply(1:n_cancer, function(type) {
      indicators = rbinom(length(covariates_by_cancer[[type]]), size = 1, prob = pi_start[covariates_by_cancer[[type]]])
      names(indicators) = covariates_by_cancer[[type]]
      list(indicators)
    })
    
    starting_values = list(sigma2_start = sigma2_start,
                           beta_tilde_start = beta_tilde_start,
                           lambda2_start = lambda2_start,
                           gamma_start = gamma_start, # this will get passed on to be the list that stores posterior draws
                           pi_start = pi_start)
    
    sim.out = HierarchicalNormalSpikeSlab(Covariates, Survival, Censored, 
                                          starting_values, iters, 
                                          covariates_in_model, 
                                          covariates_by_cancer, priors)
    
    # Preparing the results in the form of an inclusion indicator vector
    inc_ind = lapply(1:n_cancer, function(cancer) {
      # Taking the mean of the vector of gammas - gives proportion included and not included in each iteration
      mat_of_gammas = do.call(rbind, sim.out$gamma[[cancer]])[(iters/2 + 1):iters, ]
      inclusion_prob = colMeans(mat_of_gammas)
      as.numeric(inclusion_prob >= 0.5)
    })
    
    # setting up the true gammas vector and the inclusion vector from JAGS for comparison; adding names
    # cancers_4 = rep(names(true_gammas), each = 4) # to name each entry in the inclusion prob vector of which cancer type the betas came from
    # cancers_and_betas = paste(cancers_4, rep(1:4)) # combining the cancers with their respective betas
    cancers = sapply(1:n_cancer, function(cancer) {
      rep(cancer, length(true_gammas[[cancer]]))
    })
    cancers_and_betas = unlist(sapply(1:n_cancer, function(type) {
      paste(cancers[[type]], covariates_by_cancer[[type]])
    }))
    
    inc_ind_vec = unlist(inc_ind)
    names(inc_ind_vec) = cancers_and_betas
    
    true_gammas_vec = unlist(unname(true_gammas))
    names(true_gammas_vec) = cancers_and_betas
    
    # add the number of correctly identified inclusion probabilities to the total
    SelPerc = SelPerc + as.numeric(inc_ind_vec == true_gammas_vec)
  
    # if the inclusion indicators from the model differed in any way from the true inclusion indicators
    if (!all(inc_ind_vec == true_gammas_vec)) {
      different = which(inc_ind_vec != true_gammas_vec)
      
      # creating data frame to store where the differences were
      dt_of_diffs = matrix(nrow = length(different), ncol = 2) 
      colnames(dt_of_diffs) = c("cancer", "beta")
      
      # filling in the differences matrix
      for (k in 1:length(different)) {
        current_diff = names(different[k])
        current_diff_remodel = as.numeric(unlist(strsplit(current_diff, " ")))
        dt_of_diffs[k, ] = current_diff_remodel
      }
      Differences[[i]] = dt_of_diffs
    }
  }
  
  # To return a percentage
  SelPerc = SelPerc/nsim
  return(list(SelPerc = SelPerc, Differences = Differences))
}

testsel = SimulateSelectionAccuracy(nsim = 1000, giter = 5000)
