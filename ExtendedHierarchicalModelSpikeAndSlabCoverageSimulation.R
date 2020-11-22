# Purpose: Check coverage of extended hierarchical model WITH spike-and-slab priors (in which not every covariate exists for every cancer type)
#   by generating data using "true parameter values" and checking if the model-derived credible intervals 
#   cover the true values ~ 95% of the time. 
# This is a coverage simulation for the spike-and-slab model. 
# Author: Sarah Samorodnitsky and Eric Lock (2020)
# University of Minnesota

source("ExtendedHierarchicalModelNormalSpikeAndSlab.R")

SimulateCoverage = function(nsim, giter) {
  # Simulate the coverage of the Gibbs sampler model. 
  # nsim (int): number of simulations to run (1000 should be good)
  # giter (int): number of iterations for the Gibbs sampler to go through (> 2000 to ensure convergence)
  
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
  
  PosteriorBetasReformat = function(betas, p_gen) {
    # Reformats the posteriors generated for the betas in order to more easily 
    # calculate the credible intervals. 
    # betas (list of lists of ints): list of length = number of cancers, each inner list is a vector with the posteriors for the available betas
    I = length(betas) # number of cancer types
    Reformatted_Betas = list()
    for (i in 1:I) {
      p_i = length(betas[[i]][[1]])
      betas_current = apply(matrix(unlist(betas[i]), 
                                   ncol = p_i, byrow = TRUE)[(iters/2 + 1):iters,], 2, sort)
      colnames(betas_current) = as.character(p_gen[[i]])
      Reformatted_Betas[[i]] = betas_current
    }
    return(Reformatted_Betas)
  }
  
  CredibleInterval = function(parameter) {
    # Returns credible interval for sorted parameter
    # if sorted_parameter is a matrix, then it returns a list of credible intervals (for the betas)
    # otherwise returns a single credible interval, for sigma^2 or lambda^2 or beta_tilde
    interval = list()
    if (is.vector(parameter)) { # if credible interval for sigma^2
      interval = c(quantile(parameter, 0.025), 
                   quantile(parameter, 0.975))
    } 
    else { # if intervals for lambda^2 or \tilde\beta
      p = ncol(parameter)
      for (i in 1:p) {
        current_interval = c(quantile(parameter[,i], 0.025), 
                             quantile(parameter[,i], 0.975))
        interval[[i]] = current_interval
      }
    }
    return(interval)
  }
  
  Check_if_in = function(true_param, interval) {
    # Checks if the true parameter is contained in the credible interval
    # true_param (double)
    # interval (double vector)
    check = true_param >= interval[1] & true_param <= interval[2]
    return(as.numeric(check))
  }
  
  # The available covariates for each cancer (I chose these myself). Note that
  # each cancer has an intercept. 
  covariates_by_cancer = list(c(1, 2, 4), c(1, 2, 3, 4), c(1, 3), c(1, 2, 3), c(1, 4), c(1, 2, 4),
               c(1, 3), c(1, 2), c(1, 2, 3), c(1, 2, 3, 4), c(1, 2), c(1, 3)) 
  
  covariates_in_model = 1:4 # total number of predictors - 3 covariates and the intercept (self-chosen)
  n_cancer = length(covariates_by_cancer) # number of cancer types
  n_betas = length(covariates_in_model) 
  p = n_betas
  n_vec = sample(50:500, size = n_cancer) # sample size for each of 3 cancer types, randomly generated.
  n_params = length(unlist(covariates_by_cancer)) + 3*length(covariates_in_model) + 1 # total number of parameters, for iterating through intervals (4 beta_tildes, 4 lambdas, 1 sigma^2, and all the betas)
  Coverage = rep(0, n_params) # initialize vector to store coverage proportion for each parameter
  
  # Each iteration of the simulation: Generate data, run model, calculate credible intervals, check if model param is contained in interval
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
    true_beta_tilde = rnorm(4, 0, sqrt(betatilde_priorvar)) # generate from the prior
    
    true_sigma2 = 1/rgamma(1, priors$sigma2_priorshape, priors$sigma2_priorrate) # generate from the prior
    
    lambda2_priorshape = c(priors$lambda2_priorshape_intercept, rep(priors$lambda2_priorshape_coefficient, p-1))
    lambda2_priorrate = c(priors$lambda2_priorrate_intercept, rep(priors$lambda2_priorrate_coefficient, p-1))
    true_lambda2 = 1/rgamma(4, lambda2_priorshape, lambda2_priorrate) # generate from the prior

    true_pis = rbeta(n_betas, 1, 1)
    true_pis[1] = 1
    true_gammas = lapply(1:n_cancer, function(type) {
      gammas = rbinom(4, 1, true_pis)
      gammas
    })
    true_betas = lapply(1:n_cancer, function(type) {
      true_gammas_i = true_gammas[[type]]
      ps_i = covariates_by_cancer[[type]]
      avail = covariates_in_model %in% ps_i
      var_vec = true_lambda2
      var_vec[true_gammas_i == 0] = 1/10000
      mvrnorm(1, true_beta_tilde * true_gammas_i, diag(var_vec))[avail]
    })
    
    # Other useful values
    iters = giter
    
    # Generate the data based on the generated true values. 
    SimulatedData = DataBasedOnTrueValuesAndInclusion(true_betas, true_sigma2, n_vec, covariates_by_cancer)
    Covariates = SimulatedData$Design # The covariate values for each cancer type, each cancer may have a varying number of columns. The columns have a name indicating which predictor they correspond to. 
    Survival = SimulatedData$Response # the survival time for each observation
    Censored = SimulatedData$Censored # the censor times for each observation
    
    # Starting values
    beta_tilde_start = rnorm(n_betas, 0, 1)
    sigma2_start = 1/rgamma(1,1,1)
    lambda2_start = 1/rgamma(n_betas,5,1)
    pi_start = rbeta(n_betas,1,1)
    
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
    
    # Running the model, generating the posteriors
    Posteriors = HierarchicalNormalSpikeSlab(Covariates, Survival, Censored, 
                                             starting_values, iters, covariates_in_model, 
                                             covariates_by_cancer, priors)
    
    # Saving the posteriors for each parameter
    beta_tilde = Posteriors$beta_tilde[(iters/2 + 1):iters,]
    sigma2 = Posteriors$sigma2[(iters/2 + 1):iters]
    lambda2 = Posteriors$lambda2[(iters/2 + 1):iters,]
    betas = lapply(1:n_cancer, function(type) {
      mat = do.call(rbind, Posteriors$betas[[type]])[(iters/2):iters,]
      colnames(mat) = covariates_by_cancer[[type]]
      mat
    }) # PosteriorBetasReformat(Posteriors$betas, covariates_by_cancer) # note that the column names are such that 1 = intercept
    pi = Posteriors$pi[(iters/2 + 1):iters,]

    # Creating credible intervals
    beta_tilde_interval = CredibleInterval(beta_tilde)
    sigma2_interval = CredibleInterval(sigma2)
    lambda2_interval = CredibleInterval(lambda2)
    pi_interval = CredibleInterval(pi)
    
    # Creating credible intervals for betas
    beta_intervals = lapply(1:n_cancer, function(type) list())
    for (j in 1:n_cancer) {
      betas_cancer_j = betas[[j]]
      for (k in 1:ncol(betas_cancer_j)) {
        beta_intervals[[j]] = CredibleInterval(betas_cancer_j)
      }
    }
    
    # Removing unnecessary betas
    # true_betas_avail = list()
    # for (j in 1:n_cancer) {
    #   ps_j = covariates_by_cancer[[j]]
    #   ps_avail_j = c(covariates_in_model %in% ps_j) # vector of Ts and Fs for the covariates that are available, include TRUE in front for the intercept
    #   true_betas_avail[[j]] = true_betas[[j]][ps_avail_j]
    # }
    true_betas_avail = true_betas
    
    # Combining everything together
    Intervals = matrix(unlist(list(beta_tilde_interval, sigma2_interval, 
                                   lambda2_interval, beta_intervals, pi_interval)), ncol = 2, byrow = T)
    
    TrueValues = unlist(list(true_beta_tilde, true_sigma2,
                             true_lambda2, unlist(unname(true_betas_avail)), true_pis))
    
    for (j in 1:n_params) {
      if (j >= 1 & j <= length(beta_tilde_interval)) { # check \tilde\beta
        In = Check_if_in(TrueValues[j], Intervals[j,])
        Coverage[j] = Coverage[j] + In
      }
      else if (j >= length(beta_tilde_interval) + 1 & j < length(beta_tilde_interval) + length(sigma2_interval)) {
        In = Check_if_in(TrueValues[j], Intervals[j,])
        Coverage[j] = Coverage[j] + In
      }
      else if (j >= length(beta_tilde_interval) + length(sigma2_interval) & j < 2*length(beta_tilde_interval) + length(sigma2_interval)) { # lambda^2
        In = Check_if_in(TrueValues[j], Intervals[j,])
        Coverage[j] = Coverage[j] + In
      }
      else if (j >= 2*length(beta_tilde_interval) + length(sigma2_interval) & j < 2*length(beta_tilde_interval) + length(sigma2_interval) + length(beta_intervals)) {
        In = Check_if_in(TrueValues[j], Intervals[j,])
        Coverage[j] = Coverage[j] + In
      }
      else {
        In = Check_if_in(TrueValues[j], Intervals[j,])
        Coverage[j] = Coverage[j] + In
      }
    }
  }
  Coverage = Coverage/nsim
  return(Coverage)
}

Test = SimulateCoverage(nsim = 1000, giter = 2000)


