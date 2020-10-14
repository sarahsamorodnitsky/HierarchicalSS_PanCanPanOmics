# Purpose: Build hierarchical model in which not all predictors are shared 
# by every cancer type. Assumes survival is **normally** distributed.
# This model implements a spike-and-slab prior on the effect of each 
# principal component on survival. 
# Author: Sarah Samorodnitsky and Eric Lock (2020)
# University of Minnesota

library(MASS)
library(truncnorm)
library(EnvStats)

HierarchicalNormalSpikeSlab = function(Covariates, Survival, Censored, 
                                      starting_values, iters, covariates_in_model, 
                                       covariates_by_cancer, priors) {
  # covariates_in_model (int vec): A vector of PC identifiers that are included in the dataset
  # covariates_by_cancer (list of int vecs): A list of vectors for each of the cancer types. Each vector has the PC identifiers that are available in the data (not every cancer has the same PCs)
  # priors (list of doubles): A list of doubles that includes the prior values for variance of beta_tilde, shape and rate of lambda2, shape and rate of sigma2
  
  logSum <- function(l){
    #given l=log(v), where v is a vector,
    #computes log(sum(v))
    max(l) + log(sum(exp(l-max(l))))
  }
  
  Insert0sIntoVec = function(vec, avail, add_NA) {
    # Insert 0s or NAs into indices of the vector where no estimate is available.
    # Use for adding 0s/NAs to substitute beta's that don't exist for certain cancers. 
    # vec (int vec): vector in which to add 0s
    # avail (logical vec): indices where estimates are not available (TRUE = available, FALSE = unavailable)
    # add_NA (logical): add_NA = TRUE then add NAs to replace missing betas, add_NA = FALSE then add 0s
    if (add_NA) { # add NAs instead of 0s
      new_vec = rep(NA, length(avail))
      new_vec[avail] = vec
    } else { # add 0s instead of NAs
      new_vec = rep(0, length(avail))
      new_vec[avail] = vec
    }
    return(new_vec)
  }
  
  NumIncludedAndAvail = function(covariates_by_cancer, covariates_in_model, gamma_i) {
    # Returns the number of cancer types for each covariate in the model that have that covariate
    # both available to them and also have included it in their own model.
    # i.e. the number of cancers for which each covariate is available and included in the slab.
    # Return vector is the length of the number of covariates in the model
    number_included_and_available = sapply(covariates_in_model, function(predictor) {
      sum(unlist(sapply(1:n_cancer, function(type){
        covariates_for_cancer_i = covariates_by_cancer[[type]]
        if (predictor %in% covariates_for_cancer_i) {
          gamma_i[[type]][covariates_for_cancer_i == predictor] == 1
        }
      })))
    }) # for each covariate, count how many cancers have a measurement for it
    
    names(number_included_and_available) = covariates_in_model
    return(number_included_and_available)
  }
  
  # For testing
  # x_i = X[[j]]; y_i = Y[[j]]
  # sigma2_i = sigma2[i]
  # lambda2_i = lambda2[i,][avail]
  # beta_tilde_i = beta_tilde[i,][avail]
  
  # Functions to help calculate posteriors
  Beta_Posterior_Helper = function(x_i, y_i, gamma_ij, sigma2_i, lambda2_i, beta_tilde_i) {
    # Returns mean and variance of multivariate normal posterior for each beta_i based on location in spike or slab. 
    # gamma_i is the vector of inclusion indicators for cancer type i's available covariates
    
    x_i = as.matrix(x_i); y_i = as.matrix(y_i)
    
    # Setting up the prior for the betas given lambda2, betatilde
    # Setting up alpha_i vector
    alpha_i = (1/spike_priorvar) * lambda2_i
    alpha_i[gamma_ij == 0] = 1
    
    # Setting up the covariance matrix for the coefficients
    Variance = spike_priorvar * alpha_i
    
    # Mean vector
    beta_tilde_after_adding_0s = beta_tilde_i * gamma_ij
    
    B = solve((1/sigma2_i)*t(x_i)%*%x_i + diag(1/Variance), tol = 1e-35) # if (length(lambda2_i) == 1) solve((1/sigma2_i)*t(x_i)%*%x_i + (1/lambda2_i)) else 
    b = (1/sigma2_i) * t(x_i) %*% y_i + diag(1/Variance) %*%  beta_tilde_after_adding_0s
    mean = B%*%b; var = B
    
    return(list(mean = mean, var = var))
  }
  
  Lambda2_Posterior_Helper = function(current_betas, current_beta_tilde, gamma_i, n_cancer) {
    # Calculates posterior variance for each coefficent across all 4 cancer types.
    # Computing this posterior requires calculating a sum of squared differences between
    # \beta_{ip} and \tilde\beta_p for i=1,..., N where N = number of cancers. 
    # Since not every cancer contains a beta for each covariate, insert 0s where no
    # squared difference is available so that it does not effect the overall sum. 
    # Also don't want to include the estimated \beta_{ip} which are in the spike,
    # so I knock out non-slab terms and insert 0s instead. 
    tot = rep(0, p)
    for (type in 1:n_cancer) {
      ps_i = covariates_by_cancer[[type]]
      avail = covariates_in_model %in% ps_i 
      gamma_ij = gamma_i[[type]]
      summand_i = (current_betas[[type]] - current_beta_tilde[avail])^2 * gamma_ij
      
      # need to insert 0s where necessary so it does not contribute to the vector tot entries 
      summand_i_2 = Insert0sIntoVec(summand_i, avail, add_NA = FALSE)
      
      tot = tot + summand_i_2
    }
    return(tot)
  }
  
  Sigma2_Posterior_Rate_Helper = function(X, Y, current_betas_augmented, n_cancer) {
    # returns the rate for the posterior of sigma2
    tot = 0
    for (i in 1:n_cancer) {
      x_i = X[[i]]; y_i = Y[[i]]; current_betas_i = current_betas_augmented[[i]]
      x_i = as.matrix(x_i); y_i = as.matrix(y_i)
      tot = tot + sum((y_i - x_i %*% current_betas_i)^2)
    }
    return(tot)
  }
  
  Beta_Tilde_Posterior_Helper = function(current_betas_mean, current_lambda2, num_inc_avail, tau2) {
    # Calculates the posterior parameters for \tilde\beta
    post_beta_tilde_mean = (num_inc_avail*tau2*current_betas_mean)/(current_lambda2 + num_inc_avail*tau2)
    post_beta_tilde_var = (current_lambda2*tau2)/(current_lambda2 + num_inc_avail*tau2)
    return(list(mean = post_beta_tilde_mean, var = post_beta_tilde_var)) 
  }
  
  # Creating new variables to store the data so they can be amended without altering 
  # the originals during testing. 
  X = Covariates
  Y = Survival
  
  p = length(covariates_in_model) # the number of covariates to estimate (including the intercept)
  n_vec =  c(unlist(lapply(X, nrow))) # the number of observations per group
  n_cancer = length(n_vec) # the number of cancer types
  I_p = sapply(covariates_in_model, function(i) sum(sapply(covariates_by_cancer, function(cancer) i %in% cancer))) # for each covariate, count how many cancers have a measurement for it
  
  # Prior values
  # beta tilde
  tau2_intercept = priors$betatilde_priorvar_intercept
  tau2_coefficient = priors$betatilde_priorvar_coefficient
  # lambda2
  lambda2_priorshape_intercept = priors$lambda2_priorshape_intercept
  lambda2_priorrate_intercept = priors$lambda2_priorrate_intercept
  lambda2_priorshape_coefficient = priors$lambda2_priorshape_coefficient
  lambda2_priorrate_coefficient = priors$lambda2_priorrate_coefficient
  # sigma2
  sigma2_priorshape = priors$sigma2_priorshape
  sigma2_priorrate = priors$sigma2_priorrate
  # spike for beta
  spike_priorvar = priors$spike_priorvar
  
  # Starting values
  sigma2_start = starting_values$sigma2_start
  beta_tilde_start = starting_values$beta_tilde_start
  lambda2_start = starting_values$lambda2_start
  gamma_start = starting_values$gamma_start
  pi_start = starting_values$pi_start
  
  # Other important quantities
  total_obs = sum(n_vec) 
  
  # Censored observations
  num_censored = sapply(lapply(Y, is.na), sum) # the number of censored observations in each cancer
  which_censored = sapply(Y, is.na) # TRUE/FALSE list
  
  # initializing the objects to store the posterior samples
  betas = lapply(1:n_cancer, function(i) list()) # list to store posterior betas
  beta_tilde = matrix(ncol = p, nrow = iters + 1) # So beta_tilde is a 
  beta_tilde[1,] = beta_tilde_start
  lambda2 = matrix(ncol = p, nrow = iters + 1)
  lambda2[1,] = lambda2_start
  sigma2 = c()
  sigma2[1] = sigma2_start
  gamma = gamma_start # use the initial value object to store the posterior samples
  pi = matrix(ncol = p, nrow = iters + 1) # to store the posterior samples for the inclusion probability of each covariate
  pi[1,] = pi_start
  
  # Initialize censored observations to have their last_contact time be their survival time. 
  # Within loop, generate their survival times randomly. 
  for (k in 1:n_cancer) {
    Y[[k]][which_censored[[k]]] = Censored[[k]][which_censored[[k]]]
  }
  
  for (i in 1:iters) {
    # svMisc::progress(i/(iters/100))
    
    # First, calculating the number of cancers that have each covariate available
    # to them and are in the slab. 
    gamma_i = sapply(gamma, `[`, i) # the most up-to-date inclusion indicators
    number_included_and_available = NumIncludedAndAvail(covariates_by_cancer, covariates_in_model, gamma_i) # number of cancers for each covariate that have that covariate available and in the slab
    # print(number_included_and_available)
    
    # Posterior for the betas 
    # Note that not every beta vector will be the same length because not every cancer type 
    # has the same number of covariates. So this list of beta posteriors only includes 
    # the available betas for each cancer type. 
    for (j in 1:n_cancer) { 
      ps_j = covariates_by_cancer[[j]]
      avail = covariates_in_model %in% ps_j # logical vector: the covariates that are available for the jth cancer
      lambda2_i = lambda2[i,][avail] # the current estimates for lambda2, subsetting only the available ones
      beta_tilde_i = beta_tilde[i,][avail] # the current estimates for beta_tilde, subsetting only the available ones
      gamma_ij = gamma[[j]][[i]]
      beta_post_params = Beta_Posterior_Helper(X[[j]], Y[[j]], gamma_ij, sigma2[i], lambda2_i, beta_tilde_i)
      beta_j_gen = mvrnorm(1, mu = beta_post_params$mean, Sigma = beta_post_params$var) # posterior beta vector for jth cancer type, includes unnecessary betas
      
      betas[[j]][[i]] = beta_j_gen # for the jth cancer type, the ith iteration of the betas
    }
    
    # posterior for lambda^2
    current_betas = sapply(betas, `[`, i) # current betas for the current iteration
    current_beta_tilde = beta_tilde[i,] # current beta_tildes for the ith iteration
    W = Lambda2_Posterior_Helper(current_betas, current_beta_tilde, gamma_i, n_cancer)
    post_lambda2_shape = (number_included_and_available/2) + c(lambda2_priorshape_intercept, rep(lambda2_priorshape_coefficient, p-1)) # change 1 to 0.01
    post_lambda2_rate = c(lambda2_priorrate_intercept, rep(lambda2_priorrate_coefficient, p-1)) + (0.5)*W # change 1 to 0.01
    lambda2[i+1, ] = 1/rgamma(p, shape = post_lambda2_shape, rate = post_lambda2_rate) 
    
    # posterior for beta tilde
    current_betas_augmented =  lapply(1:n_cancer, function(k) { # add NAs so that missing betas don't contribute to the mean
      betas_i = current_betas[[k]]
      gamma_ij = gamma_i[[k]]
      betas_i[gamma_ij == 0] = NA
      ps_k = covariates_by_cancer[[k]]
      avail = c(covariates_in_model %in% ps_k)
      Insert0sIntoVec(betas_i, avail, add_NA = TRUE) 
    }) # add 0s to indices of betas that do not exist; now all beta vectors will be the same length; necessary for mean computation
    current_betas_mean = colMeans(do.call(rbind, current_betas_augmented), na.rm = TRUE)
    current_betas_mean[is.nan(current_betas_mean)] = 0 # if any of nan it means that all betas were in the spike this iteration
    current_lambda2 = lambda2[i+1,] 
    post_beta_tilde = Beta_Tilde_Posterior_Helper(current_betas_mean, current_lambda2, number_included_and_available, tau2 = c(tau2_intercept, rep(tau2_coefficient, p-1)))
    beta_tilde[i+1, ] = rnorm(p, post_beta_tilde$mean, sqrt(post_beta_tilde$var))

    # posterior for sigma^2
    post_sigma2_shape = (total_obs/2) + sigma2_priorshape # change 1 to 0.01
    post_sigma2_rate = 0.5*Sigma2_Posterior_Rate_Helper(X, Y, current_betas, n_cancer) + sigma2_priorrate # change 1 to 0.01
    sigma2[i+1] = 1/rgamma(1, post_sigma2_shape, rate = post_sigma2_rate)
    
    # posterior for pi (the inclusion probabilities)
    number_of_inclusions = NumIncludedAndAvail(covariates_by_cancer, covariates_in_model, gamma_i)
    pi[i+1,] = rbeta(p, 1 + number_of_inclusions, 1 + I_p - number_of_inclusions)
    
    # posterior for gamma (the vector of inclusion indicators)
    for (j in 1:n_cancer) {
      ps_j = covariates_by_cancer[[j]]
      avail = covariates_in_model %in% ps_j
      pi_j = pi[i+1, ][avail]
      beta_tilde_i = beta_tilde[i+1, ][avail]
      lambda2_i = lambda2[i+1,][avail]
      
      like_slab = dnorm(current_betas[[j]], mean = beta_tilde_i, sd = sqrt(lambda2_i), log = TRUE)
      like_spike = dnorm(current_betas[[j]], mean = 0, sd = sqrt(spike_priorvar), log = TRUE)
      
      probs_j = sapply(1:length(ps_j), function(k) {
        x = log(pi_j[k]) + like_slab[k]
        y = log(1 - pi_j[k]) + like_spike[k]
        exp(x - logSum(c(x,y)))
      })
      
      # Check probs_j equivalent to (pi_j*like_slab)/(pi_j*like_slab + (1-pi_j)*like_spike) for non-log like_slab and like_spike
      new_gammas = rbinom(length(ps_j), 1, probs_j)
      new_gammas[1] = 1
      gamma[[j]][[i+1]] = new_gammas
    }
    
    # Generate Ys for censored observations
    for (k in 1:n_cancer) { # n needs to be the number of cancer types
      n_gens = num_censored[[k]]
      Censored_Obs = X[[k]][which_censored[[k]],]
      Censor_Lower_Bound = Censored[[k]][which_censored[[k]]]
      
      if (n_gens != 0) { # only do this if there are censored observations
        Mu_Survival = if(length(current_betas[[k]]) == 1) Censored_Obs * current_betas[[k]] else Censored_Obs %*% current_betas[[k]] # two cases: a cancer may only have the intercept and no PCs available
        random_survival = rtruncnorm(n_gens, a = Censor_Lower_Bound,
                                     mean = Mu_Survival, sd = rep(sqrt(sigma2[i]), n_gens))
        Y[[k]][which_censored[[k]]] = random_survival
      }
    }
  }
  
  return(list(betas = betas, beta_tilde = beta_tilde, lambda2 = lambda2, sigma2 = sigma2,
              gamma = gamma, pi = pi))
}


