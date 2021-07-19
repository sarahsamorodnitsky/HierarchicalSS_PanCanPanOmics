### Script with helper functions to import each time they are used elsewhere 
### to eliminate the problem of updating a function in one script but not in
### another.
### Author: Sarah Samorodnitsky
### Date: 12/9/19

library(Matrix)

load("~/PanTCGA/SpikeAndSlabResults/SpikeAndSlabGibbsSampler/ModelWithInterceptOutsideSS/ResultsAfterFiltering_V3/XYC_V2_WithAge_StandardizedPredictors.rda") # loads the new factorization data with all PCS and age standardized

# Useful variables
cancer_types = names(Covariates) # should be 29 cancer types
covariates_by_cancer_full = sapply(Covariates, function(type) as.numeric(colnames(type)))
covariates_table = table(unlist(covariates_by_cancer_full)) # gives the number of cancer types who have data on each PC
n_vec = sapply(Covariates, nrow)
n_cancer = length(cancer_types)

# Useful functions

# General functions
InitManyMatrices = function(n, n_vec) {
  # initializing a list of matrices where the initialized matrix
  # contains as many rows as there were subjects used in the factorization
  # n (int): number of inner lists to include
  new_list = list()
  for (i in 1:n) {
    new_list[[i]] = matrix(nrow = n_vec[i], ncol = 0)
  }
  return(new_list)
}

InitManyLists = function(n) {
  # initializing a list of lists
  # n (int): number of inner lists to include
  new_list = list()
  for (i in 1:n) {
    new_list[[i]] = list()
  }
  return(new_list)
}

# For forward selection:

Generate5FoldTrainingSet = function(X, Y, cancer_types, n_vec) {
  # Generates 5 folds of training values for each cancer type
  # non-randomly because of the seed
  
  InitManyLists = function(n) {
    new_list = list()
    for (i in 1:n) {
      new_list[[i]] = list()
    }
    return(new_list)
  }
  
  set.seed(1) # setting a seed so the same training data is produced each time
  
  # subsetting Y into training set
  prop = 0.2 # proportion to go into test set
  n = length(cancer_types)
  Training_obs = InitManyLists(n)
  names(Training_obs) = cancer_types
  for (k in 1:n) { # for each cancer type
    current_n = n_vec[k] # size of current cancer type
    reorder_ind = sample(1:current_n, size = current_n, replace = F) # reorder the observations
    current_test_n = round(prop*current_n) # how many observations should go in the test set
    fold_ind = seq(1, current_n, by = current_test_n) # the indices of the start of each test set
    all_ind = 1:current_n
    
    test_folds = lapply(1:length(fold_ind), function(i) { # the indices of observations for the test set
      if (i < length(fold_ind)) {
        test_ind = seq(fold_ind[i], fold_ind[i + 1] - 1) # get a sequence of indices for the test set
        reorder_ind[test_ind] # select from the reordered observations
      } else {
        test_ind = seq(fold_ind[i], current_n) # get a sequence of indices for the test set
        reorder_ind[test_ind] # select from the reordered observations
      }
    })
    
    if (length(fold_ind) > 5) { # if the set of obs was not evenly split into 5 groups
      extra_obs = test_folds[[6]]
      for (i in 1:length(extra_obs)) {
        test_folds[[i %% 5]] = c(test_folds[[i %% 5]], extra_obs[i]) # add an extra obs to another fold in a circular fashion
      }
    }
    test_folds = test_folds[-6] # get rid of the extra set of obs
    
    training_folds = lapply(1:length(test_folds), function(i) { # indices for observations in the training set
      all_ind[!(all_ind %in% test_folds[[i]])]
    })
    
    Training_obs[[k]] = training_folds
  }
  return(list(Training_obs = Training_obs))
}

# To subset out a specific set of PCs
SubsetNPCs = function(X, PCs) {
  # Subsets entire dataset to just a selected set of PCs
  n = length(X)
  PCs = as.character(PCs) # add the intercept to ensure it is not lost
  Xnew = X 
  for (i in 1:n) {
    current_type = X[[i]]
    current_PCs = colnames(current_type)
    current_available = current_PCs %in% PCs
    current_subset = as.matrix(current_type[, current_available])
    colnames(current_subset) = current_PCs[current_available]
    Xnew[[i]] = current_subset
  }
  return(Xnew)
}

# For credible intervals:

PosteriorBetasReformat = function(betas, p_gen) {
  # Reformats the posteriors generated for the betas in order to more easily 
  # betas (list of lists of ints): list of length = number of cancers, each inner list is a vector with the posteriors for the available betas
  I = length(betas) # number of cancer types
  Reformatted_Betas = list()
  for (i in 1:I) {
    # Reorganizes the betas for each cancer type into a matrix 
    p_i = length(betas[[i]][[1]])
    
    if (p_i > 1) { # if a cancer type has more than just the intercept
      betas_current = matrix(unlist(betas[i]), ncol = p_i, byrow = TRUE)[(iters/2 + 1):iters,]
      colnames(betas_current) = as.character(p_gen[[i]])
    } else { # if a cancer type has just the intercept 
      betas_current = matrix(unlist(betas[i]), ncol = p_i, byrow = TRUE)[(iters/2 + 1):iters,,drop = F]
    }
    colnames(betas_current) = as.character(p_gen[[i]])
    Reformatted_Betas[[i]] = betas_current
  }
  return(Reformatted_Betas)
}

CredibleInterval = function(parameter) {
  # Returns credible interval for parameter
  # if parameter is a matrix, then it returns a list of credible intervals (for the betas)
  # otherwise returns a single credible interval, for sigma^2 or lambda^2 or beta_tilde
  interval = list()
  if (is.vector(parameter)) { # if credible interval for sigma^2
    interval = c(quantile(parameter, 0.025), quantile(parameter, 0.975))
  } 
  else { # if intervals for lambda^2 or \tilde\beta
    p = ncol(sorted_parameter)
    for (i in 1:p) {
      current_interval = c(quantile(parameter[,i], 0.025), quantile(parameter[,i], 0.975))
      interval[[i]] = current_interval
    }
  }
  return(interval)
}

PlotCredibleIntervalsForBetas = function(betas, filename) {
  # Given the reformatted raw output from the Gibbs sampler, constructs credible intervals for the betas
  # with a specified burn-in. This output should be after burn-in using the ReformatBetas function. 
  # betas (list of doubles): for i=1,...,29, betas[[i]] = matrix where each column j contains the MCMC results for 
  # the jth PC available for cancer i. 
  # filename (character string): the name of the pdf file to be generated that will contain the 
  # plot images. 
  
  # Necessary functions:
  PlotCovariate = function(IntervalsForThisParam) {
    # Returns a ggplot with the intervals plotted in orange and 0 marked in a green horizontal line. 
    intervals_covariate = matrix(unlist(IntervalsForThisParam), ncol = 2, byrow = T)
    intervals_median = lapply(betas, function(x) {
      present = param %in% colnames(x)
      if (present) {
        median(x[,which(param == colnames(x))])
      }
    })
    intervals_median = intervals_median[cancers_avail]
    intervals_covariate = cbind(unlist(intervals_median), intervals_covariate)
    colnames(intervals_covariate) = c("median", "L", "U")
    
    # adding a column that indicates whether the interval contains 0 or not
    containszero = as.numeric(!(intervals_covariate[,2] > 0 | intervals_covariate[,3] < 0)) # if lower bound of interval is greater than 0 or upper bound is less than 0, set to 0
    containszero[containszero == 0] = -1 # if interval does not contain 0, set containszero to -1
    intervals_covariate = cbind(intervals_covariate, containszero) # 1 if it contains 0, -1 if it does not contain 0
    rownames(intervals_covariate) = cancer_types[cancers_avail]
    
    plot_covariate = ggplot(data.frame(intervals_covariate), aes(x = rownames(intervals_covariate), y = median)) +
      geom_point(size = 2, colour = "dark grey") +
      geom_errorbar(aes(ymax = U, ymin = L, colour = cut(containszero, c(-Inf, 0, Inf)))) +
      scale_color_manual(name = "Interval Contains Zero",
                         values = c("(-Inf,0]" = "#D55E00",
                                    "(0, Inf]" = "black"),
                         labels=c("Does not contain zero", "Contains zero")) + 
      theme(axis.text.x = element_text(angle = 75, vjust = 0.5, size = 11),
            legend.position="none",
            axis.text.y = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            title = element_text(size = 13)) + 
      xlab("Cancer Type") + ylab("Credible Interval") +
      ggtitle(paste("Credible Intervals for Effect of PC", param, "on Survival")) +
      geom_hline(yintercept = 0, color = "#009E73")
    print(plot_covariate)
  }
  
  I = length(betas) # should be the number of cancer types 
  pcs = sort(unique(as.numeric(unlist(sapply(betas, colnames))))) # the pcs used in this model
  
  # Constructing the credible intervals 
  beta_intervals = InitManyLists(I)
  for (i in 1:I) { # for each cancer type
    betas_cancer_i = betas[[i]]
    beta_intervals_i = CredibleInterval(betas_cancer_i)
    names(beta_intervals_i) = colnames(betas_cancer_i)
    beta_intervals[[i]] = beta_intervals_i
  }
  
  pdf(filename) 
  
  # Iterating through all the PCs
  for (param in pcs) {  
    
    # Extracting the credible intervals across all cancer types for just the current param
    IntervalsForThisParam = lapply(beta_intervals, function(type) {
      pcs_type = as.numeric(names(type))
      avail = param %in% pcs_type
      if (avail) {
        type[which(param == pcs_type)]
      }
    })
    
    # Not every cancer type will have a credible interval -> removing those cancers. 
    cancers_avail = !sapply(IntervalsForThisParam, is.null)
    IntervalsForThisParam = IntervalsForThisParam[cancers_avail]
    names(IntervalsForThisParam) = cancer_types[cancers_avail]
    
    PlotCovariate(IntervalsForThisParam)
  }
  dev.off()
}

# For simulations
GenerateTrueValues = function(p, covariates_by_cancer, covariates_in_model, n_cancer, priors, condition, seed) {
  # p (int): the number of covariates in the model
  # n_cancer (int): number of cancers in the model
  # priors (list): contains the hyperparameters for all the prior distributions
  # condition (character): the condition under which to generate the true inclusion indicators
  # condition == (half_present_across_all, 0.1_present_across_all, independent_0.5, independent_0.1, full_all_included, null_none_included)
  # all other parameters generated from their priors
  # seed is the seed to set so the true values are generated in the same way
  
  # Beta tilde
  betatilde_priorvar = c(priors$betatilde_priorvar_intercept, rep(priors$betatilde_priorvar_coefficient, p-1))
  true_beta_tilde = rnorm(p, 0, sqrt(betatilde_priorvar)) # generate from the prior
  
  # Sigma2
  true_sigma2 = 1/rgamma(1, priors$sigma2_priorshape, priors$sigma2_priorrate) # generate from the prior
  
  # Lambda2
  lambda2_priorshape = c(priors$lambda2_priorshape_intercept, rep(priors$lambda2_priorshape_coefficient, p-1))
  lambda2_priorrate = c(priors$lambda2_priorrate_intercept, rep(priors$lambda2_priorrate_coefficient, p-1))
  true_lambda2 = 1/rgamma(p, lambda2_priorshape, lambda2_priorrate) # generate from the prior
  
  # Gammas - true inclusion indicator
  if (condition == "half_present_across_all") { # if predictor k is in slab, then in slab for all cancers. 
    # generate 1s and 0s for each covariate so each will be included with certainty or not across
    # all cancers
    true_pis = rbinom(p, size = 1, prob = 0.5)
    true_pis[1] = 1
    names(true_pis) = covariates_in_model
    
    true_gammas = lapply(1:n_cancer, function(type) {
      # subset the covariates for cancer i
      ps_i = covariates_by_cancer[[type]]
      avail = covariates_in_model %in% ps_i

      # generate indicators using only the pis available for this cancer
      indicators = rbinom(length(covariates_by_cancer[[type]]), size = 1, prob = true_pis[avail])
      indicators[1] = 1
      names(indicators) = covariates_by_cancer[[type]]
      indicators
    })

    # check to make sure that the same inclusion indicator generated for the same covariate across
    # all cancer types
    # all_same = c()
    # i=1
    # for (k in covariates_in_model) {
    #   k=as.character(k)
    #   inclusion_inds_for_k = unlist(sapply(true_gammas, function(type) type[names(type) == k]))
    #   all_same[i] = all((inclusion_inds_for_k == 0)) | all((inclusion_inds_for_k == 1))
    #   i=i+1
    # }
    # all(all_same)
    
  } else if (condition == "0.1_present_across_all") {
    true_pis = rbinom(p, size = 1, prob = 0.1)
    true_pis[1] = 1
    names(true_pis) = covariates_in_model
    
    true_gammas = lapply(1:n_cancer, function(type) {
      # subset the covariates for cancer i
      ps_i = covariates_by_cancer[[type]]
      avail = covariates_in_model %in% ps_i
      
      # generate indicators using only the pis available for this cancer
      indicators = rbinom(length(covariates_by_cancer[[type]]), size = 1, prob = true_pis[avail])
      indicators[1] = 1
      names(indicators) = covariates_by_cancer[[type]]
      indicators 
    })
    
    # check to make sure that the same inclusion indicator generated for the same covariate across
    # all cancer types
    # all_same = c()
    # i=1
    # for (k in covariates_in_model) {
    #   k=as.character(k)
    #   inclusion_inds_for_k = unlist(sapply(true_gammas, function(type) type[names(type) == k]))
    #   all_same[i] = all((inclusion_inds_for_k == 0)) | all((inclusion_inds_for_k == 1))
    #   i=i+1
    # }
    # all(all_same)
    
  } else if (condition == "independent_0.5") { # each predictor with probability 0.5 is present in each cancer type 
    true_gammas = lapply(1:n_cancer, function(type) {
      # subset the covariates for cancer i
      ps_i = covariates_by_cancer[[type]]
      avail = covariates_in_model %in% ps_i
      
      # generate indicators using only the pis available for this cancer
      indicators = rbinom(length(covariates_by_cancer[[type]]), size = 1, prob = 0.5)
      indicators[1] = 1
      names(indicators) = covariates_by_cancer[[type]]
      indicators 
    })
  } else if (condition == "independent_0.1") { # each predictor with probability 0.1 present in each cancer 
    true_gammas = lapply(1:n_cancer, function(type) {
      # subset the covariates for cancer i
      ps_i = covariates_by_cancer[[type]]
      avail = covariates_in_model %in% ps_i
      
      # generate indicators using only the pis available for this cancer
      indicators = rbinom(length(covariates_by_cancer[[type]]), size = 1, prob = 0.1)
      indicators[1] = 1
      names(indicators) = covariates_by_cancer[[type]]
      indicators 
    })
  } else if (condition  == "full_all_included") { # all predictors are included
    true_pis = rep(1, p)
    names(true_pis) = covariates_in_model
    
    true_gammas = lapply(1:n_cancer, function(type) {
      # subset the covariates for cancer i
      ps_i = covariates_by_cancer[[type]]
      avail = covariates_in_model %in% ps_i
      
      # generate indicators using only the pis available for this cancer
      indicators = rbinom(length(covariates_by_cancer[[type]]), size = 1, prob = true_pis[avail])
      indicators[1] = 1
      names(indicators) = covariates_by_cancer[[type]]
      indicators 
    })
  } else if (condition == "null_none_included") { # no predictors are included
    true_pis = rep(0, p)
    true_pis[1] = 1 # intercept always included
    
    true_gammas = lapply(1:n_cancer, function(type) {
      # subset the covariates for cancer i
      ps_i = covariates_by_cancer[[type]]
      avail = covariates_in_model %in% ps_i
      
      # generate indicators using only the pis available for this cancer
      indicators = rbinom(length(covariates_by_cancer[[type]]), size = 1, prob = true_pis[avail])
      indicators[1] = 1
      names(indicators) = covariates_by_cancer[[type]]
      indicators 
    })
  }
  
  # Betas
  true_betas = lapply(1:n_cancer, function(type) {
      true_gammas_i = true_gammas[[type]]
      ps_i = covariates_by_cancer[[type]]
      avail = covariates_in_model %in% ps_i
      var_vec = true_lambda2[avail]
      mvrnorm(1, true_beta_tilde[avail], diag(var_vec)) * true_gammas_i # set to 0 betas that are not included
    })
  
  # Check all the betas that are 0 are the ones that are not included
  # vec=c()
  # for(i in 1:n_cancer){vec[i] = all(which(true_betas[[i]]==0) == which(true_gammas[[i]]==0))}
  # all(vec)
  
  # gathering all the true values together
  true_values = list(true_betas = true_betas, true_gammas = true_gammas, 
                     true_beta_tilde = true_beta_tilde, true_lambda2 = true_lambda2,
                     true_sigma2 = true_sigma2)
  
  return(true_values)
}

DataBasedOnTrueValues = function(true_betas, true_sigma2, n_vec, covariates_by_cancer) {
  # Generate data in the same format as the TCGA data for simulations
  # Survival outcomes include 50% censored observations. 
  # true_betas (list of vectors of floats): beta values for each covariate, will have unnecessary betas
  # true_sigma2 (float): true variance of survival times across cancer types
  # n_vec (int vector): sample size of each cancer type, should have same length of number of cancers
  # covariates_by_cancer (list of int vectors): the available covariates for each cancer type. 0 = intercept, 1 = 1st predictor, 2 = 2nd predictor, ..., k = kth predictor
  
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
    
    X = matrix(nrow = n_i, ncol = num_covariates) # initializing an empty design matrix for cancer i
    colnames(X) = ps_i
    
    for (j in 1:length(ps_i)) { # iterate through the covariates
      if (j == 1) { # adding an intercept
        X[,j] = rep(1, n_i)
      } else {
        X[,j] = rnorm(n_i, 0, 1)
      }
    }
    
    Y = matrix(rlnorm(n_i, mean = X %*% beta_coefs, sd = sqrt(true_sigma2)), nrow = n_i, ncol = 1) # generating survival times
    C = matrix(rlnorm(n_i, mean = X %*% beta_coefs, sd = sqrt(true_sigma2)), nrow = n_i, ncol = 1) # generating censor times
    
    # Decide which observations to censor: 
    # if censor time comes before survival time, then subject is censored and we don't know how long they survived until
    num_censored = sum(Y > C)
    num_not_censored = length(Y) - num_censored
    which_to_censor = Y > C # the observations which are censored
    which_to_not_censor = Y <= C # the observations which are not censored
    
    Y[which_to_censor,] = rep(NA, num_censored) # censor all observations whose survival times exceeded last contact times
    C[which_to_not_censor,] = rep(NA, num_not_censored) # ignore all censored times for observations whose survival times were less than their censor times
    
    # Checking that for every censored observation there is a censor time and vice versa
    # correct_censor=c()
    # for (k in 1:nrow(Y)) {
    #   correct_censor[k] <- (is.na(Y[k,]) & !is.na(C[k,])) | (!is.na(Y[k,]) & is.na(C[k,]))
    # }
    # all(correct_censor)
    
    Design[[i]] = X
    Response[[i]] = Y
    Censored[[i]] = C
  }
  return(list(Design = Design, Response = Response, Censored = Censored))
}

SumOfSquaredDeviations = function(true_gammas, Posteriors, iters, mod) {
  # Computes the sum of squared deviations between the true inclusion indicator
  # and the posterior inclusion probability. Also includes a burn-in so only
  # takes the inclusion indicators after the burn-in. 
  
  if (mod != "separate_model") {
    # Saving the number of cancers used in the model
    n_cancer <- length(Posteriors$gamma) 
    
  } else {
    n_cancer <- length(true_gammas)
    Posteriors_temp <- lapply(1:1, function(type) list()); names(Posteriors_temp) <- "gamma" # creating a list with just one element "gamma" to store the resulting inclusions for each cancer
    Posteriors_temp[[1]] <- lapply(1:n_cancer, function(type) {
      Posteriors[[type]]$gamma[[1]]
    })
    
    # Saving the resulting inclusion indicators for every cancer type back
    # into the Posteriors variable
    Posteriors <- Posteriors_temp
  }
   
  # Computing the posterior inclusion probability
  thinned_iters = seq(1,(iters/2), by = 10) # for thinning 
  gamma.gibbs = lapply(1:n_cancer, function(type) { # combine the inclusion indicators into a list of matrices for each cancer type so we can take column means
    do.call(rbind, Posteriors$gamma[[type]])[(iters/2):iters, , drop = FALSE][thinned_iters,]
  })
  
  # Computing the posterior inclusion probability
  if (!is.null(dim(gamma.gibbs[[1]]))) { # if there is more than 1 covariate in the model
    posterior_inclusion = lapply(gamma.gibbs, function(type) colMeans(type)) 
  } else { # if we fit the null model
    posterior_inclusion_intercept_only = lapply(gamma.gibbs, function(type) mean(type))
    posterior_inclusion = lapply(1:n_cancer, function(type) {
      c(posterior_inclusion_intercept_only[[type]], rep(0, length(true_gammas[[type]])-1))
    })
  }
  
  # For the joint model, duplicate the resulting inclusion indicator vector for
  # all 29 cancers. Select only the covariates available for that cancer type
  if (mod == "joint_model") { 
    posterior_inclusion_joint = lapply(1:length(true_gammas), function(type) {
      avail = covariates_in_model_full %in% covariates_by_cancer[[type]]
      posterior_inclusion_type = posterior_inclusion[[1]][avail]
      names(posterior_inclusion_type) <- names(true_gammas[[type]])
      posterior_inclusion_type
    })
    
    posterior_inclusion = posterior_inclusion_joint
    
    # Checking that the posterior inclusion length matches the true gammas length
    # all(sapply(1:n_cancer, function(type) length(true_gammas[[type]]) == length(posterior_inclusion[[type]])))
  }
  
  # Computing the sum of squared deviations
  sum_of_squared_devs = sum(unlist(sapply(1:length(true_gammas), function(type) {
    (true_gammas[[type]] - posterior_inclusion[[type]])^2
  })))
  
  # dividing by the total number of covariates in the model (minus all the intercepts, of which there are 29)
  n_covariates = length(unlist(posterior_inclusion)) - 29
  sum_of_squared_devs_mean = sum_of_squared_devs/n_covariates
  
  return(sum_of_squared_devs_mean)
}

ComputeSumOfSquaredDeviationsForEveryModel = function(true_gammas, PosteriorsForEveryModel, iters, model_types) {
  # Returns the sum of squared deviations for all 5 models being considered
  # model_types (character): vector with the names of all the models
  
  # Number of models to consider
  n_models = length(model_types)
  
  # Vector to store all squared deviations in
  squared_devs = c()
  
  for (i in 1:n_models) {
    mod = model_types[i]
    
    PosteriorsForModel_i = PosteriorsForEveryModel[[mod]] # for testing: Posteriors = PosteriorsForModel_i
    squared_devs[i] = SumOfSquaredDeviations(true_gammas, PosteriorsForModel_i, iters, mod)
  }
  
  names(squared_devs) = model_types
  
  return(squared_devs)
}

GenerateStartingValues = function(p, covariates_in_model, covariates_by_cancer) {
  # Returns starting values for each of the 5 model types
  # p (int): number of covariates in the model
  
  # Starting values that will be used for every model
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
  
  # For the hierarchical spike/slab model
  starting_values = list(sigma2_start = sigma2_start,
                         beta_tilde_start = beta_tilde_start,
                         lambda2_start = lambda2_start,
                         gamma_start = gamma_start, # this will get passed on to be the list that stores posterior draws
                         pi_start = pi_start)
  
  # For fixed at 0.5 model
  # We can initialize just one shared value - it will be replicated across the posteriors
  # matrix for every covariate. 
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
  # Generate just a single parameter for the beta tildes and lambda2s and pis
  covariates_in_model_null = covariates_in_model[1]
  p_null = length(covariates_in_model_null)
  covariates_by_cancer_null = lapply(covariates_by_cancer, function(type) type[1]) # just storing the intercept for every cancer
  
  # Generating starting values for the null model
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
  
  # For the joint model
  gamma_start_joint <- list(list(rep(1, length(covariates_in_model)))) # since there is only "1" cancer type here, we only need one list for the inclusion indicators. We thus have a list of one list. 
  starting_values_joint <- starting_values
  starting_values_joint$gamma_start <- gamma_start_joint
  
  # For the separate model
  starting_values_separate <- lapply(1:n_cancer, function(i) list())
  names(starting_values_separate) <- cancer_types
  
  for (type in 1:n_cancer) {
    # Adjusting these parameters to reflect just a single cancer type
    covariates_by_cancer_type <- covariates_by_cancer[type]
    covariates_in_model_type <- unlist(covariates_by_cancer_type)
    
    # Selecting just the starting values that correspond to this cancer type
    avail <- covariates_in_model %in% covariates_by_cancer_type[[1]]
    starting_values_type <- list(sigma2_start = starting_values$sigma2_start,
                                 beta_tilde_start = starting_values$beta_tilde_start[avail],
                                 lambda2_start = starting_values$lambda2_start[avail],
                                 gamma_start = list(list(starting_values$gamma_start[[type]][[1]])),
                                 pi_start = starting_values$pi_start[avail])
    
    # Storing the starting values in the list
    starting_values_separate[[type]] <- starting_values_type
  }
  
  
  # Gathering all the starting values together
  starting_values_for_every_model = list(starting_values = starting_values,
                                         starting_values_fixed_at_0.5 = starting_values_fixed_at_0.5,
                                         starting_values_fixed_at_1.0 = starting_values_fixed_at_1.0,
                                         starting_values_shared_across_betas = starting_values_shared_across_betas,
                                         starting_values_null = starting_values_null,
                                         starting_values_joint = starting_values_joint,
                                         starting_values_separate = starting_values_separate)
  
  return(starting_values_for_every_model)
}

PosteriorLikelihoodForSimulations = function(posteriors, X, Y, Censored, iters, burnin, joint = FALSE) {
  # Computes the posterior likelihood where X, Y, and Censored are simulated test data
  # posteriors is the result from the model computation
  # burnin is the number of iterations from the start to throw away
  
  logSum <- function(l) {
    # given l = log(v), where v is a vector,
    # computes log(sum(v))
    return(max(l) + log(sum(exp(l-max(l)))))
  }
  
  LogLikelihood = function(x_ij, y_ij, y_ijc, betas_it, sigma2_it) {
    # calculate the likelihood for just one set of posteriors
    mu_ij = sum(x_ij * betas_it) 
    
    LogNormalLikelihood = function(y, mu, sig2) {
      # Calculates log density of normal 
      # nl = (1/sqrt(2*pi*sig2))*exp(-(1/(2*sig2))*(y - mu)^2)
      lnl = (-0.5)*log(2*pi*sig2) - (1/(2*sig2))*(y - mu)^2
      return(lnl)
    }
    
    if (!is.na(y_ij)) { # y_ij is not censored
      post_ij = LogNormalLikelihood(log(y_ij), mu = mu_ij, sig2 = sigma2_it) # dnorm(log(y_ij), mean = mu_ij1, sd = sqrt(sigma21_it)) # sample from log-normal
      
    } else { # if y_ij is censored, use censoring time
      post_ij = log(1 - pnorm((log(y_ijc) - mu_ij)/sqrt(sigma2_it), mean = 0, sd = 1)) 
    }
    
    return(post_ij)
  }
  
  # Extracting the necessary information
  n_vec = sapply(X, nrow) # number of patients per cancer type
  betas = posteriors$betas
  sigma2 = posteriors$sigma2
  
  # If this is the joint model, combine the test data together:
  if (joint) {
    # Combine the test data into one list for the joint model
    X <- list(do.call(rbind, CombineCancerDataTogether(X, covariates_in_model_full)))
    Y <- list(unlist(Y))
    Censored <- list(unlist(Censored))
  } else {
    # Converting each entry in Y and Censored to be vectors so that the column names stick in the cbind portion following
    Y = lapply(Y, function(type) as.vector(type))
    Censored = lapply(Censored, function(type) as.vector(type))
  }

  # Obtaining test observations and their times of last contact, etc
  # combine test data into one list
  List_XY = lapply(seq(length(X)), function(i) cbind(X[[i]], "Survival" = Y[[i]], "Censored" = Censored[[i]]))
  
  posterior.vec = list()
  ind_to_use = seq(burnin, iters, by = 10)
  
  for (j in ind_to_use) { ##EFL: changed index 'i' to 'j' here
    
    # for each cancer type
    # for each row in that cancer type
    # apply Likelihood with the corresponding cancer type's posterior parameters
    # at the jth iteration
    #EFL: also changed here 'i' to 'j' on left-hand side below
    posterior.vec[[j]] = unlist(lapply(seq(length(X)), 
                                       function(k) sapply(seq(nrow(List_XY[[k]])), # for each cancer type
                                                          function(i) { # for each subject in cancer type k
                                                            y_col = which(colnames(List_XY[[k]]) %in% "Survival")  # ncol(List_XY[[k]]) - 1 # the column in List_XY that contains the survival time
                                                            yc_col = which(colnames(List_XY[[k]]) %in% "Censored")   # ncol(List_XY[[k]]) # the column in List_XY that contains the censor time
                                                            
                                                            LogLikelihood(x_ij = List_XY[[k]][i, -c(y_col, yc_col)], 
                                                                          y_ij = List_XY[[k]][i, y_col],  # for each patient in the test set
                                                                          y_ijc = List_XY[[k]][i, yc_col], # ith observation, censor time
                                                                          betas_it = betas[[k]][[j]],  # kth cancer type, jth iteration posterior
                                                                          sigma2_it = sigma2[j])
                                                          }
                                       )))
    
  }
  
  #EFL: here is revised code to aggregate the log-likelihoods
  #EFL: organize log-likelihoods into M by N matrix, where M is number of test samples and N is number of iterations
  logLike.mat <- matrix(nrow = sum(n_vec), ncol = length(ind_to_use))
  for (j in ind_to_use) { # corresponds to every 10th sample 
    k = which(j == ind_to_use) # need to explicitly state the column number because the values of j do not correspond
    logLike.mat[,k] <- unlist(posterior.vec[[j]])
  }
  #EFL: compute joint log-likelihoods, at each iteration, by summing the logs
  logLike.joint <- colSums(logLike.mat)
  #EFL: Compute log of full posterior predictive likelihood, using logSum
  PL_Survival <- -log(length(ind_to_use)) + logSum(logLike.joint) # -log(length(posterior.vec))+logSum(logLike.joint)
  
  return(PL_Survival)
}

JointLogLikelihoodForSeparateModelForSimulations = function(posteriors, X, Y, Censored, burnin) {
  # For computing the joint log-likelihood for the separated model
  # posteriors = posterior samples for a single cancer type
  # X = Covariates data for a single cancer type
  # Y = Survival times for a single cancer type
  # Censored = list of censored times for each cancer type
  logSum <- function(l) {
    # given l = log(v), where v is a vector,
    # computes log(sum(v))
    return(max(l) + log(sum(exp(l-max(l)))))
  }
  
  LogLikelihood = function(x_ij, y_ij, y_ijc, betas_it, sigma2_it) {
    # calculate the likelihood for just one set of posteriors
    mu_ij = sum(x_ij * betas_it) 
    
    LogNormalLikelihood = function(y, mu, sig2) {
      # Calculates log density of normal 
      # nl = (1/sqrt(2*pi*sig2))*exp(-(1/(2*sig2))*(y - mu)^2)
      lnl = (-0.5)*log(2*pi*sig2) - (1/(2*sig2))*(y - mu)^2
      return(lnl)
    }
    
    if (!is.na(y_ij)) { # y_ij is not censored
      post_ij = LogNormalLikelihood(log(y_ij), mu = mu_ij, sig2 = sigma2_it) # dnorm(log(y_ij), mean = mu_ij1, sd = sqrt(sigma21_it)) # sample from log-normal
      
    } else { # if y_ij is censored, use censoring time
      post_ij = log(1 - pnorm((log(y_ijc) - mu_ij)/sqrt(sigma2_it), mean = 0, sd = 1)) 
    }
    
    return(post_ij)
  }
  
  # Extracting the necessary information
  n_vec = sapply(X, nrow) # number of patients per cancer type
  betas = posteriors$betas
  sigma2 = posteriors$sigma2
  
  # Combine the test data together
  Y <- lapply(Y, as.vector)
  Censored <- lapply(Censored, as.vector)
  List_XY = lapply(seq(length(X)), function(i) cbind(X[[i]], "Survival" = Y[[i]], "Censored" = Censored[[i]]))
  
  posterior.vec = list()
  ind_to_use = seq(burnin, iters, by = 10)
  
  for (j in ind_to_use) { ##EFL: changed index 'i' to 'j' here
    
    # for each cancer type
    # for each row in that cancer type
    # apply Likelihood with the corresponding cancer type's posterior parameters
    # at the jth iteration
    #EFL: also changed here 'i' to 'j' on left-hand side below
    posterior.vec[[j]] = unlist(lapply(seq(length(X)), # change from X to X_test
                                       function(k) sapply(seq(nrow(List_XY[[k]])), # for each cancer type
                                                          function(i) { # for each test subject in cancer type k
                                                            y_col = which(colnames(List_XY[[k]]) %in% "Survival")  # ncol(List_XY[[k]]) - 1 # the column in List_XY that contains the survival time
                                                            yc_col = which(colnames(List_XY[[k]]) %in% "Censored")   # ncol(List_XY[[k]]) # the column in List_XY that contains the censor time
                                                            
                                                            LogLikelihood(x_ij = List_XY[[k]][i, -c(y_col, yc_col)], 
                                                                          y_ij = List_XY[[k]][i, y_col],  # for each patient in the test set
                                                                          y_ijc = List_XY[[k]][i, yc_col], # ith observation, censor time
                                                                          betas_it = betas[[k]][[j]],  # kth cancer type, jth iteration posterior
                                                                          sigma2_it = sigma2[j])
                                                          }
                                       )))
    
  }
  
  #EFL: here is revised code to aggregate the log-likelihoods
  #EFL: organize log-likelihoods into M by N matrix, where M is number of test samples and N is number of iterations
  logLike.mat <- matrix(nrow = sum(n_vec), ncol = length(ind_to_use))
  for (j in ind_to_use) { # corresponds to every 10th sample 
    k = which(j == ind_to_use) # need to explicitly state the column number because the values of j do not correspond
    logLike.mat[,k] <- unlist(posterior.vec[[j]])
  }
  #EFL: compute joint log-likelihoods, at each iteration, by summing the logs
  logLike.joint <- colSums(logLike.mat)
  
  return(logLike.joint)
}


ComputePosteriorLikelihoodForEveryModel = function(PosteriorsForEveryModel, TestDesign, TestResponse, TestCensored, 
                                                   iters, burnin, model_types) {
  # Computes the posterior likelihood for every model in the simulations
  n_cancer <- length(TestDesign)
  
  # Number of models to consider
  n_models = length(model_types)
  
  # Vector to store all squared deviations in
  pls = c()
  
  for (i in 1:n_models) {
    mod = model_types[i]
    
    PosteriorsForModel_i = PosteriorsForEveryModel[[mod]] # For testing: posteriors = PosteriorsForModel_i
    
    if (mod == "null_model") { # if we are computing the posterior likelihood for the null model
      TestDesign_Null = lapply(TestDesign, function(type) type[, 1, drop = FALSE])
      
      pls[i] = PosteriorLikelihoodForSimulations(PosteriorsForModel_i, TestDesign_Null, TestResponse, TestCensored,
                                                 iters, burnin = iters/2)
    } else if (mod == "joint_model") {
      pls[i] = PosteriorLikelihoodForSimulations(PosteriorsForModel_i, TestDesign, TestResponse, TestCensored,
                                                 iters, burnin = iters/2, joint = TRUE)
      
    } else if (mod == "separate_model") {
      joint_ll_for_each_cancer <- lapply(1:n_cancer, function(type) {
        JointLogLikelihoodForSeparateModelForSimulations(posteriors = PosteriorsForModel_i[[type]],
                                           X = TestDesign[type],
                                           Y = TestResponse[type],
                                           Censored = TestCensored[type], 
                                           burnin = iters/2)
      })
      iters_used_in_ll <- (iters/2/10) + 1
      
      pls[i] <- -log(iters_used_in_ll) + logSum(Reduce('+', joint_ll_for_each_cancer))
      
    } else if (mod == "horseshoe_model") {
      ind_for_betas_by_cancer <- seq(1, p*n_cancer + p, by = p)
      betas_horseshoe <- lapply(1:n_cancer, function(type) { # for each cancer type
        betas_type <- PosteriorsForModel_i$BetaSamples[ind_for_betas_by_cancer[type]:(ind_for_betas_by_cancer[(type+1)]-1),]
        betas_type_list <- as.list(data.frame(betas_type)) # converting each column to a list element
        
        # Checking that converting the matrix to a list worked
        # match <- c(); for (i in 1:(iters/2)) match[i] <- all.equal(betas_type_list[[i]], betas_type[,i]); all(match)
        
        names(betas_type_list) <- NULL
        betas_type_list
      })
      pls[i] <- PosteriorLikelihoodForSimulations(posteriors = list(betas = betas_horseshoe, 
                                                                    sigma2 = PosteriorsForModel_i$Sigma2Samples),
                                                  X = CombineCancerDataTogether(TestDesign, covariates_in_model), 
                                                  Y = TestResponse, 
                                                  Censored = TestCensored, iters = iters/2, burnin = 1)
      
    }
    else { # for any other model
      pls[i] = PosteriorLikelihoodForSimulations(PosteriorsForModel_i, TestDesign, TestResponse, TestCensored,
                                                 iters, burnin = iters/2)
      
    }
  }
  
  names(pls) = model_types
  
  return(pls)
}

CreateInclusionMatrixForCovariates = function(Posteriors, covariates_by_cancer, covariates_in_model, iters) {
  # Returns an inclusion matrix for the covariates based on the posterior draws from the 
  # Gibbs sampler. 
  
  thinned_iters = seq(1,(iters/2), by = 10) # for thinning 
  gamma.gibbs = lapply(1:n_cancer, function(type) { # combine the inclusion indicators into a list of matrices for each cancer type so we can take column means
    do.call(rbind, Posteriors$gamma[[type]])[(iters/2):iters,][thinned_iters,]
  })
  
  inclusion.gibbs = matrix(nrow = length(cancer_types), ncol = p)
  for (i in 1:length(cancer_types)) {
    cov = covariates_by_cancer[[i]]
    n_cov = length(cov)
    for (j in 1:n_cov) {
      k = which(covariates_in_model == cov[j]) # index for covariate in overall inclusion matrix
      inclusion.gibbs[i, k] = mean(gamma.gibbs[[i]][,j]) # entry for current covariate in current cancer type in inclusion matrix
    }
  }
  
  ##adding column and row names to each inclusion matrix
  colnames(inclusion.gibbs) = covariates_in_model
  rownames(inclusion.gibbs) = cancer_types
  
  return(inclusion.gibbs)
}

CredibleIntervalForEveryParamGibbsSampler = function(Posteriors, iters, covariates_in_model) {
  # Returns credible intervals for each parameter in the model and a heatmap for the inclusion indicators
  
  CredibleIntervals = lapply(1:length(Posteriors), function(i) list())
  names(CredibleIntervals) = names(Posteriors)
  
  # Save the parameters
  thinned_iters = seq(1,(iters/2), by = 10) # for thinning 
  sigma2.gibbs = Posteriors$sigma2[(iters/2):iters][thinned_iters]
  lambda2.gibbs = Posteriors$lambda2[(iters/2):iters,][thinned_iters,]
  betatilde.gibbs = Posteriors$beta_tilde[(iters/2):iters,][thinned_iters,]
  betas.gibbs = lapply(1:n_cancer, function(type) {
    do.call(rbind, Posteriors$betas[[type]])[(iters/2):iters,][thinned_iters,]
  })
  gamma.gibbs = lapply(1:n_cancer, function(type) {
    do.call(rbind, Posteriors$gamma[[type]])[(iters/2):iters,][thinned_iters,]
  })
  
  #Sigma2
  CredibleIntervals$sigma2 = CredibleInterval(sigma2.gibbs)
  
  #Lambda2
  Lambda2CredibleIntervals = data.frame(GibbsLower = double(),
                                        GibbsUpper = double())
  for (i in 1:p) {
    intervals_vec = CredibleInterval(lambda2.gibbs[,i])
    Lambda2CredibleIntervals[i,] = intervals_vec
  }
  
  CredibleIntervals$lambda2 = Lambda2CredibleIntervals
  
  #BetaTilde
  BetaTildeCredibleIntervals = data.frame(GibbsLower = double(),
                                          GibbsUpper = double())
  for (i in 1:p) {
    intervals_vec = CredibleInterval(betatilde.gibbs[,i])
    BetaTildeCredibleIntervals[i,] = intervals_vec
  }
  
  CredibleIntervals$beta_tilde = BetaTildeCredibleIntervals
  
  #Betas
  BetaCredibleIntervals = lapply(1:length(cancer_types), function(i) list())
  names(BetaCredibleIntervals) = cancer_types
  
  for (i in 1:n_cancer) {
    covariates_i = covariates_by_cancer[[i]]
    n_covariate_i = length(covariates_i)
    BetaGibbsCredibleInterval_Cancer_i = t(sapply(1:n_covariate_i, function(k) {
      CredibleInterval(betas.gibbs[[i]][,k])
    }))
    colnames(BetaGibbsCredibleInterval_Cancer_i) = NULL
    BetaCredibleIntervals[[i]] = BetaGibbsCredibleInterval_Cancer_i
  }
  
  CredibleIntervals$betas = BetaCredibleIntervals
  
  #Gammas
  inclusion.gibbs = CreateInclusionMatrixForCovariates(Posteriors, covariates_in_model)
  
  ##convert inclusion matrix into a tibble for ggplot
  ##make a heatmap of the estimated inclusion probability
  
  inclusion.gibbs = as_tibble(inclusion.gibbs)
  x.gibbs = gather(data = inclusion.gibbs, key = pc, value = inclusion, 1:p)
  x.gibbs$cancer = cancer_types
  x.gibbs$pc = factor(x.gibbs$pc, levels = covariates_in_model)
  
  GammaPlot <- ggplot(data = x.gibbs, mapping = aes(x = pc, y = cancer, fill = inclusion)) +
    geom_tile() +
    xlab(label = "PC") +
    ggtitle("Inclusion Heatmap, Gibbs Sampler, Intercept Outside SS") + 
    geom_text(aes(label = round(inclusion,1)))
  
  CredibleIntervals$gamma = GammaPlot
  
  return(CredibleIntervals)
}

PosteriorLikelihood = function(posteriors, X, Y, Censored, current_train_ind, burnin, joint = FALSE) {
  # posteriors is the result from the model computation
  # X and Y are the datasets are Covariates and Survival datasets, respectively
  # Censored is the censor time for each subject
  # burnin is the number of iterations from the start to throw away
  # joint is a boulean indicator for whether this is for the joint model or not 
  logSum <- function(l) {
    # given l = log(v), where v is a vector,
    # computes log(sum(v))
    return(max(l) + log(sum(exp(l-max(l)))))
  }
  
  LogLikelihood = function(x_ij, y_ij, y_ijc, betas_it, sigma2_it) {
    # calculate the likelihood for just one set of posteriors
    mu_ij = sum(x_ij * betas_it) 
    
    LogNormalLikelihood = function(y, mu, sig2) {
      # Calculates log density of normal 
      # nl = (1/sqrt(2*pi*sig2))*exp(-(1/(2*sig2))*(y - mu)^2)
      lnl = (-0.5)*log(2*pi*sig2) - (1/(2*sig2))*(y - mu)^2
      return(lnl)
    }
    
    if (!is.na(y_ij)) { # y_ij is not censored
      post_ij = LogNormalLikelihood(log(y_ij), mu = mu_ij, sig2 = sigma2_it) # dnorm(log(y_ij), mean = mu_ij1, sd = sqrt(sigma21_it)) # sample from log-normal
      
    } else { # if y_ij is censored, use censoring time
      post_ij = log(1 - pnorm((log(y_ijc) - mu_ij)/sqrt(sigma2_it), mean = 0, sd = 1)) 
    }
    
    return(post_ij)
  }

  # Extracting the necessary information
  n_vec = sapply(X, nrow) # number of patients per cancer type
  iters = length(posteriors$betas[[1]]) # changed this from nrow(beta_tilde) because not every model has a beta_tilde
  Training_obs_cv_iter = current_train_ind
  betas = posteriors$betas
  sigma2 = posteriors$sigma2
  Test_obs = lapply(seq(length(X)), function(i) seq(n_vec[i])[!(seq(n_vec[i]) %in% Training_obs_cv_iter[[i]])]) # for each cancer type (seq(length(X)) gives number of cancers since it is a list), obtain the subject indices not contained in the training data set indices
  
  # Obtaining test observations and their times of last contact, etc
  Y_test = lapply(seq(length(X)), function(i) Y[[i]][Test_obs[[i]]]) 
  X_test = lapply(seq(length(X)), function(i) X[[i]][Test_obs[[i]], ])
  Censored_test = lapply(seq(length(X)), function(i) Censored[[i]][Test_obs[[i]]])
  
  # If this is the joint model, combine the test data together:
  if (joint) {
    # Combine the test data into one list for the joint model
    X_test <- list(do.call(rbind, CombineCancerDataTogether(X_test, covariates_in_model_full)))
    Y_test <- list(unlist(Y_test))
    Censored_test <- list(unlist(Censored_test))
  } 
  
  # Combine the test data together
  List_XY = lapply(seq(length(X_test)), function(i) cbind(X_test[[i]], "Survival" = Y_test[[i]], "Censored" = Censored_test[[i]]))
  
  posterior.vec = list()
  ind_to_use = seq(burnin, iters, by = 10)
  
  for (j in ind_to_use) { ##EFL: changed index 'i' to 'j' here
    
    # for each cancer type
    # for each row in that cancer type
    # apply Likelihood with the corresponding cancer type's posterior parameters
    # at the jth iteration
    #EFL: also changed here 'i' to 'j' on left-hand side below
    posterior.vec[[j]] = unlist(lapply(seq(length(X_test)), # change from X to X_test
                                       function(k) sapply(seq(nrow(List_XY[[k]])), # for each cancer type
                                       function(i) { # for each test subject in cancer type k
                                         y_col = which(colnames(List_XY[[k]]) %in% "Survival")  # ncol(List_XY[[k]]) - 1 # the column in List_XY that contains the survival time
                                         yc_col = which(colnames(List_XY[[k]]) %in% "Censored")   # ncol(List_XY[[k]]) # the column in List_XY that contains the censor time
                                         
                                         LogLikelihood(x_ij = List_XY[[k]][i, -c(y_col, yc_col)], 
                                                       y_ij = List_XY[[k]][i, y_col],  # for each patient in the test set
                                                       y_ijc = List_XY[[k]][i, yc_col], # ith observation, censor time
                                                       betas_it = betas[[k]][[j]],  # kth cancer type, jth iteration posterior
                                                       sigma2_it = sigma2[j])
                                         }
                                       )))
    
  }
  
  #EFL: here is revised code to aggregate the log-likelihoods
  #EFL: organize log-likelihoods into M by N matrix, where M is number of test samples and N is number of iterations
  logLike.mat <- matrix(nrow = length(unlist(Test_obs)), ncol = length(ind_to_use))
  for (j in ind_to_use) { # corresponds to every 10th sample 
    k = which(j == ind_to_use) # need to explicitly state the column number because the values of j do not correspond
    logLike.mat[,k] <- unlist(posterior.vec[[j]])
  }
  #EFL: compute joint log-likelihoods, at each iteration, by summing the logs
  logLike.joint <- colSums(logLike.mat)
  #EFL: Compute log of full posterior predictive likelihood, using logSum
  PL_Survival<- -log(length(ind_to_use)) + logSum(logLike.joint) 
  
  return(PL_Survival)
}

logSum <- function(l) {
  # given l = log(v), where v is a vector,
  # computes log(sum(v))
  return(max(l) + log(sum(exp(l-max(l)))))
}

JointLogLikelihoodForSeparateModel = function(posteriors, X, Y, Censored, current_train_ind, burnin) {
  # For computing the joint log-likelihood for the separated model
  # posteriors = posterior samples for a single cancer type
  # X = Covariates data for a single cancer type
  # Y = Survival times for a single cancer type
  # Censored = list of censored times for each cancer type
  logSum <- function(l) {
    # given l = log(v), where v is a vector,
    # computes log(sum(v))
    return(max(l) + log(sum(exp(l-max(l)))))
  }
  
  LogLikelihood = function(x_ij, y_ij, y_ijc, betas_it, sigma2_it) {
    # calculate the likelihood for just one set of posteriors
    mu_ij = sum(x_ij * betas_it) 
    
    LogNormalLikelihood = function(y, mu, sig2) {
      # Calculates log density of normal 
      # nl = (1/sqrt(2*pi*sig2))*exp(-(1/(2*sig2))*(y - mu)^2)
      lnl = (-0.5)*log(2*pi*sig2) - (1/(2*sig2))*(y - mu)^2
      return(lnl)
    }
    
    if (!is.na(y_ij)) { # y_ij is not censored
      post_ij = LogNormalLikelihood(log(y_ij), mu = mu_ij, sig2 = sigma2_it) # dnorm(log(y_ij), mean = mu_ij1, sd = sqrt(sigma21_it)) # sample from log-normal
      
    } else { # if y_ij is censored, use censoring time
      post_ij = log(1 - pnorm((log(y_ijc) - mu_ij)/sqrt(sigma2_it), mean = 0, sd = 1)) 
    }
    
    return(post_ij)
  }
  
  # Extracting the necessary information
  n_vec = sapply(X, nrow) # number of patients per cancer type
  iters = length(posteriors$betas[[1]]) # changed this from nrow(beta_tilde) because not every model has a beta_tilde
  Training_obs_cv_iter = current_train_ind
  betas = posteriors$betas
  sigma2 = posteriors$sigma2
  Test_obs = lapply(seq(length(X)), function(i) seq(n_vec[i])[!(seq(n_vec[i]) %in% Training_obs_cv_iter[[i]])]) # for each cancer type (seq(length(X)) gives number of cancers since it is a list), obtain the subject indices not contained in the training data set indices
  
  # Obtaining test observations and their times of last contact, etc
  Y_test = lapply(seq(length(X)), function(i) Y[[i]][Test_obs[[i]]]) 
  X_test = lapply(seq(length(X)), function(i) X[[i]][Test_obs[[i]], ])
  Censored_test = lapply(seq(length(X)), function(i) Censored[[i]][Test_obs[[i]]])
  
  # Combine the test data together
  List_XY = lapply(seq(length(X_test)), function(i) cbind(X_test[[i]], "Survival" = Y_test[[i]], "Censored" = Censored_test[[i]]))
  
  posterior.vec = list()
  ind_to_use = seq(burnin, iters, by = 10)
  
  for (j in ind_to_use) { ##EFL: changed index 'i' to 'j' here
    
    # for each cancer type
    # for each row in that cancer type
    # apply Likelihood with the corresponding cancer type's posterior parameters
    # at the jth iteration
    #EFL: also changed here 'i' to 'j' on left-hand side below
    posterior.vec[[j]] = unlist(lapply(seq(length(X_test)), # change from X to X_test
                                       function(k) sapply(seq(nrow(List_XY[[k]])), # for each cancer type
                                                          function(i) { # for each test subject in cancer type k
                                                            y_col = which(colnames(List_XY[[k]]) %in% "Survival")  # ncol(List_XY[[k]]) - 1 # the column in List_XY that contains the survival time
                                                            yc_col = which(colnames(List_XY[[k]]) %in% "Censored")   # ncol(List_XY[[k]]) # the column in List_XY that contains the censor time
                                                            
                                                            LogLikelihood(x_ij = List_XY[[k]][i, -c(y_col, yc_col)], 
                                                                          y_ij = List_XY[[k]][i, y_col],  # for each patient in the test set
                                                                          y_ijc = List_XY[[k]][i, yc_col], # ith observation, censor time
                                                                          betas_it = betas[[k]][[j]],  # kth cancer type, jth iteration posterior
                                                                          sigma2_it = sigma2[j])
                                                          }
                                       )))
    
  }
  
  #EFL: here is revised code to aggregate the log-likelihoods
  #EFL: organize log-likelihoods into M by N matrix, where M is number of test samples and N is number of iterations
  logLike.mat <- matrix(nrow = length(unlist(Test_obs)), ncol = length(ind_to_use))
  for (j in ind_to_use) { # corresponds to every 10th sample 
    k = which(j == ind_to_use) # need to explicitly state the column number because the values of j do not correspond
    logLike.mat[,k] <- unlist(posterior.vec[[j]])
  }
  #EFL: compute joint log-likelihoods, at each iteration, by summing the logs
  logLike.joint <- colSums(logLike.mat)
  
  return(logLike.joint)
}

CombineCancerDataTogether = function(Covariates, covariates_in_model) {
  # Fill in 0s for all missing values. 
  
  # Initializing a new Covariates list to store the filled-in data.
  CovariatesFilledIn <- Covariates
  names(CovariatesFilledIn) <- cancer_types
  
  # Iterating through the cancers and filling in the missing covariates with 0s. 
  for (cancer in 1:n_cancer) {
    # Storing the current cancer
    current_cancer <- CovariatesFilledIn[[cancer]]
    
    # Checking which covariates are available for this cancer time
    current_avail <- covariates_in_model %in% colnames(current_cancer)
    
    # Creating a new matrix to store the filled in data in
    current_cancer_filled_in <- matrix(nrow = nrow(current_cancer), ncol = length(covariates_in_model))
    colnames(current_cancer_filled_in) <- covariates_in_model
    rownames(current_cancer_filled_in) <- rownames(current_cancer)
    
    # Filling in the 0s
    current_cancer_filled_in[, !current_avail] <- 0
    
    # Filling in the available data
    current_cancer_filled_in[, current_avail] <- current_cancer
    
    # Saving to the new list
    CovariatesFilledIn[[cancer]] <- current_cancer_filled_in
  }

  # Return
  CovariatesFilledIn
}

PrepareCancerDataForPanCanVarSel <- function(Covariates, covariates_in_model) {
  # Prepares the covariates in the way that the PanCanVarSel package expects
  # Initializing a new Covariates list to store the filled-in data.
  CovariatesFilledIn <- Covariates
  names(CovariatesFilledIn) <- cancer_types
  
  for (cancer in 1:n_cancer) {
    # Storing the current cancer
    current_cancer <- CovariatesFilledIn[[cancer]]
    
    # Checking which covariates are available for this cancer time
    current_avail <- covariates_in_model %in% colnames(current_cancer)
    
    # Creating a new matrix to store the filled in data in
    current_cancer_filled_in <- matrix(nrow = nrow(current_cancer), ncol = length(covariates_in_model))
    colnames(current_cancer_filled_in) <- covariates_in_model
    rownames(current_cancer_filled_in) <- rownames(current_cancer)
    
    # Filling in the 0s
    current_cancer_filled_in[, !current_avail] <- 0
    
    # Filling in the available data
    current_cancer_filled_in[, current_avail] <- current_cancer
    
    # Saving to the new list
    CovariatesFilledIn[[cancer]] <- current_cancer_filled_in
  }
  
  # Return
  CovariatesFilledIn.BDiag <- bdiag(CovariatesFilledIn)
  as.matrix(CovariatesFilledIn.BDiag)
}
