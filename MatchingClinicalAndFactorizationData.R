### PanTCGA Data Manipulation: Using the patient IDs from the pan.fac.results.rda file, 
### match the patients with their survival times and censor times. Store these data
### in lists to be used in the hierarchical model. 
### Note: the Age column is labeled "0.5" so the numeric order is maintained for the covariates
### Author: Sarah Samorodnitsky and Eric Lock (2020)
### University of Minnesota

source("HelperFunctions.R")

# Loading in the data
clinical_data = read.csv("TCGA-CDR.csv", header = T) # the clinical data
load("pan.fac.results_v2.rda")
load("mod.pcs.v3.rda") # this is the filtered data based on the method Eric and I discussed

cancer_types = sort(unique(unlist(tum.type.list))) # should be 29 cancer types

# Matching the patient IDs from the integrative factorization data to the clinical data from TCGA
MatchBarCodes = function(sample.id.list, mod.pcs, tum.type.list, clinical_data, cancer_types) {
  # sample.id.list (list of characters): contains a list of patient IDs for each module
  # clinical_data (matrix): contains the clinical data from TCGA on each patient, including a column with patient IDs
  # mod.pcs (list of vectors or list of matrices): a list with length equal to the number of modules where each
  # entry in the list is the PCs for that module to be used in prediction. 
  
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
  
  # Extracting the clinical data that corresponds to the patients used in the factorization
  unique.ids = unique(unlist(sample.id.list)) # Should be 6973 IDs - matches the number of tumor samples listed in the manuscript
  clinical_data$bcr_patient_barcode = as.character(clinical_data$bcr_patient_barcode) # changing from factor -> character to match the type of the factorization data IDs
  clin.data.for.ids = clinical_data[clinical_data$bcr_patient_barcode %in% unique.ids, ] # the clinical data for the IDs found in the factorization
  
  n_mod = length(sample.id.list) # the number of modules
  n_cancer = length(cancer_types) # number of cancers used in analysis
  n_vec = table(tum.type.list[[1]]) # use the first entry in tum.type.list because that module used every patient
  
  # iterate through the modules
  # for each module 
  #     for each cancer type in that module: 
  #         cbind the PC values for the patients with that cancer type to the data matrix 
  #         store the PC values in the corresponding list
  
  X = InitManyMatrices(n_cancer, n_vec); names(X) = cancer_types
  Y = InitManyLists(n_cancer); names(Y) = cancer_types
  C = InitManyLists(n_cancer); names(C) = cancer_types
  
  # Gathering all the PCs together for each cancer type
  for (i in 1:n_mod) { # for each module
    cancer_i = unique(tum.type.list[[i]]) # cancers used in module i
    tum.type.i = tum.type.list[[i]] # the cancer type of each patient used in module i
    n_cancer_i = length(cancer_i) # number of cancers used in module i
    pcs_i = mod.pcs[[i]] # the principal components for module i
    ids.i = sample.id.list[[i]] # patient IDs used to create this module
    
    for (j in 1:n_cancer_i) { # for each cancer in module i
      cancer_ij = cancer_i[j] 
      pcs_cancer_ij = pcs_i[tum.type.i == cancer_ij,, drop = FALSE] # the entries in the PCs corresponding to the jth cancer type used in this module
      temp_X = cbind(X[[cancer_ij]], pcs_cancer_ij) # append the corresponding entries for the jth cancer type in the PCs for the ith module. Store temporarily here
      rownames(temp_X) = ids.i[tum.type.i == cancer_ij] # name each entry with the corresponding patient ID for matching with the clinical data
      X[[cancer_ij]] = temp_X
    }
  }
  
  # Reordering the patients so the barcodes match between the covariates and the clinical data. 
  # Also adding an intercept to the covariate matrices for each cancer type
  for (i in 1:n_cancer) {
    # Reordering the patients so they match between X and Y
    cancer_i = cancer_types[i] # the ith cancer type
    
    if (cancer_i == "CORE") { # contraction of COAD and READ
      temp_YC = clin.data.for.ids[clin.data.for.ids$type %in% c("COAD", "READ"), ] # survival, censor data for the ith cancer type
      inds = c()
      ids.fac = rownames(X[[i]])
      clin.ids = temp_YC$bcr_patient_barcode # ids in the clinical data
      for (j in 1:length(ids.fac)) { # manually sorting the patient IDs
        inds[j] = which(ids.fac[j] == clin.ids) # locating the factorization id in the clinical data
      }
      temp_YC.2 = temp_YC[inds, ] # now the clinical data is in the same order as the original order of the factorization data
      
    } else {
      temp_YC = clin.data.for.ids[clin.data.for.ids$type %in% cancer_i, ] # survival, censor data for the ith cancer type
      ids.x = rownames(X[[i]]) # the order of patient IDs in the covariate matrix (the factorization data)
      ord.ind = match(temp_YC$bcr_patient_barcode, ids.x) # ordering the patient IDs in the clinical data to match the patient IDs used in factorization (i.e. the order of the data in the PCs) (gives positions of matches to patient IDs in covariates to patient IDs in clinical data so we can reorder the clinical data to match the covariate data)
      temp_YC.2 = temp_YC[ord.ind, ] # reordering the clinical data so that patient bar codes match that of the factorization data
    }
    
    # Centering age
    # This may throw a warning because some age values will be NAs
    age.i = as.numeric(as.character(temp_YC.2$age_at_initial_pathologic_diagnosis)) 
    
    # Adding an intercept and age as a column
    Xi_WithAge = cbind(1, age.i, X[[i]])
    colnames(Xi_WithAge) = c("0", "0.5", colnames(Xi_WithAge)[-(1:2)]) # adding column names, getting rid of the first two column names to begin with because they are meaningless
    
    # Storing the new covariates matrix back into the X list
    X[[i]] = Xi_WithAge
    
    # will give NA warning because of censoring
    Y[[cancer_i]] = as.numeric(as.character(temp_YC.2$death_days_to)) # survival time
    names(Y[[cancer_i]]) = temp_YC.2$bcr_patient_barcode
    C[[cancer_i]] = as.numeric(as.character(temp_YC.2$last_contact_days_to)) # censor time
    names(C[[cancer_i]]) = temp_YC.2$bcr_patient_barcode
  }
  
  # Removing observations that are missing survival times and censor times,
  # obs that have negative survival time,
  # obs that have 0 survival time,
  # obs that have negative censor time,
  # and obs that have a missing age value 
  for (i in 1:n_cancer) {
    X_i = X[[i]]
    Y_i = Y[[i]]
    C_i = C[[i]]
    
    no_surv_cens = which(is.na(Y_i) & is.na(C_i))
    if (length(no_surv_cens) != 0) { # if there are observations missing both a survival time and censor time
      X_i = X_i[-no_surv_cens, ]
      Y_i = Y_i[-no_surv_cens]
      C_i = C_i[-no_surv_cens]
    }
    
    neg_surv = which(Y_i <= 0)
    if (length(neg_surv) != 0) { # if there are observations with a negative or 0 survival time
      X_i = X_i[-neg_surv, ]
      Y_i = Y_i[-neg_surv]
      C_i = C_i[-neg_surv]
    }
    
    neg_cens = which(C_i <= 0) 
    if (length(neg_cens) != 0) {
      X_i = X_i[-neg_cens, ]
      Y_i = Y_i[-neg_cens]
      C_i = C_i[-neg_cens]
    }
    
    missing_age = which(is.na(X_i[, "0.5"]))
    if (length(missing_age) != 0) { # if any observations have missing age value
      X_i = X_i[-missing_age, ]
      Y_i = Y_i[-missing_age]
      C_i = C_i[-missing_age]
    }
    
    # return updated data
    X[[i]] = X_i
    Y[[i]] = Y_i
    C[[i]] = C_i
  }
  
  # Standardizing the variables after removing the observations that will not be kept in the analysis
  for (i in 1:n_cancer) {
    X_i = X[[i]]
    
    # standardizing age
    age.centered.i = (X_i[, "0.5"] - mean(X_i[, "0.5"]))/sd(X_i[, "0.5"])
    
    # standardizing the factorization covariates
    factorization.pcs.i = X_i[, !(colnames(X_i) %in% c("0", "0.5"))]
    factorization.pcs.scaled.i = apply(factorization.pcs.i, 2, function(col) {
      (col - mean(col))/sd(col)
    })
    
    # putting the scaled covariates back together
    Xi_WithScaledPredictors = cbind(X_i[, 1], age.centered.i, factorization.pcs.scaled.i)
    colnames(Xi_WithScaledPredictors) = c("0", "0.5", colnames(Xi_WithScaledPredictors)[-(1:2)]) # adding column names, getting rid of the first two column names to begin with because they are meaningless
    
    # Restoring them in X
    X[[i]] = Xi_WithScaledPredictors
    
  }
  
  return(list(X = X, Y = Y, C = C))
}

# Will give warnings about introducing NAs

XY = MatchBarCodes(sample.id.list, mod.pcs, tum.type.list, clinical_data, cancer_types)
Covariates = XY$X
Survival = XY$Y 
Censored = XY$C

save(Covariates, Survival, Censored, file = "XYC_V2_WithAge_StandardizedPredictors.rda")
