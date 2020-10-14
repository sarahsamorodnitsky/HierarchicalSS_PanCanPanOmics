### Checking that the results from matching the clinical and factorization data worked properly
### Author: Sarah Samorodnitsky and Eric Lock (2020)
### University of Minnesota

# Either source this script
source("MatchingClinicalAndFactorizationData.R") # will give warnings because of NAs 
# or
load("mod.pcs.v3.rda")
load("XYC_V2_WithAge.rda") # this is the data from MatchingClinicalAndFactorizationData but non-standardized so we can check if it matches. Comment out the standardization portion of that script to do this
load('pan.fac.results_v2.rda') # comparing the reordered data to the original
clinical_data = read.csv("TCGA-CDR.csv", header = T) # the clinical data
source("HelperFunctions.R")

### Checking the way the data was matched (this is the body of the matching function): 

# Extracting the clinical data that corresponds to the patients used in the factorization
unique.ids = unique(unlist(sample.id.list)) # Should be 6973 IDs - matches the number of tumor samples listed in the manuscript
clinical_data$bcr_patient_barcode = as.character(clinical_data$bcr_patient_barcode) # changing from factor -> character to match the type of the factorization data IDs
clin.data.for.ids = clinical_data[clinical_data$bcr_patient_barcode %in% unique.ids, ] # the clinical data for the IDs found in the factorization

n_mod = length(sample.id.list) # the number of modules
n_cancer = length(cancer_types) # number of cancers used in analysis
n_vec = table(tum.type.list[[1]]) # use the first entry in tum.type.list because that module used every patient

# initializing these variables
X = InitManyMatrices(n_cancer, n_vec); names(X) = cancer_types
Y = InitManyLists(n_cancer); names(Y) = cancer_types
C = InitManyLists(n_cancer); names(C) = cancer_types

# Gathering all the PCs together for each cancer type
# The PCs are organized by module. We need to organize them by cancer type. 
# So we extract each value in each PC by the cancer type it originated from
# and place that in its own separate matrix. This matrix will have a consistent
# number of rows because the same patients for each cancer type were used to 
# create the corresponding modules. The PC is the right singular vector which
# has the same length as the number of patients times the corresponding singular
# value. 
for (i in 1:n_mod) { # for each module
  cancer_i = unique(tum.type.list[[i]]) # cancers used in module i
  tum.type.i = tum.type.list[[i]] # the cancer type of each patient used in module i
  n_cancer_i = length(cancer_i) # number of cancers used in module i
  pcs_i = mod.pcs[[i]] # the principal components for module i
  ids.i = sample.id.list[[i]] # patient IDs used to create this module
  
  for (j in 1:n_cancer_i) { # for each cancer in module i
    cancer_ij = cancer_i[j] 
    pcs_cancer_ij = pcs_i[tum.type.i == cancer_ij,, drop = FALSE] # the entries in the PCs corresponding to the jth cancer type used in this module, drop keeps this variable as a column
    temp_X = cbind(X[[cancer_ij]], pcs_cancer_ij) # append the corresponding entries for the jth cancer type in the PCs for the ith module. Store temporarily here. If X[[cancer_ij]] is empty, then it will just create the first column
    rownames(temp_X) = ids.i[tum.type.i == cancer_ij] # name each entry with the corresponding patient ID for matching with the clinical data
    X[[cancer_ij]] = temp_X
  }
}

# Reordering the patients so the barcodes match between the covariates and the clinical data. 
# Also adding an intercept to the covariate matrices for each cancer type

# One thing I've changed here - I don't center age to make sure it matches with the original
# clinical data. 

same_order = c()

for (i in 1:n_cancer) {
  # Reordering the patients so they match between X and Y
  cancer_i = cancer_types[i] # the ith cancer type
  
  if (cancer_i == "CORE") {
    temp_YC = clin.data.for.ids[clin.data.for.ids$type %in% c("COAD", "READ"), ] # survival, censor data for the ith cancer type
    inds = c()
    ids.fac = rownames(X[[i]])
    clin.ids = temp_YC$bcr_patient_barcode
    for (j in 1:length(ids.fac)) { # manually sorting the patient IDs. Could also have used match()
      inds[j] = which(ids.fac[j] == clin.ids)
    }
    temp_YC.2 = temp_YC[inds, ]
    same_order[i] = all(temp_YC.2$bcr_patient_barcode == ids.fac)
    
  } else {
    temp_YC = clin.data.for.ids[clin.data.for.ids$type %in% cancer_i, ] # survival, censor data for patients in the ith cancer type
    ids.x = rownames(X[[i]]) # the order of patient IDs in the covariate matrix (the factorization data)
    ord.ind = match(ids.x, temp_YC$bcr_patient_barcode) # ordering the patient IDs in the clinical data to match the patient IDs used in factorization (i.e. the order of the data in the PCs) (gives positions of matches to patient IDs in covariates to patient IDs in clinical data so we can reorder the clinical data to match the covariate data)
    temp_YC.2 = temp_YC[ord.ind, ] # reordering the clinical data so that patient bar codes match that of the factorization data
    same_order[i] = all(temp_YC.2$bcr_patient_barcode == ids.x)
  }
  
  # Saving age
  # This may throw a warning because some age values will be NAs
  age.i = as.numeric(as.character(temp_YC.2$age_at_initial_pathologic_diagnosis)) 
  
  # Adding an intercept and age as a column
  Xi_WithAge = cbind(1, age.i, X[[i]])
  colnames(Xi_WithAge) = c("0", "Age", colnames(Xi_WithAge)[-(1:2)]) # adding column names, getting rid of the first two column names to begin with because they are meaningless
  
  # Storing the new covariates matrix back into the X list
  X[[i]] = Xi_WithAge
  
  Y[[cancer_i]] = as.numeric(as.character(temp_YC.2$death_days_to)) # survival time, will give warning about NAs
  names(Y[[cancer_i]]) = temp_YC.2$bcr_patient_barcode
  C[[cancer_i]] = as.numeric(as.character(temp_YC.2$last_contact_days_to)) # censor time
  names(C[[cancer_i]]) = temp_YC.2$bcr_patient_barcode
}

all(same_order)

# At this point, all patient ids in Y_i should match the factorization ids

# Removing observations that are missing survival times and censor times, obs that have negative survival time, 
# and obs that have 0 survival time. Then checking that the lengths match. Not checking order here

num_surv_NA_matches = c()
num_cens_NA_matches = c()
surv_time_matches = c()
cens_time_matches = c()
lengths_match = c()
age_matches = c()
no_na_age = c()

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
  
  missing_age = which(is.na(X_i[, "Age"]))
  if (length(missing_age) != 0) { # if any observations have missing age value
    X_i = X_i[-missing_age, ]
    Y_i = Y_i[-missing_age]
    C_i = C_i[-missing_age]
  }
  
  # return updated data
  X[[i]] = X_i
  Y[[i]] = Y_i
  C[[i]] = C_i
  
  # making sure the lengths match
  lengths_match[i] = (length(Y[[i]]) == length(C[[i]])) & (length(Y[[i]]) == nrow(X[[i]]))
  
  # Checking matches
  cancer_i = cancer_types[i] # the ith cancer type
  
  if (cancer_i == "CORE") {
    temp_YC = clin.data.for.ids[clin.data.for.ids$type %in% c("COAD", "READ"), ] # survival, censor data for the ith cancer type
  } else {
    temp_YC = clin.data.for.ids[clin.data.for.ids$type %in% cancer_i, ] # survival, censor data for the ith cancer type
  }
  
  dat = temp_YC[temp_YC$bcr_patient_barcode %in% names(Y_i), ] # saving just the Ys that were not removed during processing
  dat = dat[match(names(Y_i), dat$bcr_patient_barcode),] # to ensure the order matches between the clinical data and the reordered data
  
  # this is the clinical data
  surv = dat$death_days_to; names(surv) = dat$bcr_patient_barcode
  cens = dat$last_contact_days_to; names(cens) = dat$bcr_patient_barcode
  age = dat$age_at_initial_pathologic_diagnosis
  
  # checking the survival times match
  num_surv_NA_matches[i] = sum(surv == "#N/A") == sum(is.na(Y_i))   # since I already checked the lengths of Y and C (and X) match, then if the number of NAs matches, the number of non-NAs will match, too. 
  surv_time_matches[i] = all(surv[surv != "#N/A"] == na.omit(Y_i)) # if these all match then the rest must be NAs
  
  # checking the censor times match
  num_cens_NA_matches[i] = sum(cens == "#N/A") == sum(is.na(C_i)) 
  cens_time_matches[i] = all(cens[cens != "#N/A"] == na.omit(C_i))
  
  # checking the ages match
  age_matches[i] = all(age[!is.na(age)] == X_i[, "Age"])
  
  # checking that no ages are missing
  no_na_age[i] = sum(is.na(X_i[,"Age"])) == 0
  
  # Checking that the censor and survival time for a patient ID matches in X, Y, C matches the clinical data
}

all(num_surv_NA_matches)
all(num_cens_NA_matches) 
all(lengths_match)

all(cens_time_matches)
all(surv_time_matches)

all(age_matches)
all(no_na_age)

### Checking that the data looks proper:

# 1. Checking that the dimensions of the data match 
dim_match = c()
for (i in 1:length(X)) {
  dim_match[i] = (nrow(X[[i]]) == length(Y[[i]])) & (nrow(X[[i]]) == length(C[[i]]))
}
all(dim_match)

# 2. Checking that every observation has either a survival time or a censor time and the survival times are not negative or the censor times are not negative if they are censored
prop_surv_cens = c()
all_have_one = c() # all have a non-NA survival time and NA censor time OR NA survival time and non-NA censor time
for (i in 1:length(Y)) {
  prop_surv_cens[i] = any(na.omit(Y[[i]]) <= 0) 
  all_have_one[i] = any(is.na(Y[[i]]) & is.na(C[[i]])) # any((is.na(Y[[i]]) & !is.na(C[[i]])) | (!is.na(Y[[i]]) & is.na(C[[i]])))
}
any(prop_surv_cens) # should all be false
any(all_have_one) # should all be false

# 3. Very explicitly checking that for each patient ID in the final amended dataset, their censor time
# and survival time were not changed by doing as.numeric, etc. 

order_matches = c()
surv_time_correct = c()
cens_time_correct = c()
lengths_match = c()
age_correct = c()

for (i in 1:length(cancer_types)) {
  X_i = X[[i]]
  Y_i = Y[[i]]
  C_i = C[[i]]
  
  cancer_i = cancer_types[i]
  
  # checking that the order of the patient IDs matches across all lists
  order_matches[i] = all((rownames(X_i) == names(Y_i))) & all(rownames(X_i) == names(C_i))
  
  # checking that the lengths match
  lengths_match[i] = (length(Y[[i]]) == length(C[[i]])) & (length(Y[[i]]) == nrow(X[[i]]))
  
  if (cancer_i == "CORE") {
    temp_YC = clin.data.for.ids[clin.data.for.ids$type %in% c("COAD", "READ"), ] # survival, censor data for the ith cancer type
  } else {
    temp_YC = clin.data.for.ids[clin.data.for.ids$type %in% cancer_i, ] # survival, censor data for the ith cancer type
  }

  n_cancer_i = nrow(X_i)
  ids_cancer_i = names(Y_i)
  
  surv_time_matches_i = c()
  cens_time_matches_i = c()
  age_matches_i = c()

  for (j in 1:n_cancer_i) {
    id_j = ids_cancer_i[j]
    surv_time_j = Y_i[j]
    cens_time_j = C_i[j]
    age_j = X_i[j, "Age"]
    
    clin_surv_j = temp_YC[temp_YC$bcr_patient_barcode == id_j, ]$death_days_to
    clin_cens_j = temp_YC[temp_YC$bcr_patient_barcode == id_j, ]$last_contact_days_to
    clin_age_j = temp_YC[temp_YC$bcr_patient_barcode == id_j, ]$age_at_initial_pathologic_diagnosis
    
    clin_surv_j.2 = if(clin_surv_j == "#N/A") NA else clin_surv_j
    clin_cens_j.2 = if(clin_cens_j == "#N/A") NA else clin_cens_j
    clin_age_j.2 = if(clin_age_j == "#N/A") NA else clin_age_j
    
    surv_time_matches_i[j] = (surv_time_j == clin_surv_j.2) | (is.na(surv_time_j) & is.na(clin_surv_j.2))
    cens_time_matches_i[j] = (cens_time_j == clin_cens_j.2) | (is.na(cens_time_j) & is.na(clin_cens_j.2))
    age_matches_i[j] = (age_j == clin_age_j.2) | (is.na(age_j) & is.na(clin_age_j.2))
  }
  
  surv_time_correct[i] = all(surv_time_matches_i)
  cens_time_correct[i] = all(cens_time_matches_i)
  age_correct[i] = all(age_matches_i)
}


all(order_matches)
all(lengths_match)
all(surv_time_correct)
all(cens_time_correct)
all(age_correct)

# Very explicitly checking that each observation matches
all_match = c()

for (i in 1:length(Y)) {
  all_match_i = c()
  for (j in 1:length(Y[[i]])) {
    current_y = Y[[i]][j]
    current_c = C[[i]][j]
    current_id = names(current_y)
    current_age = X[[i]][j,"Age"]
    
    current_clin_data = clin.data.for.ids[clin.data.for.ids$bcr_patient_barcode == current_id, ]
    clin_y_ij = if (current_clin_data$death_days_to == "#N/A") NA else current_clin_data$death_days_to
    clin_c_ij = if (current_clin_data$last_contact_days_to == "#N/A") NA else current_clin_data$last_contact_days_to
    clin_age_ij = if(current_clin_data$age_at_initial_pathologic_diagnosis == "#N/A") NA else current_clin_data$age_at_initial_pathologic_diagnosis
    
    all_match_i[j] = (((current_y == clin_y_ij) & (is.na(current_c) &  is.na (clin_c_ij))) | ((is.na(current_y) & is.na(clin_y_ij)) & (current_c == clin_c_ij)) | ((current_y == clin_y_ij) & (current_c == clin_c_ij))) & (current_age == clin_age_ij)
  }
  all_match[i] = all(all_match_i)
}

all(all_match)


