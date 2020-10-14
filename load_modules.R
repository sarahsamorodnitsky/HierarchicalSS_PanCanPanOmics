# Code for loading in factorization modules
# Authors: Sarah Samorodnitsky and Eric Lock (2020)
# University of Minnesota

library(foreach)
library(doParallel)

##load factorization results
load('~/PanTCGA/pan.fac.results_v2.rda') # second dataset Eric sent

#data sources in each module (modules are ordered by variance explained)
source.list

#tumor types in each module
type.list

#'S.list' gives a list of the low-rank matrix for each module (for the rows and columns that it is present for)
#Module 1:
S.list[[1]]
dim(S.list[[1]])

#sample.id.list is a list of sample ids present for each module i, corresponding to the columns of S[[i]]
sample.id.list[[1]]
length(sample.id.list[[1]]) # the length matches the number of columns in the low-rank matrix

#tum.type.list is a list of tumor types for each module i, corresponding to the columns of S[[i]]
tum.type.list[[1]]

##Compute the first principal component for module 1:
SVD.1 <- svd(S.list[[1]],nu=2,nv=2) #this can take a couple minutes - the SVD for smaller modules is faster 

# first principal component scores for module 1 (to be used for prediction)
PC1 <- SVD.1$d[1]*SVD.1$v[,1] # length matches the number of patients
PC2 <- SVD.1$d[2]*SVD.1$v[,2]

##Compute the first principal component for module 1:
SVD.2 <- svd(S.list[[2]],nu=2,nv=2) #this can take a couple minutes - the SVD for smaller modules is faster 

# first principal component scores for module 1 (to be used for prediction)
PC3 <- SVD.2$d[1]*SVD.2$v[,1] # length matches the number of patients
PC4 <- SVD.2$d[2]*SVD.2$v[,2]

mod.pcs = list(cbind(PC1,PC2), cbind(PC3,PC4))

# Checking the dimensions of each data type match
n_mod = length(mod.pcs)
SListTumSampleLengthMatch = c()
for (i in 1:n_mod) {
  SListTumSampleLengthMatch[i] = (ncol(S.list[[i]]) == length(tum.type.list[[i]])) & (length(tum.type.list[[i]]) == length(sample.id.list[[i]]))
}
all(SListTumSampleLengthMatch) # the given data types all match in dimension

# Running the SVDs for each module
# Parallelizing to speed things up
cl <- makeCluster(6) 
registerDoParallel(cl)
mod.SVD = foreach (mod=1:length(source.list)) %dopar% {
  svd(S.list[[mod]]) 
  # as.matrix(SVD$d[1]*SVD$v[,1], ncol = 1)
}
stopCluster(cl)

save(mod.SVD, file = "mod.SVD.v2.rda")

# Looking at the scale of the SVD principal components
load("mod.SVD.v2.rda") # loads mod.SVD

# First looking at the eigenvalues alone
mod.eigenvalues = lapply(1:length(mod.SVD), function(mod) {
  mod.SVD[[mod]]$d^2
}) 

### 6/18/20 Filtering method: 
###     1. Include first PC from each module
###     2. Consider only the singular values that are greater than 0
###         2.1 Of these, select PCs that have eigenvalue/variance of pre-fac data ratio > 0.01

load("pan.fac.data_v2.rda") # loading in the pre-factorization data 

MaxNumSingValsNonZero = max(sapply(sapply(mod.eigenvalues, function(mod) which(round(mod, 5) >0)), max)) # == 37, take just the first 50 values

svd.ratios = lapply(1:length(mod.SVD), function(i) list())
for (j in 1:length(mod.SVD)) { # for the jth module
  SVD.j = mod.SVD[[j]]
  modulej.ratios = c()
  for (i in 1:50) { # for each component in module j, consider only the first 50 singular values
    modulej.ratios[i] = SVD.j$d[i]^2/sum(X0[p.ind.list.order[[j]],n.ind.list.order[[j]]]^2) # data is normalized to have mean 0
  }
  svd.ratios[[j]] = modulej.ratios
}

svd.ratios

save(svd.ratios, file = "svd.ratios.rda")

load("svd.ratios.rda")

mod.pcs = lapply(1:length(mod.SVD), function(mod) list()) # create list to store principal components from each module
for (i in 1:length(mod.SVD)) {
  svd.ratios.i = svd.ratios[[i]] # ratios for ith module
  ratios_g0.01 = svd.ratios.i >= 0.01 # select the ratios > 0.01
  ratios_g0.01[1] = TRUE # ensure we always pick the first PC 
  inds_g0.01 = which(ratios_g0.01)
  mod.pcs.i = as.matrix(mod.SVD[[i]]$d[inds_g0.01]*mod.SVD[[i]]$v[, inds_g0.01])
  colnames(mod.pcs.i) = paste(i, ".", inds_g0.01, sep = "") # module number "." PC number
  mod.pcs[[i]] = mod.pcs.i
}

# Using this filtering approach, I selected the same PCs as before. 

save(mod.pcs, file = "mod.pcs.v3.rda")
