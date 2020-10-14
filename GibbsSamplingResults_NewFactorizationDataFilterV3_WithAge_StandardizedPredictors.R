# Running the hierarchical, spike-and-slab model (intercept outside S/S) on the
# new factorization data and including age which is standardized as are the other predictors. 
# Author: Sarah Samorodnitsky and Eric Lock (2020)
# University of Minnesota

# Loading in the important files
source("ExtendedHierarchicalModelLogNormalSpikeAndSlabInterceptOutsideSS.R")
source("HelperFunctions.R")
load("XYC_V2_WithAge_StandardizedPredictors.rda")

# Miscellaneous important values
covariates_by_cancer_full = lapply(Covariates, function(cancer) as.numeric(colnames(cancer))) # the covariates each cancer type has (each should have a 1 for the intercept)
covariates_table = table(unlist(covariates_by_cancer_full)) # gives the number of cancer types who have data on each PC
covariates_in_model_full = as.numeric(names(covariates_table))

covariates_in_model = covariates_in_model_full
p = length(covariates_in_model)
covariates_by_cancer = covariates_by_cancer_full # the name of this variable in the HelperFunctions.R file
iters = 1e5
n_cancer = length(Covariates)

cancer_types = names(Covariates)

# Creating the priors argument
priors = list(betatilde_priorvar_intercept = 10^2,    # \tilde\beta_0 ~ N(0,100). both beta tildes are centered at 0
              betatilde_priorvar_coefficient = 1,     # -> \tilde\beta_p ~ N(0,1), p > 0
              lambda2_priorshape_intercept = 1,       # -> \lambda^2_0 ~ IG(1,1)
              lambda2_priorrate_intercept = 1,         
              lambda2_priorshape_coefficient = 5,     # -> \lambda^2_p ~ IG(5,1) p > 0
              lambda2_priorrate_coefficient = 1,
              sigma2_priorshape = .01,                # -> \sigma^2 ~ IG(0.01, 0.01)
              sigma2_priorrate = .01,
              spike_priorvar = 1/10000)

##Generating starting values
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

# Gathering all the starting values together
starting_values = list(sigma2_start = sigma2_start,
                       beta_tilde_start = beta_tilde_start,
                       lambda2_start = lambda2_start,
                       gamma_start = gamma_start, # this will get passed on to be the list that stores posterior draws
                       pi_start = pi_start)

# Running the model
start <- Sys.time()
Posteriors_100k_WithAge_StandardizedPredictors = HierarchicalLogNormalSpikeSlab(Covariates, Survival, Censored, 
                                            starting_values,
                                            iters, covariates_in_model, 
                                            covariates_by_cancer, priors, 
                                            pi_generation = "shared_across_cancers")
end <- Sys.time()
end - start

save(Posteriors_100k_WithAge_StandardizedPredictors, file = "GibbsSamplingResults_100kiters_Chain1_WithAge_StandardizedPredictors.rda")


load(paste(datawd, "GibbsSamplingResults_100kiters_Chain1_WithAge_StandardizedPredictors.rda", sep = ""))

# Looking at the results for the model with age and standardized predictors
thinned_iters = seq(1,(iters/2), by = 10) # for thinning 
sigma2.gibbs = Posteriors_100k_WithAge_StandardizedPredictors$sigma2[(iters/2):iters][thinned_iters]
lambda2.gibbs = Posteriors_100k_WithAge_StandardizedPredictors$lambda2[(iters/2):iters,][thinned_iters,]
betatilde.gibbs = Posteriors_100k_WithAge_StandardizedPredictors$beta_tilde[(iters/2):iters,][thinned_iters,]
betas.gibbs = lapply(1:n_cancer, function(type) {
  do.call(rbind, Posteriors_100k_WithAge_StandardizedPredictors$betas[[type]])[(iters/2):iters,][thinned_iters,]
})
gamma.gibbs = lapply(1:n_cancer, function(type) {
  do.call(rbind, Posteriors_100k_WithAge_StandardizedPredictors$gamma[[type]])[(iters/2):iters,][thinned_iters,]
})
pi.gibbs = Posteriors_100k_WithAge_StandardizedPredictors$pi[(iters/2):iters,][thinned_iters, ]

# Plotting of the results

# Visualizing the inclusion indicators

library(RColorBrewer)
library(dplyr)
library(tidyr)
library(ggplot2)

##taking the mean of the gammas after burn-in for each model
inclusion.gibbs = matrix(nrow = length(cancer_types), ncol = p)
for (i in 1:length(cancer_types)) {
  cov = covariates_by_cancer[[i]]
  n_cov = length(cov)
  for (j in 1:n_cov) {
    k = which(covariates_in_model %in% cov[j])
    inclusion.gibbs[i, k] = mean(gamma.gibbs[[i]][,j])
  }
}

##adding column and row names to each inclusion matrix
colnames(inclusion.gibbs) = covariates_in_model
rownames(inclusion.gibbs) = cancer_types

##convert inclusion matrix into a tibble for ggplot
##make a heatmap of the estimated inclusion probability

inclusion.gibbs = as_tibble(inclusion.gibbs)
x.gibbs = gather(data = inclusion.gibbs, key = pc, value = inclusion, 1:p)
x.gibbs$cancer = rep(cancer_types, p) 
x.gibbs$pc = factor(x.gibbs$pc, levels = covariates_in_model)

pdf("GibbsInclusionHeatmapNewDataNewFiltering100kChain1_WithAge_StandardizedPredictors.pdf",
    width = 22, height = 22)
ggplot(data = x.gibbs, mapping = aes(x = pc, y = cancer, fill = inclusion)) +
  geom_tile() +
  xlab(label = "PC") +
  ggtitle("Inclusion Heatmap, Gibbs Sampler, Intercept Outside SS, Age Included, Predictors Standardized") + 
  geom_text(aes(label = round(inclusion,1))) 
dev.off()


########################################################################
############## Heatmap with changed x-axis labels ######################
########################################################################

# Removing the intercept to exclude it from heatmap
x.gibbs.without.intercept = x.gibbs[x.gibbs$pc != 0,]

pdf("GibbsInclusionHeatmapNewDataNewFiltering100kChain1_WithAge_StandardizedPredictors_FINAL.pdf",
    width = 40, height = 20)
ggplot(data = x.gibbs.without.intercept, mapping = aes(x = pc, y = cancer, fill = inclusion)) +
  geom_tile() +
  xlab(label = "PC") +
  ggtitle("Heatmap of Posterior Inclusion Probabilities by Cancer and Predictor") + 
  geom_text(aes(label = round(inclusion,1)), size=5) +
  scale_x_discrete(labels=c("0.5" = "Age")) +
  labs(y="Cancer Type", x = "Predictor") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0), size = 15),
        axis.title.x = element_text(margin = margin(t = 15, r = 20, b = 0, l = 0), size = 15),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 0.6),
        title = element_text(size = 16))
dev.off()
