# -----------------------------------------------------------------------------
# Comparing Bayesian hierarchical spike-and-slab survival modeling to other 
# Bayesian survival models. 
# Comparison based on log posterior likelihood. 
# Author: Sarah Samorodnitsky
# -----------------------------------------------------------------------------

# Loading in the data:
currentwd <- "/Users/sarahsamorodnitsky/Documents/PanCancerOmics/HierarchicalSS_PanCanPanOmics/"
load(paste0(currentwd, "XYC_V2_WithAge_StandardizedPredictors.rda"))

# Loading in the helper functions
source(paste0(currentwd, "HelperFunctions.R"))

# -----------------------------------------------------------------------------
# Comparing to hierarchical horseshoe prior from Maity et al. (2020)
# -----------------------------------------------------------------------------

# Loading in the package with the model. 
library("PanCanVarSel")
