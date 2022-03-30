# A Hierarchical Spike-and-Slab Model for Pan-Cancer Survival Using Pan-Omic Data

This is the repository for code used for the project titled "A Hierarchical Spike-and-Slab Model for Pan-Cancer Survival Using Pan-Omic Data." The purpose of this project was to demonstrate prediction using multiple sources of data and multiple sample sets. We used results of a method for bidimensional integration of multi-source, multi-way data called BIDIFAC+ by Lock et al. (https://arxiv.org/abs/2002.02601) to predict overall patient survival from the Cancer Genome Atlas database. Our predictive model is a Bayesian hierarchical model with hierarchical spike-and-slab predictors that allow for the borrowing of information across sample sets when determining the sparsity structure of the data. This repository also provides simulation results for this hierarchical model and others under different data-generating conditions. 

The following are instructions on how to reproduce our results. 

## Performing the BIDIFAC+ Analysis

To reproduce the modules derived from the TCGA data using BIDIFAC+, `TCGA_BIDIFAC_plus.R` loads in the BIDIFAC+ code from `bidifac_plus.R.` At the end of the `TCGA_BIDIFAC_plus.R` script, the modules are saved. These are then loaded in through the `load_modules.R` script in the next step. Due to the large file size of the TCGA data, we do not store these data here. Please email us at (samor007@umn.edu \& elock@umn.edu) if you would like access to those files.

## Reading in and processing the BIDIFAC+ results

`load_modules.R` contains code for loading in the BIDIFAC+ modules and computing their SVDs. The scores that we selected based on our filtering criteria will be saved at the end. This script will save a file called `mod.pcs.v3.rda` which contains the selected features from the BIDIFAC+ results. These features are in a variable named `mod.pcs`. Due to the large size of the BIDIFAC+ modules and the scores resulting from our filtering criteria, we don't share them here. Please email us at (samor007@umn.edu \& elock@umn.edu) if you would like access to those files. 

Then, `MatchingClinicalAndFactorizationData.R` matches the features of the BIDIFAC+ module SVDs to TCGA clinical data. You will have to change all `load()` statements to reflect your directory. In the end, this script will save a file called `XYC_V2_WithAge_StandardizedPredictors.rda` which contains the features, the survival data, and the censored data for all the subjects. This file is provided in this repository. We also provide a file called `XYC_V2_WithAge.rda` which contains the same data as `XYC_V2_WithAge_StandardizedPredictors.rda` only with non-standardized predictors. `XYC_V2_WithAge.rda` is only used for the next step of our analysis, which is to check that the data matches properly. 

Next, `CheckingDataMatchingResults.R` checks proper matching of subjects from modules to TCGA clinical data. This ensures data is of proper form (no missing, positive survival times, etc). This script begins by deconstructing the matching function from the previous script and then performs an explicit checker to ensure the data matches properly. In order to properly check that the age column matches, we manually change the following line out of the `MatchingClinicalAndFactorizationData.R`script from:

```
age.centered.i = (X_i[, "0.5"] - mean(X_i[, "0.5"]))/sd(X_i[, "0.5"]
```

to

```
age.centered.i = X_i[, "0.5"]
```

because the clinical data has the original age values, not the standardized version. We saved this version of the data in a file called `XYC_V2_WithAge.rda` to match our naming in `CheckingDataMatchingResults.R`. `XYC_V2_WithAge.rda` is only used for checking the data was cleaned properly and will not be used for the rest of the analysis. 

Again, we do not provide the data to run these last three scripts because it is too big to store locally. Please email me if you'd like the data we used and it can be sent through another avenue. We save both `XYC_V2_WithAge.rda` and `XYC_V2_WithAge_StandardizedPredictors.rda` to this repository so you can begin your analysis here. 

The last two scripts and most of the later ones source `HelperFunctions_V2.R` which contains helper functions used throughout project and analysis.

## The Gibbs sampler 

`ExtendedHierarchicalModelLogNormalSpikeAndSlabInterceptOutsideSS.R` contains code for the Gibbs sampling algorithm used in analysis. Assumes *log-normally* distributed survival and applies hierarchical spike-and-slab priors to coefficients, excluding the intercept. We use a version of this model that assumes *normally* distributed survival to check proper coverage for simplicitly. Code for the normal-likelihood model can be found in `ExtendedHierarchicalModelNormalSpikeAndSlab.R`. This model is only used for coverage checking purposes. 

The coverage simulation can be found in `ExtendedHierarchicalModelSpikeAndSlabCoverageSimulation.R` which ensures 95% coverage of credible intervals generated by Gibbs sampling algorithm. `ExtendedHierarchicalModelSpikeAndSlabSelectingTrueParametersSimulation.R` contains code for the selection accuracy simulation to ensure Gibbs sampling algorithm correctly identifies predictors to include in model. This simulation should give nominal accuracy, around 90%. 

## TCGA data application

`ModelComparison_V2.R` compares eight different possible models for best fit on TCGA data. We conclude the hierarchical spike-and-slab model we propose fits the data best so we apply it to the TCGA data in `GibbsSamplingResults_NewFactorizationDataFilterV3_WithAge_StandardizedPredictors.R`. This script runs the Gibbs sampling algorithm on the TCGA data for the log-normal model. Results are saved for further analysis. This script also creates inclusion heatmap we include in our article. `ScoresVsSubtypes.R` contains code to create all additional figures in the article. To create these figures, we used data from the TCGA-CDR (Liu and others, 2018) found in `TCGA-CDR.csv`, data from an analysis by TCGA Research Network (2015) found in `LGGClinicalDataSubtypes.csv`, and data from the analysis done in Ricketts and others (2018) found in `mmc2.csv`. The links to these three articles are found at https://pubmed.ncbi.nlm.nih.gov/29625055/, https://www.nejm.org/doi/full/10.1056/nejmoa1402121, and https://pubmed.ncbi.nlm.nih.gov/29617669/, respectively. 

## Model simulations

We run our model simulations to compare the performance of our model with seven other Bayesian hierarchical frameworks in `ModelSimulations_V2.R`. We compare these models under six different data-generating scenarios to assess the performance and flexibility of each. This script sources `RunFiveModelTypes_V2.R` which runs the eight different possible models. At the end of the simulations script, pairwise t-tests are performed to compare the results from each model under each condition. 

