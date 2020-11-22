# This script explores the relationship between components selected by our model and the cancer subtypes. 
# We will see if selected components appear to explain subtype differences.

# Save the important working directories
currentwd = "~/PanTCGA/SpikeAndSlabResults/SpikeAndSlabGibbsSampler/ModelWithInterceptOutsideSS/ResultsAfterFiltering_V3/ScoreVsSubtype/"
PanTCGAwd = "~/PanTCGA/"
datawd = "~/PanTCGA/SpikeAndSlabResults/SpikeAndSlabGibbsSampler/ModelWithInterceptOutsideSS/ResultsAfterFiltering_V3/"

# load in the Covariates, Survival, Censored data to generate fake data with the same structure
load(paste(datawd, "XYC_V2_WithAge_StandardizedPredictors.rda", sep = ""))

# load in helper functions
source(paste(PanTCGAwd, "HelperFunctions.R", sep = ""))

# The cancers we wish to explore
UCEC.Covariates = Covariates$UCEC
LGG.Covariates = Covariates$LGG
Kidney.Covariates = Covariates[c("KICH", "KIRC", "KIRP")]

# Using the TCGA-CDR data
clinical_data <- read.csv(paste(PanTCGAwd, "TCGA-CDR.csv", sep = ""))

#################################################################################
##################### First: UCEC and component 16.1   ##########################
#################################################################################

# Ensure clinical data contains the same patient barcodes
# more patients were used in the analysis than are in the clinical data
# using TCGA-CDR data.
ucec.clinical <- clinical_data[clinical_data$type == "UCEC", ]

# Extracting the names available in both datasets
ucec.patient.ids.analysis = rownames(UCEC.Covariates) 
# Using the clinical data resource.
ucec.patient.ids.clinical <- ucec.clinical$bcr_patient_barcode

# More clinical ids that analysis ids
length(ucec.patient.ids.analysis) < length(ucec.patient.ids.clinical)

# The ids we will use to assess the relationship between subtype and score
UCEC.Covariates.Avail = UCEC.Covariates # no need to subset since the ids are already a subset of the TCGA CDR
ucec.patient.ids.avail = ucec.patient.ids.analysis # the ones used in the analysis are the available ids.

# Match histological subtype and 16.1 together using patient barcode. 
# Just the patient ids used in the analysis
ucec.clin.data.available = ucec.clinical[ucec.patient.ids.clinical %in% ucec.patient.ids.avail, ]

# Matching the order to the data used in the analysis
ucec.clin.data.available.reorder = ucec.clin.data.available[match(ucec.patient.ids.avail, ucec.clin.data.available$bcr_patient_barcode), ]

# Check matching
nrow(ucec.clin.data.available.reorder) == nrow(UCEC.Covariates.Avail) # should be true
all(ucec.clin.data.available.reorder$bcr_patient_barcode == ucec.patient.ids.avail) # should be true and the avail ids are in the same order as the analysis ids because they came from the analysis dataset
all(ucec.clin.data.available.reorder$bcr_patient_barcode == rownames(UCEC.Covariates.Avail)) # should be true and the avail ids are in the same order as the analysis ids because they came from the analysis dataset

# combine the histological data from the clinical dataset with the covariates data for UCEC
UCEC.Covariates.Avail.Histological = cbind.data.frame(UCEC.Covariates.Avail, histology = ucec.clin.data.available.reorder$histological_type)
UCEC.Covariates.Avail.Histological$histology <- ucec.clinical <- droplevels(UCEC.Covariates.Avail.Histological$histology)

# Plotting the scores with the histological type
# The three histological types are Endometrioid, Mixed, and Serous
library(viridis)

par_col <- viridis::viridis(3) # colors to use
par_col[par_col == "#FDE725FF"] <- "#A87C07" # Changing the yellow to an ochre - easier to see. 
par_col_ucec <- par_col[UCEC.Covariates.Avail.Histological$histology] # selecting colors

# Checking that the color matches the in the legend to the plot
# with this website: https://htmlcolorcodes.com/ copy and paste the 
# color value output and make sure it matches with the plot
unique(par_col_ucec[UCEC.Covariates.Avail.Histological$histology == "Endometrioid endometrial adenocarcinoma"]) # "#440154FF" == purple
unique(par_col_ucec[UCEC.Covariates.Avail.Histological$histology == "Serous endometrial adenocarcinoma"]) # changed to "#A87C07" which is orange. "#FDE725FF" == yellow
unique(par_col_ucec[UCEC.Covariates.Avail.Histological$histology == "Mixed serous and endometrioid"]) # "#21908CFF" == teal


###############################################################################
#################################  UCEC KDE Plots #############################
###############################################################################

# Library for plotting kernel densities
library(lattice)

# density plot for PC 16.1 stratified by cancer subtype
png(paste(currentwd, "UCEC_16.1_KDEPlot_TCGA_CDR.png", sep = ""))
densityplot(~c(UCEC.Covariates.Avail.Histological[, "16.1"]), 
                 groups = UCEC.Covariates.Avail.Histological[, "histology"], 
                 xlab = "Component 16.1",
                 main = "Component 16.1 Scores for UCEC \n by Histological Subtype",
                 par.settings = list(superpose.line = list(col=par_col)),
                 auto.key = list(rows = 3),
                 plot.points = FALSE)
dev.off()

# density plot for PC 16.1 stratified by cancer subtype
png(paste(currentwd, "UCEC_5.1_KDEPlot_TCGA_CDR.png", sep = ""))
densityplot(~c(UCEC.Covariates.Avail.Histological[, "5.1"]), 
            groups = UCEC.Covariates.Avail.Histological[, "histology"], 
            xlab = "Component 5.1",
            main = "Component 5.1 Scores for UCEC \n by Histological Subtype",
            par.settings = list(superpose.line = list(col=par_col)),
            auto.key = list(rows = 3),
            plot.points = FALSE)
dev.off()

# Check the legend matches with the automatically generated one
# densityplot(~c(UCEC.Covariates.Avail.Histological[, "16.1"]),
#             groups = UCEC.Covariates.Avail.Histological[, "histology"],
#             xlab = "Scores for Module 16.1",
#             main = "UCEC Scores by Subtype",
#             auto.key = TRUE)

########################################################################
################ UCEC Scatterplot ######################################
########################################################################

png(paste(currentwd, "UCEC_Scatterplot_TCGA_CDR.png", sep = ""))
plot(UCEC.Covariates.Avail.Histological$`5.1`, UCEC.Covariates.Avail.Histological$`16.1`, 
     col = par_col_ucec,
     xlab = "Component 5.1", ylab = "Component 16.1",
     main = "Sample Scores for UCEC, \n Labeled by Histological Type")
legend("bottomleft", legend = levels(UCEC.Covariates.Avail.Histological$histology),
       col=par_col,pch=1,
       cex = 0.75)
dev.off()

###############################################################################
##############################  UCEC Survival Plots ###########################
###############################################################################

library(survival)

# Storing the survival time for UCEC
UCEC.Survival = Survival$UCEC
UCEC.Censored = Censored$UCEC

# Using the CDR data
UCEC.Survival.Avail = UCEC.Survival
UCEC.Censored.Avail = UCEC.Censored

# Boolean vector for censored or non-censored (TRUE = non-censored, had the event; FALSE = censored)
UCEC.EventIndicator = !is.na(UCEC.Survival.Avail)

# Filling in the censored observations with their censor times
UCEC.Survival.Avail[is.na(UCEC.Survival.Avail)] = UCEC.Censored.Avail[is.na(UCEC.Survival.Avail)]

# Making sure the IDs match between the analysis and clinical data
all(names(UCEC.Survival.Avail) == ucec.clin.data.available.reorder$bcr_patient_barcode)

# Using the TCGA CDR
Survival.Histology.Avail = cbind.data.frame(survival = UCEC.Survival.Avail, 
                                            histology = ucec.clin.data.available.reorder$histological_type,
                                            event = UCEC.EventIndicator)
Survival.Histology.Avail$histology <- droplevels(Survival.Histology.Avail$histology)

# Making sure the survival times still match with the survival times from before
all(Survival.Histology.Avail[,"survival"] == UCEC.Survival.Avail)

# PLOT
png(paste(currentwd, "UCEC_Subtype_SurvivalCurves_TCGA_CDR.png", sep = ""))
plot(survfit(Surv(survival,event, type = "right") ~ histology, 
             data = Survival.Histology.Avail),
     xlab = "Days",
     ylab = "Overall Survival Probability",
     col = par_col,
     main = "Kaplan-Meier Curves for UCEC Patients \n by Histological Subtype")
legend("bottomleft", legend = levels(Survival.Histology.Avail$histology), 
       col=par_col, lty=1,
       cex = 0.75)

dev.off()

#########################################################################################
##################### Creating supplementary figures for UCEC 5.1 #######################
#########################################################################################

########################################################
### Dividing subjects into groups based on 5.1 score ###
########################################################
all(names(UCEC.Survival.Avail) == rownames(UCEC.Covariates.Avail)) # checking order

# Splitting subjects into groups depending on whether 5.1 is greater or less than 0
UCEC.5.1.Groups <- UCEC.Covariates.Avail[,"5.1"] >= 0
UCEC.5.1.Groups <- factor(ifelse(UCEC.5.1.Groups, "Greater Than 0", "Less Than 0"))

# Putting the data together
Survival.5.1Groups.Avail = cbind.data.frame(survival = UCEC.Survival.Avail, 
                                            group = UCEC.5.1.Groups,
                                            event = UCEC.EventIndicator)

# Making a survival plot
png(paste(currentwd, "UCEC_5.1Groups_Survival.png", sep = ""))
plot(survfit(Surv(survival,event, type = "right") ~ group, 
             data = Survival.5.1Groups.Avail),
     xlab = "Days",
     ylab = "Overall Survival Probability",
     col = par_col,
     main = "Kaplan-Meier Curves for UCEC Patients \n by Component 5.1 Scores")
legend("bottomleft", legend = levels(Survival.5.1Groups.Avail$group), 
       col=par_col, lty=1,
       cex = 0.75)
dev.off()

################################
### Didn't end up using these ##
################################
# Looking at histological grade
UCEC.Covariates.Avail.Histological.Grade <- cbind.data.frame(UCEC.Covariates.Avail.Histological, histological.grade = ucec.clin.data.available.reorder$histological_grade)
UCEC.Covariates.Avail.Histological.Grade$histological.grade <- droplevels(UCEC.Covariates.Avail.Histological.Grade$histological.grade)

# Creating the figure
par_col.2 <- viridis::viridis(nlevels(UCEC.Covariates.Avail.Histological.Grade$histological.grade)) # colors to use
par_col.2[par_col.2 == "#FDE725FF"] <- "#A87C07" 
par_col_ucec.2 <- par_col.2[UCEC.Covariates.Avail.Histological.Grade$histological.grade] # selecting colors

densityplot(~c(UCEC.Covariates.Avail.Histological.Grade[, "5.1"]), 
            groups = UCEC.Covariates.Avail.Histological.Grade[, "histological.grade"], 
            xlab = "Component 5.1",
            main = "UCEC Scores by Subtype",
            par.settings = list(superpose.line = list(col=par_col.2)),
            auto.key = list(rows = 3),
            plot.points = FALSE)

# The intersection between grade and histological type
UCEC.Covariates.Avail.Histological.Grade$histxgrade <- as.factor(paste(UCEC.Covariates.Avail.Histological.Grade$histology,
                                                                       " (",
                                                                       UCEC.Covariates.Avail.Histological.Grade$histological.grade,
                                                                       ")", sep = ""))
# Setting the colors
par_col_histxgrade <- viridis(nlevels(Survival.Histology.Avail$histxgrade))
par_col_histxgrade[par_col_histxgrade == "#FDE725FF"] <- "#A87C07" 

png(paste(currentwd, "UCEC_5.1_KDE_Grade_TCGA_CDR.png", sep = ""))
densityplot(~c(UCEC.Covariates.Avail.Histological.Grade[, "5.1"]), 
            groups = UCEC.Covariates.Avail.Histological.Grade[, "histxgrade"], 
            xlab = "Component 5.1",
            main = "UCEC Scores by Histological Grade",
            par.settings = list(superpose.line = list(col=par_col_histxgrade)),
            auto.key = list(rows = 8),
            plot.points = FALSE)
dev.off()

# Making a survival figure
Survival.Histology.Avail = cbind.data.frame(Survival.Histology.Avail,
                                            histological.grade = ucec.clin.data.available.reorder$histological_grade,
                                            histxgrade = as.factor(paste(ucec.clin.data.available.reorder$histological_type,
                                                                         " (",
                                                                         ucec.clin.data.available.reorder$histological_grade,
                                                                         ")", sep = "")))

png(paste(currentwd, "UCEC_Grade_SurvivalCurves_TCGA_CDR.png", sep = ""),
    width = 700, height = 700, units = "px", pointsize = 16)
plot(survfit(Surv(survival,event, type = "right") ~ histxgrade, 
             data = Survival.Histology.Avail),
     xlab = "Days",
     ylab = "Overall Survival Probability",
     col = par_col_histxgrade,
     main = "Kaplan-Meier Curves for UCEC Patients \n by Histological Grade")
legend("bottomleft", legend = levels(Survival.Histology.Avail$histxgrade), 
       col=par_col_histxgrade, lty=1,
       cex = 0.6)
dev.off()

# Just checking the survival curves match
# Survival.Histology.Avail.Endo = cbind.data.frame(survival = UCEC.Survival.Avail[ucec.clin.data.available.reorder$histological_type == "Endometrioid endometrial adenocarcinoma"],
#                                             histology = ucec.clin.data.available.reorder$histological_type[ucec.clin.data.available.reorder$histological_type == "Endometrioid endometrial adenocarcinoma"],
#                                             event = UCEC.EventIndicator[ucec.clin.data.available.reorder$histological_type == "Endometrioid endometrial adenocarcinoma"])
# 
# plot(survfit(Surv(survival,event, type = "right") ~ histology,
#              data = Survival.Histology.Avail.Endo),
#      xlab = "Days",
#      ylab = "Overall Survival Probability",
#      main = "Kaplan-Meier Curves for UCEC Patients by Histological Subtype")
# 
# Survival.Histology.Avail.Serous = cbind.data.frame(survival = UCEC.Survival.Avail[ucec.clin.data.available.reorder$histological_type == "Serous endometrial adenocarcinoma"],
#                                                  histology = ucec.clin.data.available.reorder$histological_type[ucec.clin.data.available.reorder$histological_type == "Serous endometrial adenocarcinoma"],
#                                                  event = UCEC.EventIndicator[ucec.clin.data.available.reorder$histological_type == "Serous endometrial adenocarcinoma"])
# 
# plot(survfit(Surv(survival,event, type = "right") ~ histology,
#              data = Survival.Histology.Avail.Serous),
#      xlab = "Days",
#      ylab = "Overall Survival Probability",
#      main = "Kaplan-Meier Curves for UCEC Patients by Histological Subtype")
# 
# Survival.Histology.Avail.Mixed = cbind.data.frame(survival = UCEC.Survival.Avail[ucec.clin.data.available.reorder$histological_type == "Mixed serous and endometrioid"],
#                                                    histology = ucec.clin.data.available.reorder$histological_type[ucec.clin.data.available.reorder$histological_type == "Mixed serous and endometrioid"],
#                                                    event = UCEC.EventIndicator[ucec.clin.data.available.reorder$histological_type == "Mixed serous and endometrioid"])
# 
# plot(survfit(Surv(survival,event, type = "right") ~ histology,
#              data = Survival.Histology.Avail.Mixed),
#      xlab = "Days",
#      ylab = "Overall Survival Probability",
#      main = "Kaplan-Meier Curves for UCEC Patients by Histological Subtype")


#################################################################################
##################### Next: LGG and component 7.2, 12.3   #######################
#################################################################################

# Load in clinical data from TCGA (2015) 
lgg_subtypes = read.csv(paste(currentwd, "LGGClinicalDataSubtypes.csv", sep = ""))

# Ensure clinical data contains the same patient barcodes
length(rownames(LGG.Covariates)) > nrow(lgg_subtypes) # more patients were used in the analysis than are in the clinical data, just like before

# extracting the names available in both datasets
lgg.patient.ids.analysis = rownames(LGG.Covariates) 
lgg.patient.ids.clinical = lgg_subtypes$Tumor

# Are any in the clinical data? If so, what percentage is available?
any(lgg.patient.ids.clinical %in% lgg.patient.ids.analysis)
sum(lgg.patient.ids.clinical %in% lgg.patient.ids.analysis)/nrow(lgg_subtypes)

# the ids we will use to assess the relationship between subtype and score
LGG.Covariates.Avail = LGG.Covariates[lgg.patient.ids.analysis %in% lgg.patient.ids.clinical, ]
lgg.patient.ids.avail = lgg.patient.ids.analysis[lgg.patient.ids.analysis %in% lgg.patient.ids.clinical]

# Match histological subtype and 7.2 together using patient barcode. 
# Just the patient ids used in the analysis
lgg.clin.data.available = lgg_subtypes[lgg_subtypes$Tumor %in% lgg.patient.ids.avail, ]

# Matching the order to the data used in the analysis
lgg.clin.data.available.reorder = lgg.clin.data.available[match(lgg.patient.ids.avail, lgg.clin.data.available$Tumor), ]

# Check matching
nrow(lgg.clin.data.available.reorder) == nrow(LGG.Covariates.Avail) # should be true
all(lgg.clin.data.available.reorder$Tumor == lgg.patient.ids.avail) # should be true

# combine the histological data from the clinical dataset with the covariates data for UCEC
LGG.Covariates.Avail.IDH = cbind.data.frame(LGG.Covariates.Avail, 
                                            idh_subtype = lgg.clin.data.available.reorder$IDH.1p19q.Subtype)

nrow(LGG.Covariates.Avail.IDH) # checking how many there were before removing NAs
sum(is.na(LGG.Covariates.Avail.IDH$idh_subtype)) # checking how many NAs there are

# Removing the subjects who are missing an IDH label.
LGG.Covariates.Avail.IDH <- LGG.Covariates.Avail.IDH[!is.na(LGG.Covariates.Avail.IDH$idh_subtype),]

# Checking how many subjects there are now to make sure the proper number were removed
nrow(LGG.Covariates.Avail.IDH) ## 263 == 274 - 11

# Just checking
any(is.na(LGG.Covariates.Avail.IDH$idh_subtype)) # should be false

# Storing a separate LGG.Covariates.Avail.IDH to make sure that the relabeling was done correctly
LGG.Covariates.Avail.IDH.test <- LGG.Covariates.Avail.IDH

# Relabeling IDH mutation groups to have more interpretable labels.
levels(LGG.Covariates.Avail.IDH$idh_subtype)[levels(LGG.Covariates.Avail.IDH$idh_subtype) == "IDHwt"] <- "Wildtype IDH"
levels(LGG.Covariates.Avail.IDH$idh_subtype)[levels(LGG.Covariates.Avail.IDH$idh_subtype) == "IDHmut-codel"] <- "IDH Mutation, 1p/19q Codeletion"
levels(LGG.Covariates.Avail.IDH$idh_subtype)[levels(LGG.Covariates.Avail.IDH$idh_subtype) == "IDHmut-non-codel"] <- "IDH Mutation, No 1p/19q Codeletion"

# Checking the relabeling was done correctly
all(which(LGG.Covariates.Avail.IDH$idh_subtype == "Wildtype IDH") == which(LGG.Covariates.Avail.IDH.test$idh_subtype == "IDHwt"))
all(which(LGG.Covariates.Avail.IDH$idh_subtype == "IDH Mutation, 1p/19q Codeletion") == which(LGG.Covariates.Avail.IDH.test$idh_subtype == "IDHmut-codel"))
all(which(LGG.Covariates.Avail.IDH$idh_subtype == "IDH Mutation, No 1p/19q Codeletion") == which(LGG.Covariates.Avail.IDH.test$idh_subtype == "IDHmut-non-codel"))

# Plotting the scores with the IDH mutation subtype
par_col_lgg <- par_col[LGG.Covariates.Avail.IDH$idh_subtype] # selecting colors

png(paste(currentwd, "LGG_SelectedComponents_Scatterplot.png", sep = ""))
plot(LGG.Covariates.Avail.IDH$`12.3`, LGG.Covariates.Avail.IDH$`7.2`, 
     col = par_col_lgg,
     xlab = "Component 12.3", ylab = "Component 7.2",
     main = "Sample Scores for LGG, \n Labeled by IDH Mutation Subtype")
legend("topleft", legend = levels(LGG.Covariates.Avail.IDH$idh_subtype), 
       col=par_col, pch=1,
       cex = 0.75)
dev.off()

# Checking that the colors match - checking visually here
# plot(LGG.Covariates.Avail.IDH$`12.3`[LGG.Covariates.Avail.IDH$idh_subtype == "IDHwt"], 
#      LGG.Covariates.Avail.IDH$`7.2`[LGG.Covariates.Avail.IDH$idh_subtype == "IDHwt"], 
#      xlab = "PC 12.3", ylab = "PC 7.2",
#      main = "LGG Selected Components, Labeled by IDH Mutation Subtype")
# 
# plot(LGG.Covariates.Avail.IDH$`12.3`[LGG.Covariates.Avail.IDH$idh_subtype == "IDHmut-non-codel"], 
#      LGG.Covariates.Avail.IDH$`7.2`[LGG.Covariates.Avail.IDH$idh_subtype == "IDHmut-non-codel"], 
#      xlab = "PC 12.3", ylab = "PC 7.2",
#      main = "LGG Selected Components, Labeled by IDH Mutation Subtype")
# 
# plot(LGG.Covariates.Avail.IDH$`12.3`[LGG.Covariates.Avail.IDH$idh_subtype == "IDHmut-codel"], 
#      LGG.Covariates.Avail.IDH$`7.2`[LGG.Covariates.Avail.IDH$idh_subtype == "IDHmut-codel"], 
#      xlab = "PC 12.3", ylab = "PC 7.2",
#      main = "LGG Selected Components, Labeled by IDH Mutation Subtype")

###############################################################################
#################################  LGG KDE Plots  #############################
###############################################################################

# density plot for PC 7.2 stratified by IDH mutation status
png(paste(currentwd, "LGG_7.2_KDE.png", sep = ""))
densityplot(~c(LGG.Covariates.Avail.IDH[, "7.2"]), 
            groups = LGG.Covariates.Avail.IDH[, "idh_subtype"], 
            xlab = "Component 7.2",
            main = "Component 7.2 Scores for LGG \n by IDH Mutation Subtype",
            par.settings = list(superpose.line = list(col=par_col)),
            auto.key = list(rows = 3),
            plot.points = FALSE)
dev.off()

# Compare with automatically generated colors to make sure the densities are correctly labeled
# densityplot(~c(LGG.Covariates.Avail.IDH[, "7.2"]),
#             groups = LGG.Covariates.Avail.IDH[, "idh_subtype"],
#             xlab = "Scores for Module 7.2",
#             main = "LGG Scores by Subtype",
#             auto.key = TRUE)

png(paste(currentwd, "LGG_12.3_KDE.png", sep = ""))
densityplot(~c(LGG.Covariates.Avail.IDH[, "12.3"]), 
            groups = LGG.Covariates.Avail.IDH[, "idh_subtype"], 
            xlab = "Component 12.3",
            main = "Component 12.3 Scores for LGG \n by IDH Mutation Subtype",
            par.settings = list(superpose.line = list(col=par_col)),
            auto.key = list(rows = 3),
            plot.points = FALSE)
dev.off()

# Compare with automatically generated colors to make sure they are correct
# densityplot(~c(LGG.Covariates.Avail.IDH[, "12.3"]),
#             groups = LGG.Covariates.Avail.IDH[, "idh_subtype"],
#             xlab = "Scores for Module 7.2",
#             main = "LGG Scores by Subtype",
#             auto.key = TRUE)

###############################################################################
##############################  LGG Survival Plots ############################
###############################################################################

library(survival)

# Storing the survival time for UCEC
LGG.Survival = Survival$LGG
LGG.Censored = Censored$LGG

# Selecting just the subjects who were available in both the analysis data and the clinical data
LGG.Survival.Avail = LGG.Survival[lgg.patient.ids.analysis %in% lgg.patient.ids.clinical] # all(names(LGG.Survival) == LGG.patient.ids.analysis) so can use either names or this vector
LGG.Censored.Avail = LGG.Censored[lgg.patient.ids.analysis %in% lgg.patient.ids.clinical] 

# Boolean vector for censored or non-censored (TRUE = non-censored, had the event; FALSE = censored)
LGG.EventIndicator = !is.na(LGG.Survival.Avail)

# Filling in the censored observations with their censor times
LGG.Survival.Avail[is.na(LGG.Survival.Avail)] = LGG.Censored.Avail[is.na(LGG.Survival.Avail)]

# Making sure the IDs match between the analysis and clinical data
all(names(LGG.Survival.Avail) == lgg.clin.data.available.reorder$bcr_patient_barcode)

# Adding the histology subtype
LGG.Survival.Histology.Avail = cbind.data.frame(survival = LGG.Survival.Avail, 
                                            idh_subtype = lgg.clin.data.available.reorder$IDH.1p19q.Subtype,
                                            event = LGG.EventIndicator)

# Making sure the survival times still match with the survival times from before
all(LGG.Survival.Histology.Avail[,"survival"] == LGG.Survival.Avail)

# Removing the NAs and changing the group labels
LGG.Survival.Histology.Avail <- LGG.Survival.Histology.Avail[!is.na(LGG.Survival.Histology.Avail$idh_subtype),]
levels(LGG.Survival.Histology.Avail$idh_subtype)[levels(LGG.Survival.Histology.Avail$idh_subtype) == "IDHwt"] <- "Wildtype IDH"
levels(LGG.Survival.Histology.Avail$idh_subtype)[levels(LGG.Survival.Histology.Avail$idh_subtype) == "IDHmut-codel"] <- "IDH Mutation, 1p/19q Codeletion"
levels(LGG.Survival.Histology.Avail$idh_subtype)[levels(LGG.Survival.Histology.Avail$idh_subtype) == "IDHmut-non-codel"] <- "IDH Mutation, No 1p/19q Codeletion"

# PLOT
png(paste(currentwd, "LGG_SurvivalCurves_Subtype.png", sep = ""))
plot(survfit(Surv(survival,event, type = "right") ~ idh_subtype, 
             data = LGG.Survival.Histology.Avail),
     xlab = "Days",
     ylab = "Overall Survival Probability",
     col = par_col,
     main = "Kaplan-Meier Curves for LGG Patients \n by IDH Mutation Subtype")
legend("bottomleft", legend = levels(LGG.Survival.Histology.Avail$idh_subtype), 
       col= par_col, lty=1,
       cex = 0.75)
dev.off()

# Just checking the survival curves match
# LGG.Histology.Avail.codel = cbind.data.frame(survival = LGG.Survival.Avail[lgg.clin.data.available.reorder$IDH.1p19q.Subtype == "IDHmut-codel"],
#                                             histology = lgg.clin.data.available.reorder$IDH.1p19q.Subtype[lgg.clin.data.available.reorder$IDH.1p19q.Subtype == "IDHmut-codel"],
#                                             event = LGG.EventIndicator[lgg.clin.data.available.reorder$IDH.1p19q.Subtype == "IDHmut-codel"])
# 
# plot(survfit(Surv(survival,event, type = "right") ~ histology,
#              data = LGG.Histology.Avail.codel),
#      xlab = "Days",
#      ylab = "Overall Survival Probability",
#      main = "Kaplan-Meier Curves for UCEC Patients by Histological Subtype")
# 
# LGG.Histology.Avail.noncodel = cbind.data.frame(survival = LGG.Survival.Avail[lgg.clin.data.available.reorder$IDH.1p19q.Subtype == "IDHmut-non-codel"],
#                                              histology = lgg.clin.data.available.reorder$IDH.1p19q.Subtype[lgg.clin.data.available.reorder$IDH.1p19q.Subtype == "IDHmut-non-codel"],
#                                              event = LGG.EventIndicator[lgg.clin.data.available.reorder$IDH.1p19q.Subtype == "IDHmut-non-codel"])
# 
# plot(survfit(Surv(survival,event, type = "right") ~ histology,
#              data = LGG.Histology.Avail.noncodel),
#      xlab = "Days",
#      ylab = "Overall Survival Probability",
#      main = "Kaplan-Meier Curves for UCEC Patients by Histological Subtype")
# 
# LGG.Histology.Avail.wt = cbind.data.frame(survival = LGG.Survival.Avail[lgg.clin.data.available.reorder$IDH.1p19q.Subtype == "IDHwt"],
#                                                 histology = lgg.clin.data.available.reorder$IDH.1p19q.Subtype[lgg.clin.data.available.reorder$IDH.1p19q.Subtype == "IDHwt"],
#                                                 event = LGG.EventIndicator[lgg.clin.data.available.reorder$IDH.1p19q.Subtype == "IDHwt"])
# 
# plot(survfit(Surv(survival,event, type = "right") ~ histology,
#              data = LGG.Histology.Avail.wt),
#      xlab = "Days",
#      ylab = "Overall Survival Probability",
#      main = "Kaplan-Meier Curves for UCEC Patients by Histological Subtype")

#################################################################################
####################### Supplementary LGG plots #################################
#################################################################################

#############################################################
### Dividing subjects up into groups based on 12.3 scores ###
#############################################################

# Splitting subjects into groups depending on whether 5.1 is greater or less than 0
LGG.12.3.Groups <- LGG.Covariates.Avail[,"12.3"] >= 0
LGG.12.3.Groups <- factor(ifelse(LGG.12.3.Groups, "Greater Than 0", "Less Than 0"))

# Putting the data together
Survival.12.3Groups.Avail = cbind.data.frame(survival = LGG.Survival.Avail, 
                                             group = LGG.12.3.Groups,
                                             event = LGG.EventIndicator)

# Making a survival plot
png(paste(currentwd, "LGG_12.3Groups_Survival.png", sep = ""))
plot(survfit(Surv(survival,event, type = "right") ~ group, 
             data = Survival.12.3Groups.Avail),
     xlab = "Days",
     ylab = "Overall Survival Probability",
     col = par_col,
     main = "Kaplan-Meier Curves for LGG Patients \n by Component 12.3 Scores")
legend("bottomleft", legend = levels(Survival.12.3Groups.Avail$group), 
       col=par_col, lty=1,
       cex = 0.75)
dev.off()


#################################################################################
############################  LGG Loadings Plot for 12.3  #######################
#################################################################################

# Didn't end up using this plot. 

load("~/PanTCGA/mod.SVD.v2.rda") # loads mod.SVD 
# load("~/PanTCGA/pan.fac.results_v2.rda") # loads in the factorization data

# Saving just the 12th module
mod.12 = mod.SVD[[12]]

# Plotting the 3rd right singular vector (which is under mod.12$v in the 3rd column)
# This is plotting all the right singular vectors for all cancer types
# but this module was only based on LGG. 
hist(mod.12$v[,3], breaks = 100)

#################################################################################
##################### Next: Kidney Cancers          #############################
#################################################################################

# Loading in the new dataset from Ricketts (2018) TCGA Comprehensive Renal Cell Carcinoma
kidney_subtypes2 = read.csv(paste(currentwd, "mmc2.csv", sep = ""))
kidney_subtypes2.2 = read.csv(paste(currentwd, "mmc2.csv", sep = "")) # loading in a second copy to check that changing the labels didn't change the data

# Changing the labeling scheme for KIRP so that it matches the language used in the article. 
levels(kidney_subtypes2$PanKidney.Pathology)[levels(kidney_subtypes2$PanKidney.Pathology) == "PRCC T1"] <- "KIRP: Type I"
levels(kidney_subtypes2$PanKidney.Pathology)[levels(kidney_subtypes2$PanKidney.Pathology) == "PRCC T2"] <- "KIRP: Type II"
levels(kidney_subtypes2$PanKidney.Pathology)[levels(kidney_subtypes2$PanKidney.Pathology) == "KIRP CIMP"] <- "KIRP: CIMP"
levels(kidney_subtypes2$PanKidney.Pathology)[levels(kidney_subtypes2$PanKidney.Pathology) == "PRCC Unc"] <- "KIRP: Unclassified "

# Checking that the entire dataset is not changed because these factor levels changed. 
# Removing the column for pathology because that is obviously different.
pathology.col.num <- which(colnames(kidney_subtypes2) == "PanKidney.Pathology")
all.equal(kidney_subtypes2[kidney_subtypes2$PanKidney.Pathology == "KIRP: Type I", -pathology.col.num],
          kidney_subtypes2.2[kidney_subtypes2.2$PanKidney.Pathology == "PRCC T1", -pathology.col.num])
all.equal(kidney_subtypes2[kidney_subtypes2$PanKidney.Pathology == "KIRP: Type II", -pathology.col.num],
          kidney_subtypes2.2[kidney_subtypes2.2$PanKidney.Pathology == "PRCC T2", -pathology.col.num])
all.equal(kidney_subtypes2[kidney_subtypes2$PanKidney.Pathology == "KIRP: CIMP", -pathology.col.num],
          kidney_subtypes2.2[kidney_subtypes2.2$PanKidney.Pathology == "KIRP CIMP", -pathology.col.num])
all.equal(kidney_subtypes2[kidney_subtypes2$PanKidney.Pathology == "PRCC Unc", -pathology.col.num],
          kidney_subtypes2.2[kidney_subtypes2.2$PanKidney.Pathology == "KIRP: Unclassified", -pathology.col.num])

# Checking that the same subjects with the previous label have the new label. Should all be true. 
all(which(kidney_subtypes2$PanKidney.Pathology == "KIRP: Type I") == which(kidney_subtypes2.2$PanKidney.Pathology == "PRCC T1"))
all(which(kidney_subtypes2$PanKidney.Pathology == "KIRP: Type II") == which(kidney_subtypes2.2$PanKidney.Pathology == "PRCC T2"))
all(which(kidney_subtypes2$PanKidney.Pathology == "KIRP: CIMP") == which(kidney_subtypes2.2$PanKidney.Pathology == "KIRP CIMP"))
all(which(kidney_subtypes2$PanKidney.Pathology == "KIRP: Unclassified") == which(kidney_subtypes2.2$PanKidney.Pathology == "PRCC Unc"))

# Subsetting the clinical data by TCGA cancer type. 
kidney_subtypes2.kich = kidney_subtypes2[kidney_subtypes2$Original.TCGA.project == "KICH", ]
kidney_subtypes2.kirc = kidney_subtypes2[kidney_subtypes2$Original.TCGA.project == "KIRC", ]
kidney_subtypes2.kirp = kidney_subtypes2[kidney_subtypes2$Original.TCGA.project == "KIRP", ]

# Analysis ids by cancer type
kidney.analysis.ids = sapply(Kidney.Covariates, function(type) rownames(type))
kidney.analysis.ids.kich = kidney.analysis.ids$KICH
kidney.analysis.ids.kirc = kidney.analysis.ids$KIRC
kidney.analysis.ids.kirp = kidney.analysis.ids$KIRP

# Clinical ids by cancer type. 
kidney.clinical.ids = kidney_subtypes2$bcr_patient_barcode
kidney.clinical.ids.kich = kidney_subtypes2.kich$bcr_patient_barcode
kidney.clinical.ids.kirc = kidney_subtypes2.kirc$bcr_patient_barcode
kidney.clinical.ids.kirp = kidney_subtypes2.kirp$bcr_patient_barcode

# The IDs that are available in both datasets by cancer type
length(kidney.analysis.ids.kich)
length(kidney.clinical.ids.kich)
kidney.avail.ids.kich = kidney.analysis.ids.kich[kidney.analysis.ids.kich %in% kidney.clinical.ids.kich]

length(kidney.analysis.ids.kirc)
length(kidney.clinical.ids.kirc)
kidney.avail.ids.kirc = kidney.analysis.ids.kirc[kidney.analysis.ids.kirc %in% kidney.clinical.ids.kirc]

length(kidney.analysis.ids.kirp)
length(kidney.clinical.ids.kirp)
kidney.avail.ids.kirp = kidney.analysis.ids.kirp[kidney.analysis.ids.kirp %in% kidney.clinical.ids.kirp]

# Combining all the IDs together.
kidney.avail.ids = list(KICH = kidney.avail.ids.kich,
                        KIRC = kidney.avail.ids.kirc,
                        KIRP = kidney.avail.ids.kirp)

# Subsetting analysis data to include just the available ids in both datasets
Kidney.Covariates.Avail = sapply(c("KICH", "KIRC", "KIRP"), function(type) {
  Kidney.Covariates[[type]][rownames(Kidney.Covariates[[type]]) %in% kidney.avail.ids[[type]], ]
})

# Making sure the naming of the list is retained
names(Kidney.Covariates.Avail)

# Subsetting the clinical data for just the IDs avilable in both datasets
kidney_subtypes2.kich.avail = kidney_subtypes2.kich[kidney_subtypes2.kich$bcr_patient_barcode %in% kidney.avail.ids.kich, ]
kidney_subtypes2.kirc.avail = kidney_subtypes2.kirc[kidney_subtypes2.kirc$bcr_patient_barcode %in% kidney.avail.ids.kirc, ]
kidney_subtypes2.kirp.avail = kidney_subtypes2.kirp[kidney_subtypes2.kirp$bcr_patient_barcode %in% kidney.avail.ids.kirp, ]

# Checking the same number of IDs are in the subsetted data
nrow(kidney_subtypes2.kich.avail) == nrow(Kidney.Covariates.Avail[["KICH"]])
nrow(kidney_subtypes2.kirc.avail) == nrow(Kidney.Covariates.Avail[["KIRC"]])
nrow(kidney_subtypes2.kirp.avail) == nrow(Kidney.Covariates.Avail[["KIRP"]])

# Matching the order of the clinical data to the analysis data
kidney_subtypes2.kich.avail.reorder = kidney_subtypes2.kich.avail[match(kidney.avail.ids.kich, kidney_subtypes2.kich.avail$bcr_patient_barcode), ]
kidney_subtypes2.kirc.avail.reorder = kidney_subtypes2.kirc.avail[match(kidney.avail.ids.kirc, kidney_subtypes2.kirc.avail$bcr_patient_barcode), ]
kidney_subtypes2.kirp.avail.reorder = kidney_subtypes2.kirp.avail[match(kidney.avail.ids.kirp, kidney_subtypes2.kirp.avail$bcr_patient_barcode), ]

# Checking that order matches within the dataset. Should all be true.
all(kidney_subtypes2.kich.avail.reorder$bcr_patient_barcode == rownames(Kidney.Covariates.Avail$KICH))
all(kidney_subtypes2.kirc.avail.reorder$bcr_patient_barcode == rownames(Kidney.Covariates.Avail$KIRC))
all(kidney_subtypes2.kirp.avail.reorder$bcr_patient_barcode == rownames(Kidney.Covariates.Avail$KIRP))

# Making a list of all the clinical datasets
kidney_subtypes2.avail.reorder = list(KICH = kidney_subtypes2.kich.avail.reorder,
                                      KIRC = kidney_subtypes2.kirc.avail.reorder,
                                      KIRP = kidney_subtypes2.kirp.avail.reorder)

# Adding subtypes to covariates data
Kidney.Covariates.Avail.Subtypes = lapply(c("KICH", "KIRC", "KIRP"), function(type) {
  kidney_data = cbind.data.frame(Kidney.Covariates.Avail[[type]], 
                                 pathology = kidney_subtypes2.avail.reorder[[type]]$PanKidney.Pathology,
                                 cancer = type)
  kidney_data$pathology = droplevels(kidney_data$pathology)
  kidney_data
})

names(Kidney.Covariates.Avail.Subtypes) = names(Kidney.Covariates.Avail)

############################################################################
############################# Plotting #####################################
############################################################################

# Colors to use for plotting
par_col_kidney <- viridis::viridis(4) 
# Changing the yellow to an ochre so it is more visible
par_col_kidney[par_col_kidney == "#FDE725FF"] <- "#A87C07" 

par_col_kirp <- par_col_kidney[Kidney.Covariates.Avail.Subtypes[["KIRP"]]$pathology] # selecting colors
par_col_kirc <- par_col_kidney[Kidney.Covariates.Avail.Subtypes[["KIRC"]]$pathology] # selecting colors
par_col_kich <- par_col_kidney[Kidney.Covariates.Avail.Subtypes[["KICH"]]$pathology] # selecting colors

###############################################################################
### First plotting each cancer separately with their corresponding subtypes ### 
###############################################################################

###############################################################################
################################   KIRP  ######################################
###############################################################################

png(paste(currentwd, "KIRP11.1vs22.1.png", sep = ""))
plot(Kidney.Covariates.Avail.Subtypes$KIRP[,"11.1"], 
     Kidney.Covariates.Avail.Subtypes$KIRP[,"22.1"],
     col = par_col_kirp,
     pch = 16,
     xlab = "Component 11.1",
     ylab = "Component 22.1",
     main = "Sample Scores for KIRP \n Labeled by Subtype")
legend("bottomright", legend = levels(Kidney.Covariates.Avail.Subtypes$KIRP$pathology),
       col = par_col_kidney, pch = 16, cex = 0.75)
dev.off()

### KDE plots
png(paste(currentwd, "KIRP_KDE22.1.png", sep = ""))
densityplot(~c(Kidney.Covariates.Avail.Subtypes$KIRP[, "22.1"]), 
            groups = Kidney.Covariates.Avail.Subtypes$KIRP[, "pathology"], 
            xlab = "Component 22.1",
            main = "Component 22.1 Scores for KIRP \n by Subtype",
            par.settings = list(superpose.line = list(col=par_col_kidney)),
            auto.key = TRUE,
            plot.points=FALSE)
dev.off()

png(paste(currentwd, "KIRP_KDE11.1.png", sep = ""))
densityplot(~c(Kidney.Covariates.Avail.Subtypes$KIRP[, "11.1"]), 
            groups = Kidney.Covariates.Avail.Subtypes$KIRP[, "pathology"], 
            xlab = "Component 11.1",
            main = "Component 11.1 Scores for KIRP \n by Subtype",
            par.settings = list(superpose.line = list(col=par_col_kidney)),
            auto.key = TRUE,
            plot.points=FALSE)
dev.off()

###############################################################################
###########################   Pairwise Plots   ################################
###############################################################################

###############################################################################
###########################   KIRP & KIRC   ###################################
###############################################################################

# Didn't end up using these. 

# Getting rid of the subjects who were weirdly classified as KICH
Kidney.Covariates.Avail.KIRC.Without.ChRCC = Kidney.Covariates.Avail.Subtypes$KIRC[Kidney.Covariates.Avail.Subtypes$KIRC$pathology != "ChRCC", ]
Kidney.Covariates.Avail.KIRC.Without.ChRCC$pathology = factor("KIRC")

# Combining KIRP and KIRC without ChRCC together
Kidney.Covariates.Avail.Subtypes.KIRP.KIRC = list(Kidney.Covariates.Avail.KIRC.Without.ChRCC,
                                                  Kidney.Covariates.Avail.Subtypes$KIRP)

# Selecting just the covariates of interest
Kidney.Covariates.Avail.Subtypes.KIRP.KIRC.Selected = lapply(Kidney.Covariates.Avail.Subtypes.KIRP.KIRC, function(type) {
  type[,c("11.1", "22.1", "pathology", "cancer")]
})

# Combining everything into one matrix
Kidney.Covariates.Avail.Subtypes.KIRP.KIRC.Selected.Combined = do.call(rbind, Kidney.Covariates.Avail.Subtypes.KIRP.KIRC.Selected)
# Getting rid of levels that don't exist
Kidney.Covariates.Avail.Subtypes.KIRP.KIRC.Selected.Combined$pathology = droplevels(Kidney.Covariates.Avail.Subtypes.KIRP.KIRC.Selected.Combined$pathology)

# Setting colors
par_col_kidney_5 <- viridis::viridis(5) 
# Changing the yellow to an ochre
par_col_kidney_5[par_col_kidney_5 == "#FDE725FF"] <- "#A87C07" 

# Assigning colors to subgroups
par_col_KIRC_KIRP_plotting <- par_col_kidney_5[Kidney.Covariates.Avail.Subtypes.KIRP.KIRC.Selected.Combined$pathology]
pch_KIRC_KIRP_plotting <- c(8, 19)[as.numeric(Kidney.Covariates.Avail.Subtypes.KIRP.KIRC.Selected.Combined$cancer)]

png(paste(currentwd, "KIRPKIRC_11.1vs22.1.png", sep = ""), 
    width = 500, height = 500, units = "px")
plot(Kidney.Covariates.Avail.Subtypes.KIRP.KIRC.Selected.Combined[,"11.1"], 
     Kidney.Covariates.Avail.Subtypes.KIRP.KIRC.Selected.Combined[,"22.1"],
     col = par_col_KIRC_KIRP_plotting,
     xlab = "11.1",
     ylab = "22.1",
     pch = pch_KIRC_KIRP_plotting,
     main = "KIRP & KIRC Scores by Subtype")
legend("bottomright", 
       legend = levels(Kidney.Covariates.Avail.Subtypes.KIRP.KIRC.Selected.Combined$pathology),
       col = par_col_kidney_5, pch = c(8, rep(19, 4)), cex = 0.65)
dev.off()


###############################################################################
###########################   KIRP & KICH   ###################################
###############################################################################

# Also didn't end up using these. 

# Changing the factor names for KICH subjects to KICH
Kidney.Covariates.Avail.Subtypes$KICH$pathology = "KICH"

# Combining KIRP and KICH data together
Kidney.Covariates.Avail.Subtypes.KIRP.KICH = list(Kidney.Covariates.Avail.Subtypes$KIRP,
                                                  Kidney.Covariates.Avail.Subtypes$KICH)

# Selecting just the BIDIFAC factors of interest
Kidney.Covariates.Avail.Subtypes.KIRP.KICH.Selected = lapply(Kidney.Covariates.Avail.Subtypes.KIRP.KICH, function(type) {
  type[,c("24.1", "22.1", "pathology", "cancer")]
})

# Putting everything into one matrix
Kidney.Covariates.Avail.Subtypes.KIRP.KICH.Selected.Combined = do.call(rbind, Kidney.Covariates.Avail.Subtypes.KIRP.KICH.Selected)

# Setting the colors
par_col_KIRP_KICH_plotting <- par_col_kidney_5[Kidney.Covariates.Avail.Subtypes.KIRP.KICH.Selected.Combined$pathology]
pch_KIRP_KICH_plotting <- c(19, 8)[as.numeric(Kidney.Covariates.Avail.Subtypes.KIRP.KICH.Selected.Combined$cancer)]

# Plotting
png(paste(currentwd, "KIRPKICH_22.1vs24.1.png", sep = ""), 
    width = 500, height = 500, units = "px")
plot(Kidney.Covariates.Avail.Subtypes.KIRP.KICH.Selected.Combined[,"22.1"], 
     Kidney.Covariates.Avail.Subtypes.KIRP.KICH.Selected.Combined[,"24.1"],
     col = par_col_KIRP_KICH_plotting,
     xlab = "22.1",
     ylab = "24.1",
     pch = pch_KIRP_KICH_plotting,
     main = "KIRP & KICH Scores by Subtype")
legend("topleft", 
       legend = levels(Kidney.Covariates.Avail.Subtypes.KIRP.KICH.Selected.Combined$pathology),
       col = par_col_kidney_5, pch = c(rep(19, 4), 8), cex = 0.65)
dev.off()

###############################################################################
###########################   KIRC & KICH   ###################################
###############################################################################

# Also didn't end up using these. 

Kidney.Covariates.Avail.Subtypes.KIRC.KICH = list(Kidney.Covariates.Avail.Subtypes$KIRC,
                                                  Kidney.Covariates.Avail.Subtypes$KICH)

Kidney.Covariates.Avail.Subtypes.KIRC.KICH.Selected = lapply(Kidney.Covariates.Avail.Subtypes.KIRC.KICH, function(type) {
  type[,c("14.1", "22.1", "pathology", "cancer")]
})

Kidney.Covariates.Avail.Subtypes.KIRC.KICH.Selected.Combined = do.call(rbind, Kidney.Covariates.Avail.Subtypes.KIRC.KICH.Selected)

ggplot(Kidney.Covariates.Avail.Subtypes.KIRC.KICH.Selected.Combined) + 
  geom_point(aes(x=`14.1`,y=`22.1`,color=factor(pathology),shape=cancer),size=2)


###############################################################################
###########################   All 3 together   ################################
###############################################################################
Kidney.Covariates.Avail.Subtypes.Selected = lapply(Kidney.Covariates.Avail.Subtypes, function(type) {
  type[,c("14.1", "22.1", "pathology", "cancer")]
})

Kidney.Covariates.Avail.Subtypes.Selected.Combined = do.call(rbind, Kidney.Covariates.Avail.Subtypes.Selected)

ggplot(Kidney.Covariates.Avail.Subtypes.Selected.Combined) + 
  geom_point(aes(x=`14.1`,y=`22.1`,color=factor(pathology),shape=cancer),size=2)

############################################################################
####################### Kidney Survival Plots ##############################
############################################################################

# Saving the survival times
Kidney.Survival = Survival[c("KICH", "KIRC", "KIRP")]
Kidney.Censored = Censored[c("KICH", "KIRC", "KIRP")]

# Combining the event indicator and survival times together
Kidney.Covariates.Avail.Subtypes.Surv = lapply(c("KICH", "KIRC", "KIRP"), function(type) {
  # storing the survival and censored times for the current kidney type
  type.survival = Kidney.Survival[[type]][names(Kidney.Survival[[type]]) %in% kidney.avail.ids[[type]]]
  type.censored = Kidney.Censored[[type]][names(Kidney.Censored[[type]]) %in% kidney.avail.ids[[type]]]
  
  # saving the event indicator 
  type.event = !is.na(type.survival)
  
  # storing the censored times for the censored events
  type.survival[is.na(type.survival)] = type.censored[is.na(type.survival)]
  
  kidney.data = cbind.data.frame(Kidney.Covariates.Avail.Subtypes[[type]], survival = type.survival, event = type.event)
  
  if (type == "KIRC") {
    kidney.data = kidney.data[kidney.data$pathology != "ChRCC",]
    kidney.data$pathology = factor("KIRC")
  }
  if (type == "KICH") {
    kidney.data$pathology = factor("KICH")
  }
  
  kidney.data
})

names(Kidney.Covariates.Avail.Subtypes.Surv) = names(Kidney.Covariates.Avail.Subtypes)

# Storing just the survival and cancer type data, no covariates
Kidney.Covariates.Avail.Subtypes.Surv.JustSurvival = lapply(Kidney.Covariates.Avail.Subtypes.Surv, function(type) {
  type[, c("pathology", "cancer", "survival", "event")]
})

# Combining before plotting
Kidney.Covariates.Avail.Subtypes.Surv.Combined = do.call(rbind, Kidney.Covariates.Avail.Subtypes.Surv.JustSurvival)

### Combined cancer survival plots
par_col_kidney_6 <- viridis::viridis(6) # colors to use
par_col_kidney_6[par_col_kidney_6 == "#FDE725FF"] <- "#A87C07" 

# Didn't end up using these plot. 

png(paste(currentwd, "KidneySurvivalPlots.png", sep = ""), 
    width = 500, height = 500, units = "px")
plot(survfit(Surv(survival, event, type = "right") ~ pathology, 
             data = Kidney.Covariates.Avail.Subtypes.Surv.Combined),
     xlab = "Days",
     ylab = "Overall Survival Probability",
     col = par_col_kidney_6,
     main = "Kaplan-Meier Curves for Kidney Cancer Patients")
legend("bottomright", legend = levels(Kidney.Covariates.Avail.Subtypes.Surv.Combined$pathology), 
       col=par_col_kidney_6, lty=1,
       cex = 0.5)
dev.off()

###########################################################################
############### Survival curves for just KIRP and KIRC. ###################
###########################################################################

# Combining KIRP and KIRC data together
Kidney.Covariates.Avail.Subtypes.Surv.KIRP.KIRC.Combined = do.call(rbind, Kidney.Covariates.Avail.Subtypes.Surv.JustSurvival[c("KIRP", "KIRC")])

# Setting the colors
par_col_kidney_5 <- viridis::viridis(5) # colors to use
par_col_kidney_5[par_col_kidney_5 == "#FDE725FF"] <- "#A87C07" # changing yellow to ochre

png(paste(currentwd, "KIRPKIRCSurvivalPlots.png", sep = ""), 
    width = 500, height = 500, units = "px")
plot(survfit(Surv(survival, event, type = "right") ~ pathology, 
             data = Kidney.Covariates.Avail.Subtypes.Surv.KIRP.KIRC.Combined),
     xlab = "Days",
     ylab = "Overall Survival Probability",
     col = par_col_kidney_5,
     main = "Kaplan-Meier Curves for KIRP and KIRC Patients by Subtype")
legend("bottomright", legend = levels(Kidney.Covariates.Avail.Subtypes.Surv.KIRP.KIRC.Combined$pathology), 
       col=par_col_kidney_5, lty=1,
       cex = 0.75)
dev.off()

###########################################################################
############### Survival curves for just KICH and KIRC. ###################
############### by value of scores.                     ###################
###########################################################################

KICH.22.1.Groups <- Kidney.Covariates.Avail.Subtypes.Surv[["KICH"]][,"22.1"] >= 0
KICH.22.1.Groups <- factor(ifelse(KICH.22.1.Groups, "KICH: Greater Than 0", "KICH: Less Than 0"))

KIRC.22.1.Groups <- Kidney.Covariates.Avail.Subtypes.Surv[["KIRC"]][,"22.1"] >= 0
KIRC.22.1.Groups <- factor(ifelse(KIRC.22.1.Groups, "KIRC: Greater Than 0", "KIRC: Less Than 0"))

# Putting the data together
Survival.KICH.22.1Groups.Avail = cbind.data.frame(survival = Kidney.Covariates.Avail.Subtypes.Surv[["KICH"]][["survival"]],
                                                  group = KICH.22.1.Groups,
                                                  cancer = "KICH",
                                                  event = Kidney.Covariates.Avail.Subtypes.Surv[["KICH"]][["event"]])
Survival.KIRC.22.1Groups.Avail = cbind.data.frame(survival = Kidney.Covariates.Avail.Subtypes.Surv[["KIRC"]][["survival"]],
                                                  group = KIRC.22.1.Groups,
                                                  cancer = "KIRC",
                                                  event = Kidney.Covariates.Avail.Subtypes.Surv[["KIRC"]][["event"]])
Survival.22.1Groups.Avail = rbind.data.frame(Survival.KICH.22.1Groups.Avail,
                                             Survival.KIRC.22.1Groups.Avail)


# Making a survival plot
par_col_kich_kirc <- par_col_kidney[c(1,1,2,2)]
lty_kich_kirc <- c(1,2,1,2)

png(paste(currentwd, "KICHKIRC_22.1Groups_Survival.png", sep = ""))
plot(survfit(Surv(survival,event, type = "right") ~ group, 
             data = Survival.22.1Groups.Avail),
     xlab = "Days",
     ylab = "Overall Survival Probability",
     col = par_col_kich_kirc,
     lty = lty_kich_kirc,
     main = "Kaplan-Meier Curves for KICH & KIRC Patients \n by Component 22.1 Scores")
legend("bottomleft", legend = levels(Survival.22.1Groups.Avail$group), 
       col=par_col_kich_kirc, lty=lty_kich_kirc,
       cex = 0.75)
dev.off()

# Checking that lines and colors match correct cancers and groups:
# plot(survfit(Surv(survival,event, type = "right") ~ group, 
#              data = Survival.22.1Groups.Avail),
#      xlab = "Days",
#      ylab = "Overall Survival Probability",
#      col = par_col_kidney,
#      main = "Kaplan-Meier Curves for KICH & KIRC Patients \n by Value of Component 22.1")
# legend("bottomleft", legend = levels(Survival.22.1Groups.Avail$group), 
#        col=par_col_kidney, lty=1,
#        cex = 0.75)
