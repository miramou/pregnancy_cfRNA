#Build Models

require(caret)
require(plyr)
require(dplyr)

#Set working directory
mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

setwd('/Users/miramou/Documents/Grad/Rotations/Quake/Manuscript/code_to_share')

#Source scripts
source('load_data.R')

#Source data
load('counts_gene_lists.RData')

#Process data

#For RFE want genes that we measured across all 3 cohorts hence genes_5695_6922AC
counts_all_rfe = load_process_qPCR(counts_per_mL_all_cma, genes_5695_6922AC, is_classification = FALSE, filter_5695 = FALSE) 
counts_all_rfe$cohort = factor(counts_all_rfe$cohort)
counts_to_train = counts_all_rfe %>% filter(cohort == 'Danish', delivery >= 37)

#For full model want all genes available so genes_5695_6922C
counts_all_full = load_process_qPCR(counts_per_mL_all_cma, genes_5695_6922C, is_classification = FALSE, filter_5695 = FALSE) 
counts_all_full$cohort = factor(counts_all_full$cohort)
counts_to_train_full = counts_all_full %>% filter(cohort == 'Danish', delivery >= 37)

#For part e want no cma for accurate averaging
counts_all_no_cma = load_process_qPCR(counts_per_mL_all, genes_5695_6922AC, is_classification = FALSE, filter_5695 = FALSE) 
counts_all_no_cma$cohort = factor(counts_all_no_cma$cohort)
counts_to_train_no_cma = counts_all_no_cma %>% filter(cohort == 'Danish', delivery >= 37)

#Split data
num_patient= length(unique(counts_to_train$patient))
set.seed(139)
patient_split = sample.int(num_patient, size=num_patient*0.7)
patient_train = unique(counts_to_train$patient)[patient_split]

train = counts_to_train[counts_to_train$patient %in% patient_train,] 
train_full = counts_to_train_full[counts_to_train_full$patient %in% patient_train,] 
train_no_cma = counts_to_train_no_cma[counts_to_train_no_cma$patient %in% patient_train,]

#Train random forest models for time to delivery
subsets = c(1:9, 20, 25)

newPickSizeTolerance = function(x, metric, tol=2.5, maximize) {
  return(caret::pickSizeTolerance(x, metric, tol, maximize))
}
rfFuncs$selectSize = newPickSizeTolerance

ctrl_rf = rfeControl(functions = rfFuncs, method = 'repeatedcv', repeats = 10, verbose = FALSE, saveDetails = TRUE)
ctrl_rf_full = trainControl(method = 'repeatedcv', repeats = 10, savePredictions = TRUE)

#Note on randomness and model training:
#Random forests still have randomness in how the trees are built and what splits of the data are chosen. 
#These are not controlled for here therefore this will not exactly reproduce the RF model used in the paper.
#To do so, use the saved RDS RF model provided.

#rf model for predicting time to delivery
set.seed(139)
rf_dCD = rfe(train %>% select(-starts_with('ga'), -starts_with('dCD'),-starts_with('delivery'), -starts_with('sample_id'),-starts_with('panel'), -starts_with('patient'), -starts_with('cohort')), 
             train$dCD, sizes = subsets, rfeControl = ctrl_rf)

#rf model full to compare
set.seed(139)
rf_dCD_full_matched = train(train_full %>% select(-starts_with('ga'), -starts_with('dCD'), -starts_with('sample_id'),-starts_with('delivery'), -starts_with('panel'), -starts_with('patient'), -starts_with('cohort')), 
                    train_full$dCD, method = 'rf', trControl = ctrl_rf_full)

#rf model for predicting expected delivery date
subsets = c(9) # Generate the same model as above but with no centered-moving average applied to the data
newPickSizeTolerance = function(x, metric, tol=0.5, maximize) {
  return(caret::pickSizeTolerance(x, metric, tol, maximize))
}
rfFuncs$selectSize = newPickSizeTolerance

set.seed(139)
rf_dCD_no_cma = rfe(train_no_cma %>% select(-starts_with('ga'), -starts_with('dCD'),-starts_with('delivery'), -starts_with('sample_id'),-starts_with('panel'), -starts_with('patient'), -starts_with('cohort')), 
                    train_no_cma$dCD, sizes = subsets, rfeControl = ctrl_rf)

#Save models
saveRDS(rf_dCD, 'fig_2/c_d/models/best_random_forest.rds')
saveRDS(rf_dCD_full_matched, 'fig_2/c_d/models/full_random_forest.rds')
saveRDS(rf_dCD_no_cma, 'fig_2/e/model/random_forest_no_cma.rds')