#Compare RF full vs small
require(caret)
require(plyr)
require(dplyr)
require(ggplot2)

#Set working directory
mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

#Source scripts
source('plot_theme.R')
source('load_data.R')

#Source data
load('counts_gene_lists.RData')

#Process data
counts_all = load_process_qPCR(counts_per_mL_all, genes_5695_6922C, is_classification = FALSE, filter_5695 = FALSE)
counts_all$cohort = factor(counts_all$cohort)

#Split data
counts_to_train = counts_all %>% filter(cohort == 'Danish', delivery >= 37)

num_patient= length(unique(counts_to_train$patient))
set.seed(139)
patient_split = sample.int(num_patient, size=num_patient*0.7)
patient_train = unique(counts_to_train$patient)[patient_split]

train = counts_to_train[counts_to_train$patient %in% patient_train,] 
test = counts_to_train[!(counts_to_train$patient %in% patient_train),]

#Load models
rf_dCD_small = readRDS('fig_2/c_d/models/best_random_forest.rds')
rf_dCD_full = readRDS('fig_2/c_d/models/full_random_forest.rds')

rmse_df = data.frame(model = c('Nine gene subset', 'Full 51 gene model'), 
                     RMSE = c(rf_dCD_small$results$RMSE[rf_dCD_small$results$Variables==rf_dCD_small$bestSubset],
                              rf_dCD_full$results$RMSE[rf_dCD_full$results$mtry==rf_dCD_full$bestTune$mtry]),
                     RMSESD = c(rf_dCD_small$results$RMSESD[rf_dCD_small$results$Variables==rf_dCD_small$bestSubset],
                               rf_dCD_full$results$RMSESD[rf_dCD_full$results$mtry==rf_dCD_full$bestTune$mtry]))
plot = ggplot(rmse_df) +
  geom_bar(aes(x = model, y = RMSE), stat = 'identity', fill = 'lightgray', color = 'black') +
  geom_errorbar(aes(x = model, ymax = RMSE+RMSESD, ymin = RMSE)) +
  labs(y='Root mean squared error (weeks)') +
  theme_pub() +
  theme(
    #axis.text.x = element_text(angle = 90),
    axis.title.x = element_blank())

pdf('supplement/s2/fig_s2.pdf', useDingbats = FALSE)
print(plot)
dev.off()
