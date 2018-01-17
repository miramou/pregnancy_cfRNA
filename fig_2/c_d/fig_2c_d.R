library(ggplot2)
library(caret)
library(randomForest)
library(plyr)
library(dplyr)

#Set working directory
mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

#Source scripts
source('plot_theme.R')
source('load_data.R')

#Load model
rf_dCD = readRDS('fig_2/c_d/models/best_random_forest.rds')

#Source data
load('counts_gene_lists.RData')

#Process data
counts_all = load_process_qPCR(counts_per_mL_all_cma, genes_6922AC, is_classification = FALSE, filter_5695 = FALSE)
counts_all$cohort = factor(counts_all$cohort)

#Split training and Val using same seed as model training
counts_to_train = counts_all %>% filter(cohort == 'Danish', delivery >= 37)

num_patient= length(unique(counts_to_train$patient))
set.seed(139)
patient_split = sample.int(num_patient, size=num_patient*0.7)
patient_train = unique(counts_to_train$patient)[patient_split]

train = counts_to_train[counts_to_train$patient %in% patient_train,]
test = counts_to_train[!(counts_to_train$patient %in% patient_train),]

#Organize data to plot for training and validation

#Training
rf_dCD_train = rf_dCD$pred
rf_dCD_train = rf_dCD_train %>% group_by(rowIndex) %>% filter(Variables == rf_dCD$optsize) %>% summarize_if(is.numeric, funs(mean(.)))
rf_dCD_train$ga = train$ga

#Validation
rf_dCD_pred = data.frame(pred = round(predict(rf_dCD, newdata = test)), obs = round(test$dCD), ga = test$ga, cohort = test$cohort)

#Plot training
plot_dCD_train = ggplot(rf_dCD_train) +
  geom_point(aes(x = obs, y = pred), color = 'black') +
  geom_smooth(aes(x = obs, y= pred), method = lm, fill = 'lightblue', color = 'darkblue') +
  labs(x = 'Time to delivery (weeks)', y = 'RF model value (weeks)') +
  theme(
    panel.grid.major.x = element_line(color = 'black'),
    panel.grid.major.y = element_line(color = 'black')) +
  theme_pub() +
  coord_cartesian(xlim = c(-33, 10), ylim = c(-33, 10))

pdf('fig_2/c_d/fig_2c.pdf', useDingbats = FALSE, width = 3, height = 3)
print(plot_dCD_train)
dev.off()

#Plot validation
plot_dCD_val = ggplot(rf_dCD_pred) +
  geom_point(aes(x = obs, y = pred), color = 'black') +
  geom_smooth(aes(x = obs, y= pred), method = lm, fill = 'lightblue', color = 'darkblue') +
  labs(x = 'Time to delivery (weeks)', y = 'RF model value (weeks)') +
  theme(
    panel.grid.major.x = element_line(color = 'black'),
    panel.grid.major.y = element_line(color = 'black')) +
  theme_pub() +
  coord_cartesian(xlim = c(-33, 10), ylim = c(-33, 10))

pdf('fig_2/c_d/fig_2d.pdf', useDingbats = FALSE, width = 3, height = 3)
print(plot_dCD_val)
dev.off()

#Get RMSE values
train_RMSE_trimester = c(RMSE(rf_dCD_train$pred[rf_dCD_train$ga <= 12], rf_dCD_train$obs[rf_dCD_train$ga <= 12]), #T1
                         RMSE(rf_dCD_train$pred[rf_dCD_train$ga > 12 & rf_dCD_train$ga <= 24], rf_dCD_train$obs[rf_dCD_train$ga > 12 & rf_dCD_train$ga <= 24]), #T2 
                         RMSE(rf_dCD_train$pred[rf_dCD_train$ga > 24 & rf_dCD_train$obs <= 0], rf_dCD_train$obs[rf_dCD_train$ga > 24 & rf_dCD_train$obs <= 0]), #T3
                         RMSE(rf_dCD_train$pred[rf_dCD_train$ga > 24 & rf_dCD_train$obs > 0], rf_dCD_train$obs[rf_dCD_train$ga > 24 & rf_dCD_train$obs > 0])) #PP

test_RMSE_trimester = c(RMSE(rf_dCD_pred$pred[rf_dCD_pred$ga <= 12], rf_dCD_pred$obs[rf_dCD_pred$ga <= 12]), #T1
                        RMSE(rf_dCD_pred$pred[rf_dCD_pred$ga > 12 & rf_dCD_pred$ga <= 24], rf_dCD_pred$obs[rf_dCD_pred$ga > 12 & rf_dCD_pred$ga <= 24]), #T2
                        RMSE(rf_dCD_pred$pred[rf_dCD_pred$ga > 24 & rf_dCD_pred$obs <= 0], rf_dCD_pred$obs[rf_dCD_pred$ga > 24 & rf_dCD_pred$obs <= 0]), #T3
                        RMSE(rf_dCD_pred$pred[rf_dCD_pred$ga > 24 & rf_dCD_pred$obs > 0], rf_dCD_pred$obs[rf_dCD_pred$ga > 24 & rf_dCD_pred$obs > 0])) #PP

#Get pearson correlations
cor.test(rf_dCD_train$pred, rf_dCD_train$obs)
cor.test(rf_dCD_pred$pred, rf_dCD_pred$obs)
