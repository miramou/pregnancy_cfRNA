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

#Source data
load('counts_gene_lists.RData')

#Process data
counts_all = load_process_qPCR(counts_per_mL_all, genes_5695_6922C, is_classification = FALSE, filter_5695 = FALSE)
counts_all$cohort = factor(counts_all$cohort)

#Load models for ga and rds
rf_dCD = readRDS('fig_2/c_d/models/best_random_forest.rds')

#Split 30% validation data
counts_to_train = counts_all %>% filter(cohort == 'Danish', delivery >= 37)

num_patient= length(unique(counts_to_train$patient))
set.seed(139)
patient_split = sample.int(num_patient, size=num_patient*0.7)
patient_train = unique(counts_to_train$patient)[patient_split]

train = counts_to_train[counts_to_train$patient %in% patient_train,]
test = counts_to_train[!(counts_to_train$patient %in% patient_train),]

#Predict RF values for each cohort/split of data
rf_dCD_pred = data.frame(pred = round(predict(rf_dCD, newdata = test)), dCD = round(test$dCD), ga = test$ga, cohort = test$cohort)
rf_dCD_pred_penn =  data.frame(pred = round(predict(rf_dCD, newdata = counts_all %>% filter(cohort == 'T_Penn'))), obs = round(counts_all %>% filter(cohort == 'T_Penn') %>% select(dCD)), 
                               cohort = counts_all %>% filter(cohort == 'T_Penn') %>% select(cohort))
rf_dCD_pred_uab =  data.frame(pred = round(predict(rf_dCD, newdata = counts_all %>% filter(cohort == 'T_UAB'))), obs = round(counts_all %>% filter(cohort == 'T_UAB') %>% select(dCD)), 
                              cohort = counts_all %>% filter(cohort == 'T_UAB') %>% select(cohort))
rf_dCD_pred_pt_penn =  data.frame(pred = round(predict(rf_dCD, newdata = counts_all %>% filter(cohort == 'PT_Penn'))), obs = round(counts_all %>% filter(cohort == 'PT_Penn') %>% select(dCD)), 
                                  cohort = counts_all %>% filter(cohort == 'PT_Penn') %>% select(cohort))
rf_dCD_pred_pt_uab =  data.frame(pred = round(predict(rf_dCD, newdata = counts_all %>% filter(cohort == 'PT_UAB'))), obs = round(counts_all %>% filter(cohort == 'PT_UAB') %>% select(dCD)), 
                                 cohort = counts_all %>% filter(cohort == 'PT_UAB') %>% select(cohort))

rmse_pt_t = c(RMSE(c(rf_dCD_pred_pt_penn$pred, rf_dCD_pred_pt_uab$pred), c(rf_dCD_pred_pt_penn$dCD, rf_dCD_pred_pt_uab$dCD)), #PT
              RMSE(c(rf_dCD_pred_penn$pred, rf_dCD_pred_uab$pred), c(rf_dCD_pred_penn$dCD, rf_dCD_pred_uab$dCD))) #T

#Organize data and plot
test_all = rbind(rf_dCD_pred_uab, rf_dCD_pred_pt_uab, rf_dCD_pred_pt_penn, rf_dCD_pred_penn)
test_all$class = test_all$cohort
test_all$class = revalue(test_all$class, c('T_UAB' = 'Term', 'PT_UAB' = 'Preterm', 'T_Penn' = 'Term', 'PT_Penn' = 'Preterm'))
test_all$cohort = revalue(test_all$cohort, c('T_UAB' = 'UAB', 'PT_UAB' = 'UAB', 'T_Penn' = 'UPenn', 'PT_Penn' = 'UPenn')) 

group.colors = c("Term" = "#33a02c", "Preterm" = "#386cb0")

plot_dCD_test = ggplot() +
  geom_point(data = rf_dCD_pred, aes(x = dCD, y=pred), size =1, color = 'black', alpha = 0.2, shape = 15) +
  geom_point(data = test_all, aes(x = dCD, y = pred, shape = cohort, color = class), size = 2.5) +
  scale_color_manual(values = group.colors) +
  geom_smooth(data = test_all, aes(x = dCD, y = pred),method = lm, fill = '#9ecae1', color = '#08519c') +
  labs(x = 'Time until delivery (weeks)', y = 'RF model value (weeks)', shape = 'Cohort') +
  theme(
    panel.grid.major.x = element_line(color = 'black'),
    panel.grid.major.y = element_line(color = 'black')) +
  theme_pub() +
  coord_cartesian(xlim = c(-30, 0), ylim = c(-30, 0))

pdf('supplement/s3/fig_s3.pdf',width = 4, height = 4, useDingbats = FALSE)
print(plot_dCD_test)
dev.off()