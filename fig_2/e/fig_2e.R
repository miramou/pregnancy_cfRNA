#Get ensemble predictions

require(caret)
require(randomForest)
require(plyr)
require(dplyr)
require(ggplot2)

#Set working directory
mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

#Source scripts
source('plot_theme.R')
source('load_data.R')
source('fig_2/e/fig_2e_helper.R')

#Load model
rf_dCD = readRDS('fig_2/e/models/random_forest_no_cma.rds')

#Source data
load('counts_gene_lists.RData')

#Process data
counts_all = load_process_qPCR(counts_per_mL_all, genes_6922AC, is_classification = FALSE, filter_5695 = FALSE)
counts_all$cohort = factor(counts_all$cohort)

#Split data
counts_to_train = counts_all %>% filter(cohort == 'Danish', delivery >= 37)

num_patient= length(unique(counts_to_train$patient))
set.seed(139)
patient_split = sample.int(num_patient, size=num_patient*0.7)
patient_train = unique(counts_to_train$patient)[patient_split]

train = counts_to_train[counts_to_train$patient %in% patient_train,] 
test = counts_to_train[!(counts_to_train$patient %in% patient_train),]

#Organize training data
rf_dCD_train = rf_dCD$pred %>% group_by(rowIndex) %>% filter(Variables == rf_dCD$optsize) %>% summarize_if(is.numeric, funs(mean(.))) %>% select(pred, obs)
rf_dCD_train$ga = train$ga
rf_dCD_train$est_dCD = train$ga - 40
rf_dCD_train$patient = train$patient

#Organize validation data
rf_dCD_pred = data.frame(pred = round(predict(rf_dCD, newdata = test)), obs = round(test$dCD), est_dCD = (test$ga - 40), ga = test$ga, patient = test$patient)

#Calculate expected delivery date predictions
t2 = rbind(estimate_delivery(rf_dCD_pred, trimester = 2), estimate_delivery(rf_dCD_train, trimester = 2))
t3 = rbind(estimate_delivery(rf_dCD_pred, trimester = 3), estimate_delivery(rf_dCD_train, trimester = 3))
both = rbind(estimate_delivery(rf_dCD_pred, trimester = 4), estimate_delivery(rf_dCD_train, trimester = 4))

#Organize as dataframe for plotting and stats
t2_df = data.frame(condition = c(rep('cfRNA',dim(t2)[1]), rep('Ultrasound', dim(t2)[1])), value = c(t2$delta, t2$us_delta))
t3_df = data.frame(condition = c(rep('cfRNA',dim(t3)[1]), rep('Ultrasound', dim(t3)[1])), value = c(t3$delta, t3$us_delta))
both_df = data.frame(condition = c(rep('cfRNA',dim(both)[1]), rep('Ultrasound', dim(both)[1])), value = c(both$delta, both$us_delta))

#Percentages within 1, 2, or > 2 weeks
t2_stats_cf = get_stats(t2_df, t2, 'cfRNA')
t3_stats_cf = get_stats(t3_df, t3, 'cfRNA')
both_stats_cf = get_stats(both_df, both, 'cfRNA')
both_stats_us = get_stats(both_df, both, 'Ultrasound')

#Plot
group.colors = c("cfRNA" = "#1f78b4", "Ultrasound" = "#33a02c")

plot_t2 = ggplot(t2_df) +
  geom_density(aes(x=value, fill = condition), alpha =0.4) +
  scale_fill_manual(values = group.colors) +
  theme_pub() 

pdf('fig_2/e/fig_2e_t2.pdf', useDingbats = FALSE, width = 2, height=2)
print(plot_t2)
dev.off()

plot_t3 = ggplot(t3_df) +
  geom_density(aes(x=value, fill = condition), alpha =0.4) +
  scale_fill_manual(values = group.colors) +
  theme_pub()

pdf('fig_2/e/fig_2e_t3.pdf', useDingbats = FALSE, width = 2, height=2)
print(plot_t3)
dev.off()

plot_both = ggplot(both_df) +
  geom_density(aes(x=value, fill = condition), alpha =0.4) +
  scale_fill_manual(values = group.colors) +
  theme_pub() 

pdf('fig_2/e/fig_2e_both.pdf', useDingbats = FALSE, width = 2, height=2)
print(plot_both)
dev.off()
