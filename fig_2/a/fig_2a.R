library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

#Set working directory
mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

#Source scripts
source('plot_theme.R')
source('load_data.R')
source('fig_2/a/fig_2a_helper.R')

#Source data
load('counts_gene_lists.RData')

#Process data
counts_all = load_process_qPCR(counts_per_mL_all, genes_5695_6922C, is_classification = FALSE, filter_5695 = FALSE)
counts_21 = load_process_qPCR(counts_per_mL_all, genes_6922C, is_classification = FALSE, filter_5695 = TRUE)
counts_all$cohort = factor(counts_all$cohort)
counts_danish_no_cma = counts_all %>% filter(cohort == 'Danish')
counts_21 = na.omit(counts_21)

patients = unique(counts_danish_no_cma %>% select(patient))
downsample_no_cma = data.frame()
downsample_no_cma = rbind(downsample_no_cma, counts_danish_no_cma[counts_danish_no_cma$patient %in% patients[11:31,1],]) #11-31 sampled every two weeks so no need to downsample

for (p in patients[1:10, 1]) {
  downsample = seq(2, dim(counts_danish_no_cma %>% filter(patient == p))[1], by = 2) #Downsample so that sampling scheme matches for 1-10 and 11-31
  patient_downsample = counts_danish_no_cma %>% filter(patient == p)
  patient_downsample = patient_downsample[downsample,]
  downsample_no_cma = rbind(downsample_no_cma, patient_downsample)
}

all_patient_avg_ga = rbind(avg_wrapper(downsample_no_cma, interval = 4, isGA = TRUE), avg_wrapper(counts_21, interval = 4, isGA = TRUE))
all_patient_avg_dCD = rbind(avg_wrapper(downsample_no_cma, interval = 4, isGA = FALSE), avg_wrapper(counts_21, interval = 4, isGA = FALSE))

#Plot panels
placental = c('ALPP', 'CAPN6','CGA', 'CGB','CSHL1', 'LGALS14', 'PAPPA', 'PLAC4', 'PSG7', 'ADAM12','CSH1', 'CSH2', 'GH2', 'HSD17B1', 
              'HSD3B1', 'HMGB3','HSPB8', 'PSG1', 'PSG4', 'S100P','VGLL1') 
immune = c('ANXA3', 'ARG1', 'CAMP', 'CD24', 'CEACAM8','DEFA3', 'DEFA4', 'ELANE', 'LTF', 'MMP8', 'MPO', 'PGLYRP1', 'S100A8', 'S100A9') 
fetal = c('FABP1', 'FGA', 'FGB', 'ITIH2', 'KNG1', 'OTC', 'SLC38A4')

plot_panel(all_patient_avg_ga, placental, isGA = TRUE, 'fig_2/a/fig_2a_placental.pdf')
plot_panel(all_patient_avg_ga, c(immune, fetal), isGA = TRUE, 'fig_2/a/fig_2a_immune_fetal.pdf')
