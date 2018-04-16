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
downsample_no_cma = rbind(downsample_no_cma, counts_danish_no_cma[counts_danish_no_cma$patient %in% patients[11:31,1],])

for (p in patients[1:10, 1]) {
  downsample = seq(2, dim(counts_danish_no_cma %>% filter(patient == p))[1], by = 2)
  patient_downsample = counts_danish_no_cma %>% filter(patient == p)
  patient_downsample = patient_downsample[downsample,]
  downsample_no_cma = rbind(downsample_no_cma, patient_downsample)
}

all_patient_avg_ga = rbind(avg_wrapper(downsample_no_cma, interval = 4, isGA = TRUE), avg_wrapper(counts_21, interval = 4, isGA = TRUE))
all_patient_avg_dCD = rbind(avg_wrapper(downsample_no_cma, interval = 4, isGA = FALSE), avg_wrapper(counts_21, interval = 4, isGA = FALSE))

placental = c('PLAC1', 'PTGER3') 
immune = c('BPI', 'CD160', 'CD180', 'CD2', 'CD5', 'CEACAM6', 'CNOT7', 'EPB42', 'HMGN2', 'KRT8') 

genes_all=rbind(placental, immune)

#Plot panel
plot_panel(all_patient_avg_ga, c(as.character(genes_all)), isGA = TRUE, 'supplement/s1/fig_s1.pdf')