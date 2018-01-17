#Load all csvs and merge into one Rdata object with all important info for later use
require(plyr)

#Set working directory
mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

#Source Raw qPCR Data
Danish_Panel5695 = 'raw_data/Denmark_5695.csv'
Danish_Panel6922 = 'raw_data/Denmark_6922C.csv'
Penn_UAB_Panel6922 = 'raw_data/Upenn_UAB_6922AC.csv'
ntc_Panel6922 = 'raw_data/NTC_6922C.csv'
ga_collection_delivery_info = 'raw_data/UPenn_UAB_Collection_Delivery_Dates.csv'
plasma_info = 'raw_data/Plasma_Input.csv'

#Load relevant functions
source('data_consolidation/calc_counts_per_mL.R')
source('data_consolidation/calc_centered_moving_avg.R')
source('data_consolidation/combine_ga_dCD.R')

#Calculate counts per mL
counts_per_mL_Danish_5695 = calc_copy_num_denoised_wrapper(Danish_Panel5695, ntc_Panel6922, 'Danish', sweepStart = 6, plasma_info)
counts_per_mL_Danish_6922 = calc_copy_num_denoised_wrapper(Danish_Panel6922, ntc_Panel6922, 'Danish', sweepStart = 6, plasma_info)
counts_per_mL_Penn_UAB = calc_copy_num_denoised_wrapper(Penn_UAB_Panel6922, ntc_Panel6922, 'UPenn_UAB', sweepStart = 6, plasma_info)

#Label cohort and panel
counts_per_mL_Danish_5695$panel = '5695'
counts_per_mL_Danish_5695$cohort = 'Danish'
counts_per_mL_Danish_6922$panel = '6922C'
counts_per_mL_Danish_6922$cohort = 'Danish'
counts_per_mL_Penn_UAB$panel = '6922AC'
counts_per_mL_Penn_UAB$cohort = gsub("\\_[0-9]*$", "", counts_per_mL_Penn_UAB$sample_id)

#Put Danish all together
counts_Danish = join(counts_per_mL_Danish_5695, counts_per_mL_Danish_6922, by = 'sample_id', type = 'full')

#Put all data together
counts_per_mL_all = join(counts_Danish %>% select(-starts_with('gestation_age'), -starts_with('time_to_delivery')),counts_per_mL_Penn_UAB, by = 'sample_id', type = 'full')

#Combine dCD and gestation age data
counts_per_mL_all = combineGA_dCD(ga_collection_delivery_info, counts_per_mL_all, isDanish = TRUE, counts_Danish)

#Remove any samples collected more than 9 weeks in advance of delivery (unlikely to have detectable signal)
counts_per_mL_all = counts_per_mL_all %>% filter(cohort == 'PT_Penn' & dCD >= -9 | 
                                                   cohort == 'PT_UAB' & dCD >= -9 | 
                                                   cohort == 'T_Penn' | 
                                                   cohort == 'T_UAB' | 
                                                   cohort == 'Danish', 
                                                 sample_id != 'PT_Penn_12') #pre-eclampsia

#Centered moving average for Danish cohort
count_cma_Danish = cma(counts_per_mL_all[counts_per_mL_all$cohort == 'Danish',]) %>% ungroup('patient') 
counts_Danish_pp = counts_per_mL_all %>% filter(cohort == 'Danish', dCD > 0)
counts_per_mL_all_cma = rbind(counts_per_mL_all[counts_per_mL_all$cohort != 'Danish',] , count_cma_Danish)
counts_per_mL_all_cma = rbind(counts_per_mL_all_cma, counts_Danish_pp)

genes_6922C = setdiff(colnames(counts_per_mL_Danish_6922), colnames(counts_per_mL_Danish_5695))
genes_5695 = setdiff(colnames(counts_per_mL_Danish_5695), colnames(counts_per_mL_Danish_6922))
genes_6922AC = Reduce(intersect, list(colnames(counts_per_mL_Danish_6922), colnames(counts_per_mL_Penn_UAB)))
genes_5695_6922AC = Reduce(intersect, list(colnames(counts_per_mL_Danish_6922), colnames(counts_per_mL_Danish_5695), colnames(counts_per_mL_Penn_UAB)))
genes_5695_6922C = Reduce(intersect, list(colnames(counts_per_mL_Danish_6922), colnames(counts_per_mL_Danish_5695)))

#Remove missing data
counts_per_mL_all_cma = counts_per_mL_all_cma %>% filter(dCD != 'NA')

#Save Rdata objects
save(counts_per_mL_all, counts_per_mL_all_cma, genes_6922C, genes_5695, genes_6922AC, genes_5695_6922AC, genes_5695_6922C, file = 'counts_gene_lists.RData')
