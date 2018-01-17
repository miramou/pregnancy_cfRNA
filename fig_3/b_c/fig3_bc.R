library(ggplot2)
library(reshape2)
library(exactRankTests)
library(pROC)
library(effsize)

#Set working directory
mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

#Source scripts
source('plot_theme.R')
source('load_data.R')
source('fig_3/b_c/calc_mom.R')
source('fig_3/b_c/fig3_bc_helper.R')

#Source data
load('counts_gene_lists.RData')

#Process data
counts_all = load_process_qPCR(counts_per_mL_all_cma, genes_6922AC, is_classification = FALSE, filter_5695 = TRUE)

#Calc mom
mom_counts = calc_mom_avg_wrapper(counts_all)
mom_counts$Class = ifelse(mom_counts$cohort == 'PT_UAB' | mom_counts$cohort == 'PT_Penn', 'pre', 'term')
mom_counts$Class = factor(mom_counts$Class)

#Split Penn and UAB
penn = mom_counts %>% filter(cohort == 'PT_Penn' | cohort == 'T_Penn' | cohort == 'Danish') 
uab = mom_counts %>% filter(cohort == 'PT_UAB' | cohort == 'T_UAB'| cohort == 'Danish')

#Perform fisher exact test and hedge's g calculation for every gene

#Penn statistic calculation
penn_long = melt(penn, id= c('sample_id', 'cohort', 'ga', 'delivery', 'Class'))
penn_stats = calc_eff_size_test(penn_long)
genes_penn = penn_stats %>% filter(mag == 'large', p_val_adj <= 0.05) %>% select(gene, estimate, p_val_adj)

#Reduce datasets to those genes identified in discovery dataset
idx = match(genes_penn[,1], colnames(penn))
penn = penn[,c(1:4, dim(penn)[2], idx)]
idx = match(genes_penn[,1], colnames(uab))
uab = uab[,c(1:4, dim(uab)[2], idx)]

#UAB statistic calculation
uab_long = melt(uab, id= c('sample_id', 'cohort', 'ga', 'delivery', 'Class'))
uab_stats = calc_eff_size_test(uab_long)
genes_uab = uab_stats %>% filter(mag == 'large' , p_val_adj <= 0.05) %>% select(gene, estimate, p_val_adj)

#Identify predictive combinations
combos = get_combos(genes_penn, penn)
gselect = combos %>% filter(tp>=6, combo_num == 3, fp<=1)
results = get_results(gselect, genes_penn, penn)

#Test in discovery set
penn_pred = test_results(gselect, genes_penn, penn)

#Test in validation set
uab_for_roc = uab %>% filter(cohort != 'Danish') #remove Danish samples used during discovery
uab_pred = test_results(gselect, genes_penn, uab_for_roc)

#Build roc curves
roc_penn = roc(penn_pred$Class, penn_pred$prob, levels = c('term', 'pre'), ci = TRUE, ci.method = 'bootstrap', auc=TRUE)
roc_uab = roc(uab_pred$Class, uab_pred$prob, levels = c('term', 'pre'), ci = TRUE, ci.method = 'bootstrap', auc=TRUE)

#Plot genes of interest
to_plot = unique(c(as.character(results[,1]), as.character(results[,2]), as.character(results[,3])))
all = rbind(penn, uab_for_roc)
plot_PT_T(all, genesToPlot = to_plot, 'fig_3/b_c/b_gene_plots/')

#Plot roc curve
plot_roc(roc_penn, roc_uab, 'fig_3/b_c/fig_3c.pdf')
