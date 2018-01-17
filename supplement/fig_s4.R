library(dplyr)
library(ggplot2)
library(reshape2)

#Set working directory
mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

#Source scripts
source('plot_theme.R')
source('load_data.R')

#Source data
load('counts_gene_lists.RData')
cpm_rnaseq = 'fig_3/a/data_combined_edgeR.csv'

#Load qPCR data
counts_all = load_process_qPCR(counts_per_mL_all, genes_6922C, is_classification = FALSE, filter_5695 = FALSE)
counts_all$cohort = factor(counts_all$cohort)

counts_penn = counts_all %>% filter(cohort == 'PT_Penn'|cohort == 'T_Penn') %>% select(-starts_with('dCD'), -starts_with('patient'), -starts_with('cohort'), -starts_with('ga'),
                                                                                       -starts_with('delivery'))

#Load RNASeq Data
rnaseq = read.csv(cpm_rnaseq, header=TRUE, strip.white = TRUE)

#Get sample_ids in common
idx1 = match(counts_penn$sample_id, colnames(rnaseq))
rnaseq = rnaseq[,c(1,idx1)]

#Get genes in common
rnaseq = rnaseq[rnaseq$genes %in% colnames(counts_penn),]
idx2 = match(rnaseq$genes, colnames(counts_penn))
counts_penn = counts_penn[,c(1,idx2)]

#Manipulate data to plot
rownames(rnaseq) = rnaseq$external_gene_name
rnaseq = as.data.frame(t(rnaseq[,c(-1)]))
rnaseq$sample_id = rownames(rnaseq)

rnaseq_long = melt(rnaseq, id = c('sample_id'))
counts_long = melt(counts_penn, id = c('sample_id'))
rnaseq_long$count = counts_long$value

#Get correlation
corr = cor.test(log2(rnaseq_long$value+1), log2(rnaseq_long$count+1))

#Plot data
plot = ggplot(rnaseq_long) +
  geom_point(aes(x = log2(value+1), y = log2(count+1)),  size =1, color = 'black', alpha = 0.8) +
  theme_pub() +
  labs(x = '', y = '') 

pdf('supplement/s4/fig_s4.pdf', useDingbats = FALSE, width = 3, height = 3)
print(plot)
dev.off()