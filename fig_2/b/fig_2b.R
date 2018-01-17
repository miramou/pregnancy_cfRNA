library(corrplot)
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

#Preprocess data
counts_all = load_process_qPCR(counts_per_mL_all_cma, genes_5695_6922C, is_classification = FALSE, filter_5695 = FALSE)
counts_danish = counts_all %>% filter(cohort == 'Danish') %>% select(-starts_with('panel'), -starts_with('cohort'), -starts_with('sample_id'), -starts_with('ga'), -starts_with('patient'))

counts_danish = counts_all %>% filter(cohort == 'Danish') %>% select(-starts_with('panel'), -starts_with('cohort'), -starts_with('sample_id'), -starts_with('ga'), -starts_with('patient'))
counts_danish = counts_danish %>% select(-starts_with('sample_id'), -starts_with('dCD'))

#Gene list
all_list = c('ALPP', 'CAPN6','CGA', 'CGB','CSHL1', 'LGALS14', 'PAPPA', 'PLAC4', 'PSG7', 'ADAM12','CSH1', 'CSH2', 'GH2', 'HSD17B1', 
             'HSD3B1','HSPB8', 'PSG1', 'PSG4','VGLL1', 'ANXA3', 'ARG1', 'CAMP', 'DEFA4', 'LTF', 'MMP8', 'PGLYRP1', 
             'S100A8', 'S100A9', 'S100P', 'FABP1', 'FGA', 'FGB', 'ITIH2', 'KNG1', 'OTC', 'SLC38A4')

placental = c('ALPP', 'CAPN6','CGA', 'CGB','CSHL1', 'LGALS14', 'PAPPA', 'PLAC4', 'PSG7', 'ADAM12','CSH1', 'CSH2', 'GH2', 'HSD17B1', 
              'HSD3B1','HSPB8', 'PSG1', 'PSG4','VGLL1')
immune = c('ANXA3', 'ARG1', 'CAMP', 'DEFA4', 'LTF', 'MMP8', 'PGLYRP1','S100A8', 'S100A9', 'S100P')
fetal = c('FABP1', 'FGA', 'FGB', 'ITIH2', 'KNG1', 'OTC', 'SLC38A4')

all = c(match(all_list, colnames(counts_danish)))

#Calc Corr
counts_to_plot = counts_danish[,all]
M = cor(counts_to_plot, use = 'complete.obs')
diag(M) = 0 #Mask diagonal

corr_df = as.data.frame(M)
corr_df$gene2 = rownames(corr_df)

corr_melted=  melt(corr_df, id = c('gene2'))
corr_melted$variable = factor(corr_melted$variable, levels = rev(colnames(M)))
corr_melted$gene2 = factor(corr_melted$gene2, levels = rownames(M))

#Isolate one triangle
ind = which(upper.tri(M,diag=F), arr.ind=TRUE)
corr_tri = data.frame(col=dimnames(M)[[2]][ind[,2]], row = dimnames(M)[[1]][ind[,1]], val = M[ind])

#Label samples with gene group of two genes that were used for corr in each cell
#For each gene: 0 = placental, 2 = immune, 3 = fetal
#Therefore the sum: 0 = placental/placental, 4 = immune/immune, 6 = fetal/fetal, 3 = placental, fetal

id = c()
for (idx in 1:dim(corr_tri)[1]) {
  temp = 0
  
  if (as.character(corr_tri$col[idx]) %in% immune) {
    temp = 2
  } else if (as.character(corr_tri$col[idx]) %in% fetal) {
    temp = 3
  }
  
  if (as.character(corr_tri$row[idx]) %in% immune) {
    temp = temp + 2
  } else if (as.character(corr_tri$row[idx]) %in% fetal) {
    temp = temp + 3
  }
  
  id = c(id, temp)
}
corr_tri$id = id

#Calc median corr in each gene group (placental, immune, fetal, placental-fetal)
median = c(median(as.matrix(corr_tri %>% filter(id == 0) %>% select(val))), median(as.matrix(corr_tri %>% filter(id == 4) %>% select(val))), 
           median(as.matrix(corr_tri %>% filter(id == 6) %>% select(val))), median(as.matrix(corr_tri %>% filter(id == 3) %>% select(val)))) #Placental, Immune, Fetal, Placental-Fetal
#Perform one sample t-test
pvals = c(t.test(as.matrix(corr_tri %>% filter(id == 0) %>% select(val)), mu=0)$p.value, t.test(as.matrix(corr_tri %>% filter(id == 4) %>% select(val)), mu=0)$p.value, 
          t.test(as.matrix(corr_tri %>% filter(id == 6) %>% select(val)), mu=0)$p.value, t.test(as.matrix(corr_tri %>% filter(id == 3) %>% select(val)), mu=0)$p.value) #Placental, Immune, Fetal, Placental-Fetal#Plot
pvals_adj = p.adjust(pvals, method = 'bonferroni')

#Plot
corr_plot = ggplot(corr_melted) +
  geom_tile(aes(x= gene2, y = variable, fill = value), color = 'white') +
  scale_fill_gradient2(low = '#ca0020', mid = '#f7f7f7', midpoint = 0, high = '#0571b0', space = "Lab", guide = 'colourbar', limits = c(-1,1)) +
  theme_pub() +
  theme(
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9, angle = 90),
    legend.position = 'right'
  ) +
  labs(x = '', y = '', fill = '') +
  scale_x_discrete(position = 'top') +
  guides(fill = guide_colorbar(barwidth = 0.5, barheight = 7, ticks = FALSE))

ggsave('fig_2/b/fig_2b.eps', width = 5.2, height = 4.2)
