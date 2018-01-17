library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)

#Set working directory
mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

#Source scripts
source('plot_theme.R')
source('load_data.R')

#Load data
load('counts_gene_lists.RData')

#Load data
counts_all = load_process_qPCR(counts_per_mL_all, genes_6922AC, is_classification = FALSE, filter_5695 = FALSE)
counts_all$cohort = factor(counts_all$cohort)
counts_all$dCD = round(counts_all$dCD)

#Separate cohorts
danish = counts_all %>% filter(cohort == "Danish") %>% select(sample_id, patient, ga, delivery)
penn = counts_all %>% filter(cohort == "PT_Penn" | cohort == 'T_Penn') %>% select(sample_id, ga, delivery)
penn = rbind(penn, data.frame(sample_id = c('T_Penn_16', 'T_Penn_27'), ga = c(26,29), delivery = c(40,40))) #Add two samples that are part of RNASeq but not qPCR
uab = counts_all %>% filter(cohort == 'PT_UAB' | cohort == 'T_UAB') %>% select(sample_id, ga, delivery)

#Select unique patients
patients = unique(danish %>% select(patient))
danish$patient = reorder(danish$patient, -danish$delivery)

#Order Data
for (p in patients[,1]) {
  danish$patID[danish$patient==p] = match(p, levels(danish$patient))
}

penn$patient = rownames(penn)
penn$patient = reorder(penn$patient, -penn$delivery)

for (row in penn$patient) {
  penn$patID[penn$patient==row] = match(row, levels(penn$patient))
}

uab$patient = rownames(uab)
uab$patient = reorder(uab$patient, -uab$delivery)

for (row in uab$patient) {
  uab$patID[uab$patient==row] = match(row, levels(uab$patient))
}

unique_patID = data.frame(patID = unique(danish$patID))

#Plot
plot_sampling_danish =  ggplot(danish) +
  geom_rect(aes(xmin = 0, xmax = 48, ymin = 0, ymax = 32), fill = '#a1d99b') +
  geom_segment(data = unique_patID, aes(x = 0, xend = 48, y = patID, yend = patID), color = 'black') +
  geom_segment(aes(x = 12, y = 0, xend = 12, yend = 32), linetype = 'dashed', color = '#08519c', size = 1.3) +
  geom_segment(aes(x = 24, y = 0, xend = 24, yend = 32), linetype = 'dashed', color = '#08519c', size = 1.3) +
  geom_segment(aes(x = 37, y = 0, xend = 37, yend = 32), linetype = 'dashed', color = '#08519c', size = 1.3) +
  geom_point(aes(x=delivery, y=patID+0.4), size = 6, shape = 25, color = 'black', fill = 'black') +
  geom_point(aes(x=ga, y=patID), fill = 'black', color = 'black', size = 3, shape = 15) +
  geom_text(aes(label = '0', x=0, y=-0.8), size = 8) +
  geom_text(aes(label = '12', x=12, y=-0.8), size = 8) +
  geom_text(aes(label = '24', x=24, y=-0.8), size = 8) +
  geom_text(aes(label = '37', x=37, y=-0.8), size = 8) +
  geom_text(aes(label = 'High resolution \n Denmark', x=-3, y=16), size = 8, angle = 90) +
  geom_text(aes(label = 'Gestational age (weeks)', x = 22.5, y = -3), size = 8) +
  coord_cartesian(xlim = c(-3,48), ylim = c(-3,32)) +
  theme_pub() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  labs(x='')

penn$patID = penn$patID+25

penn_uab = rbind(uab, penn)

plot_sampling_penn_uab =  ggplot(penn_uab) +
  geom_rect(aes(xmin = 0, xmax = 48, ymin = 25, ymax = 41), fill = '#6baed6') +
  geom_rect(aes(xmin = 0, xmax = 48, ymin = 0, ymax = 24), fill = '#3182bd') +
  geom_segment(aes(x = 0, xend = 48, y = penn_uab$patID, yend = penn_uab$patID), color = 'black') +
  geom_segment(aes(x = 12, y = 0, xend = 12, yend = 24), linetype = 'dashed', color = '#08519c', size = 1.3) +
  geom_segment(aes(x = 12, y = 25, xend = 12, yend = 41), linetype = 'dashed', color = '#08519c', size = 1.3) +
  geom_segment(aes(x = 24, y = 0, xend = 24, yend = 24), linetype = 'dashed', color = '#08519c', size = 1.3) +
  geom_segment(aes(x = 24, y = 25, xend = 24, yend = 41), linetype = 'dashed', color = '#08519c', size = 1.3) +
  geom_segment(aes(x = 37, y = 0, xend = 37, yend = 24), linetype = 'dashed', color = '#08519c', size = 1.3) +
  geom_segment(aes(x = 37, y = 25  , xend = 37, yend = 41), linetype = 'dashed', color = '#08519c', size = 1.3) +
  geom_point(aes(x=delivery, y=patID+0.4), size = 6, shape = 25, color = 'black', fill = 'black') +
  geom_point(aes(x=ga, y=patID), fill = 'black', color = 'black', size = 3, shape = 15) +
  geom_text(aes(label = '0', x=0, y=-0.8), size = 8) +
  geom_text(aes(label = '12', x=12, y=-0.8), size = 8) +
  geom_text(aes(label = '24', x=24, y=-0.8), size = 8) +
  geom_text(aes(label = '37', x=37, y=-0.8), size = 8) +
  geom_text(aes(label = 'Preterm 1 \n Pennsylvania', x=-3, y=35.5), size = 8, angle = 90) +
  geom_text(aes(label = 'Preterm 2 \n Alabama', x=-3, y=13.5), size = 8, angle = 90) +
  geom_text(aes(label = 'Gestational age (weeks)', x = 22.5, y = -3), size = 8) +
  coord_cartesian(xlim = c(-3,48), ylim = c(-3,43)) +
  theme_pub() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  labs(x='')

grid = plot_grid(plot_sampling_danish, plot_sampling_penn_uab, align = 'h')
pdf(file = 'fig_1/fig_1.pdf', width = 24, height = 15, useDingbats = FALSE)
print(grid)
dev.off()