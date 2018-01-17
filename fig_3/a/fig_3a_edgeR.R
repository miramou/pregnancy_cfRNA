library('edgeR')
library('limma')
library('locfit')
library('splines')
library('KernSmooth')
library('statmod')
library('MASS')
library('dplyr')

#Set working directory
mainDir = '' #CHANGE this to reflect the main directory you've selected
setwd(mainDir)

#Get raw data
rawdata = read.csv("raw_data/rnaseq_counts.csv")
group = factor(c(rep(1, 7), rep(2,8)))  #1 = Term, 2 = Preterm

y = DGEList(counts=rawdata[,2:ncol(rawdata)], genes=rawdata[,1],group=group)
keep = rowSums(cpm(y)>1) >= 2
y = y[keep,]

y$samples$lib.size = colSums(y$counts)

delivery = group
data.frame(Sample=colnames(y),delivery)

design = model.matrix(~0+delivery)
rownames(design) = colnames(y)

y = estimateGLMCommonDisp(y, design, verbose=TRUE)
y = estimateGLMTrendedDisp(y,design)
y = estimateGLMTagwiseDisp(y,design)

#Exact Test
et = exactTest(y)
et_tags = topTags(et, n=200)$table %>% filter(PValue < 0.001)

#GLM
fit = glmFit(y, design)
contrast.design= c(1,-1)
lrt = glmLRT(fit, contrast=contrast.design)
glm_tags = topTags(lrt, n = 200)$table %>% filter(PValue < 0.001)

#QL F-test
fit = glmQLFit(y, design)
contrast.design=c(1,-1)
qlf = glmQLFTest(fit, contrast=contrast.design)
qlf_tags = topTags(qlf, n=200)$table  %>% filter(PValue < 0.001)

#Merge tests
top_genes = Reduce(intersect, list(et_tags$genes, glm_tags$genes, qlf_tags$genes))

#Normalize data
rawdata_norm = as.data.frame(t(t(rawdata[,-1])/colSums(rawdata[,-1])*1000000))

#Filter data and save only genes of interest
idx = match(top_genes, rawdata$external_gene_name)
data_filtered= cbind(rawdata$external_gene_name[idx], rawdata_norm[idx,])
names(data_filtered)[names(data_filtered)=='rawdata$external_gene_name[idx]'] = 'genes'
data_filtered = data_filtered[-match('Y_RNA', data_filtered$genes),] #artifact of all female cohort
write.csv(data_filtered, 'fig_3/a/data_combined_edgeR.csv', row.names = FALSE)