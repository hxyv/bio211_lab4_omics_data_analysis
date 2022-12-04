# Prepare R packages
chooseCRANmirror()

install.packages('stringr')
install.packages('ggplotify')
install.packages('patchwork')
install.packages('cowplot')

install.packages('BiocManager')
BiocManager::install(c('edgeR', 'limma'))
# Install edgeR and limma packages for Bioconductor

# Set work path, read data and data statistics
path <- 'C:/Users/xingyu hu/Documents/BIO211_Lab4_Omics_data_analysis'
setwd(paste(path))
# Change the work path accordingly
# For example, path <- 'C://BIO211/LAB4'

rm(list=ls()) # Clean work environment
a <- read.csv('Data/LIHC_counts_lab4.csv')
dim(a)

row.names(a) <- a[,1]
a <- a[, -1] # Delete first column
group_list <- ifelse(as.numeric(substr(colnames(a), 14, 15))<10, 'tumor', 'normal')
# Separate samples into two groups
# Check the hints on TCGA Barcode for how to differentiate tumor/normal samples

group_list <- factor(group_list, levels=c("normal", "tumor"))
table(group_list)

View(a)
# Check the data matrix

# Select significant genes exhibiting differential expression
## Method 1: edgeR
library(edgeR)
dge <- DGEList(counts=a, group=group_list)
dge$samples$lib.size<-colSums(dge$counts)
design<-model.matrix(~0+group_list)
rownames(design)<-colnames(dge)
colnames(design)<-levels(group_list)
dge<-estimateGLMCommonDisp(dge,design)
dge<-estimateGLMTrendedDisp(dge,design)
dge<-estimateGLMTagwiseDisp(dge,design)
fit<-glmFit(dge,design)
fit2<-glmLRT(fit,contrast=c(-1,1))
DEG2=topTags(fit2, n=nrow(a))
DEG2=as.data.frame(DEG2)
logFC_cutoff2<-with(DEG2, mean(abs(logFC))+2*sd(abs(logFC)))
DEG2$change=as.factor(ifelse(DEG2$PValue<0.05 & abs(DEG2$logFC) > logFC_cutoff2, 
                             ifelse(DEG2$logFC > logFC_cutoff2,"UP","DOWN"),
                             "NOT"))
head(DEG2)
table(DEG2$change)

# Method 2: limma-voom
library(limma)
design<-model.matrix(~0+group_list)
colnames(design)=levels(group_list)
rownames(design)=colnames(a)
dge<-DGEList(counts=a)
dge<-calcNormFactors(dge)
logCPM<-cpm(dge,log=TRUE,prior.count=3)
v<-voom(dge,design,normalize="quantile")
fit<-lmFit(v,design)
constra=paste(rev(levels(group_list)),collapse="-")
cont.matrix<-makeContrasts(contrasts=constra,levels=design)
fit3=contrasts.fit(fit,cont.matrix)
fit3=eBayes(fit3)
DEG3=topTable(fit3,coef=constra,n=Inf)
DEG3=na.omit(DEG3)
logFC_cutoff3<-with(DEG3,mean(abs(logFC)+2*sd(abs(logFC))))
DEG3$change=as.factor(ifelse(DEG3$P.Value<0.05 & abs(DEG3$logFC) > logFC_cutoff3,
                             ifelse(DEG3$logFC > logFC_cutoff3,"UP","DOWN"),
                             "NOT"))
head(DEG3)
table(DEG3$change)
# Count the number of up- and down-regulated genes

# Save your results
edgeR_DEG <- DEG2
limma_voom_DEG <- DEG3
save(edgeR_DEG, limma_voom_DEG, group_list, file='DEG.Rdata')
