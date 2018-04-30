#Resources for installing edgeR
source("http://bioconductor.org/biocLite.R")
biocLite('edgeR')
edgeRUsersGuide ()

#Running samples from here
library(edgeR)
library(ggplot2)
library(colorspace)
library(limma)

#Set your working directory and reading in files
setwd("~/Desktop/CodingDesktop/whole root/")
x <- read.csv("WRraw_ALL_EdgeR.csv", header = TRUE, row.names=1)
x2 <- x[1:12] #experiment of interest is only these columns

#Make groups and subset data of interest.
groups=as.factor(c(rep("CM CIF+",3),rep("CM CIF-",3),rep("pWS CIF+",3),rep("pWS CIF-",3)))

#save gene names
genes=rownames(x)

#make DGElist
x <- DGEList(counts = x2, group=groups, genes = genes)

#make cpm and lcpm
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

#keep only genes that are expressed
dim(x)
keep.exprs <- rowSums(cpm>1)>=7
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x) #compare to dim(x) above

lcpm <- cpm(x, log=TRUE)
plot(density(lcpm))

#normalize data after removing low expressed genes
x <- calcNormFactors(x)

#cpm, lcpm of normalized values
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
write.csv(cpm,"Norm_CPM_WR.csv")

#account for factors in experiment
genotype=as.factor(c(rep("CM",3),rep("CM",3),rep("pWS",3),rep("pWS",3)))
peptide=as.factor(c(rep("CIF",3),rep("N",3),rep("CIF",3),rep("N",3)))

#assign factors to DGElist
x$samples$genotype=genotype
x$samples$peptide=peptide

#make design matrix of experimental factors
TS <- paste(x$samples$genotype, x$samples$peptide, sep=".")
x2=x
TS <- factor(TS, levels=c("CM.N","pWS.N","CM.CIF","pWS.CIF"))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
design #check matrix

#fit normalized dge list to model matrix
x2 <- estimateDisp(x2, design)
fit <- glmQLFit(x2, design)

#making contrast matrix for tests of interest
my.contrasts <- makeContrasts(CM.NvsCM.CIF=CM.CIF-CM.N, 
                              pWS.NvspWS.CIF=pWS.CIF-pWS.N, 
                              CM.NvspWS.N=pWS.N-CM.N, 
                              CM.CIFvspWS.CIF=pWS.CIF-CM.CIF,
                              eff=(pWS.CIF-pWS.N)-(CM.CIF-CM.N),
                              levels=design)

#QLF tests 
qlf_wtCIF <- glmQLFTest(fit, contrast=my.contrasts[,"CM.NvsCM.CIF"]) 
qlf_pWS_CIF <- glmQLFTest(fit, contrast=my.contrasts[,"pWS.NvspWS.CIF"]) 
qlf_pWS_wt <- glmQLFTest(fit, contrast=my.contrasts[,"CM.NvspWS.N"]) 
qlf_pWS_wt_CIF <- glmQLFTest(fit, contrast=my.contrasts[,"CM.CIFvspWS.CIF"])
qlf_eff <- glmQLFTest(fit, contrast=my.contrasts[,"eff"])


#how many genes are up and down in each?
summary(decideTestsDGE(qlf_wtCIF))
summary(decideTestsDGE(qlf_pWS_CIF))
summary(decideTestsDGE(qlf_pWS_wt))
summary(decideTestsDGE(qlf_pWS_wt_CIF))
summary(decideTestsDGE(qlf_eff))


#extract results and filter. 
#I filtered based on FDR<0.05 and logFC of 0.3 (about 2X fold change)
#toptags also adds F and FDR column.

#CM.NvsCM.CIF=CM.CIF-CM.N
wtCIF <- topTags(qlf_wtCIF, n = "Inf")$table
wtCIF <- subset(wtCIF, wtCIF$FDR<0.05)
wtCIF <- subset(wtCIF, abs(wtCIF$logFC)>0.3)
write.csv(wtCIF, "wtCIF_WR.csv")

#pWS.NvspWS.CIF=pWS.CIF-pWS.N
pWSCIF <- topTags(qlf_pWS_CIF, n = "Inf")$table
pWSCIF <- subset(pWSCIF, pWSCIF$FDR<0.05)
pWSCIF <- subset(pWSCIF, abs(pWSCIF$logFC)>0.3)
write.csv(pWSCIF, "pWSCIF_WR.csv")

#CM.NvspWS.N=pWS.N-CM.N
pWS_wt <- topTags(qlf_pWS_wt, n = "Inf")$table
pWS_wt <- subset(pWS_wt, pWS_wt$FDR<0.05)
pWS_wt <- subset(pWS_wt, abs(pWS_wt$logFC)>0.3)
write.csv(pWS_wt, "pWS_wt_WR.csv")

#CM.CIFvspWS.CIF=pWS.CIF-CM.CIF
pWS_wt_CIF <- topTags(qlf_pWS_wt_CIF, n="Inf")$table
pWS_wt_CIF <- subset(pWS_wt_CIF, pWS_wt_CIF$FDR<0.05)
pWS_wt_CIF <- subset(pWS_wt_CIF, abs(pWS_wt_CIF$logFC)>0.3)
write.csv(pWS_wt_CIF, "pWS_wt_CIF_WR.csv")


#eff=(pWS.CIF-pWS.N)-(CM.CIF-CM.N)
effCIF <- topTags(qlf_eff, n = "Inf")$table
effCIF <- subset(effCIF, effCIF$FDR<0.05)
effCIF <- subset(effCIF, abs(effCIF$logFC)>0.3)
write.csv(effCIF, "effCIF_WR.csv")

#######
######
##### Heatmaps 

##how are differentially expressed genes from tests above expressed in RootMap 3.0?

#Call in RootMap 3

FPKM <- read.csv("RM3_FPKM.csv", header=TRUE, row.names=1)
RM3 <- data.frame(FPKM)

#pWS+CIF vs pWS no CIF
pWSCIF_up = subset(pWSCIF, pWSCIF$logFC>0)
pWSCIF_dn = subset(pWSCIF, pWSCIF$logFC<0)

genes_up = rownames(pWSCIF_up)
genes_dn = rownames(pWSCIF_dn)

exps_up = RM3[row.names(RM3) %in% genes_up,]
exps_dn = RM3[row.names(RM3) %in% genes_dn,]

#heatmap up
pal <- choose_palette()
coul = pal(20)
dat <- as.matrix(exps_up)
heatmap(dat, col = coul)

#heatmap down
coul = pal(20)
dat <- as.matrix(exps_dn)
heatmap(dat, col = coul)

#heatmap of the positive logFC genes from the eff test

#eff (pWS.CIF-pWS.N)-(WT.CIF-WT.N)
eff_up = subset(effCIF, effCIF$logFC>0)
eff_dn = subset(effCIF, effCIF$logFC<0)

genes_up = rownames(eff_up)
genes_dn = rownames(eff_dn)

exps_up = RM3[row.names(RM3) %in% genes_up,]
exps_dn = RM3[row.names(RM3) %in% genes_dn,]

#heatmap up
coul = pal(10)
dat <- as.matrix(exps_up)
heatmap(dat, col = coul)

#heatmap down
coul = pal(10)
dat <- as.matrix(exps_dn)
heatmap(dat, col = coul, cexRow=0.5, cexCol=0.5)

####
#### also noticed root hair-related genes are down in pWS vs WT whole root

pWSvWT_dn = subset(pWS_wt, pWS_wt$logFC<0)
genes_dn = rownames(pWSvWT_dn)
exps_dn = RM3[row.names(RM3) %in% genes_dn,]

coul = pal(20)
dat <- as.matrix(exps_dn)
heatmap(dat, col = coul)
