#Resources for installing edgeR
source("http://bioconductor.org/biocLite.R")
biocLite('edgeR')
edgeRUsersGuide ()

#Running samples from here
library(edgeR)
library(ggplot2)
library(colorspace)
library(limma)

setwd("~/Desktop/")

#Set your working directory and read in files
sorted <- read.csv("RawCounts_FACs.csv", header = TRUE, row.names=1)

#do root map processing
rootmap_raw <- read.csv("CPM_MY.csv", header = TRUE, row.names=1)

#average reps of RootMap 3
rootmap_av= t(apply(rootmap_raw, 1, tapply, gl(22, 3), mean))

#name columns of reps
colnames(rootmap_av)= colnames(rootmap_raw)[seq(1, 66, 3)]


#We will merge the first two files above.
#first need toAverage counts for each rep first.
CASP_CIF=data.frame(rowSums(sorted[,1:3])/3)
CASP_none=data.frame(rowSums(sorted[,4:6])/3)
pWS_CIF=data.frame(rowSums(sorted[,7:9])/3)
pWS_none=data.frame(rowSums(sorted[,10:12])/3)

CASP_CIF=CASP_CIF/sum(CASP_CIF)*1000000
CASP_none=CASP_none/sum(CASP_none)*1000000
pWS_CIF=pWS_CIF/sum(pWS_CIF)*1000000
pWS_none=pWS_none/sum(pWS_none)*1000000

colnames(CASP_CIF)="CASP_CIF"
colnames(CASP_none)="CASP_none"
colnames(pWS_CIF)="pWS_CIF"
colnames(pWS_none)="pWS_none"

all_conds=cbind(CASP_CIF, CASP_none, pWS_CIF, pWS_none)

#finally, merge
merged=merge(rootmap_av, all_conds, by = "row.names")


#make dge
comps=c(1,3,4,5,6,8,12,13,14,15,16,18,19,20,21,23, 24, 25, 26)
merged_dge=DGEList(counts=merged[,comps+1], genes=merged[,1])


#calc norm factor
merged_dge=calcNormFactors(merged_dge)


#plotMDS with log. Maybe use this one.

x<-plotMDS(cpm(merged_dge, log=T), gene.selection="pairwise", pch=21, xlim=c(-4,6), ylim=c(-4,4))
with(merged_dge, text(x,labels=rownames(merged_dge$samples), pos=3, cex=0.3))




