#Running the edgeR R package on RNA-seq data from CASP1-mCherry sorted cells from CASP1-mCherry, CASP1-mCherry+CIF2, pWERSHRGFP/CASP1-mCherry and pWERSHRGFP/CASP1-mcherry+CIF2 samples. Samples were sequenced on the HiSeq4000, 50bp SR.
biocLite('edgeR')
library(edgeR)

#Set working directory
setwd("~/Desktop/CodingDesktop")

#From here is a modification of Julin Maloof's code from the CSHL plant course 2011

#Read in input counts
counts <- read.csv("EdgeR_raw_CED.csv", header = TRUE, row.names=1)
countsq <- read.csv("EdgeR_raw_RootMap_GT.csv", header = TRUE, row.names=1)

#check the file
head(counts)
summary(counts)
countsq <- na.omit(countsq)

#convert data to edgeR object list
#Name libraries by group and the # of reps
groups=rep(c("CASP+CIF", "CASP","pWS+CIF","pWS"), each=3)
groupsq=rep(c("315.1","E30"), each=3)
y <- DGEList(counts=counts,group=groups)
q <- DGEList(counts=countsq, group=groupsq)

#Export the normalized counts in counts per million (cpm)
normalized_counts <- cpm(q, normalized.lib.sizes=FALSE)
write.table(normalized_counts, "Normalized-CPM-MY.txt")

#normalize library 
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)


q <- calcNormFactors(q)
q <- estimateCommonDisp(q)


#Filter non-expressed genes
A <- aveLogCPM(y)
y2 <- y[A>2,]
head(y2$counts)

B <- aveLogCPM(q)
q2 <- q[B>2,]
head(q2$counts)

#compare how many genes riltered out
dim(y)
dim(y2)

dim(q)
dim(q2)

#calculate DE genes
#Tests for differences between two groups of negative binomially distributed counts

#PairA: CASP+CIF vs CASP addressing the question "What are the effects of CASP in wild-type endodermis?"
#PairB: CASP vs pWS addressing the question "What are the transcriptional changes as a result of expressing SHR in the epidermis?"
#PairC: CASP+CIF vs pWS+CIF addressing the question of "How does the presence of CIF alter the pair above?"
#PairD: pWS vs pWS+CIF addressing the question of "What is the effect of CIF2 on plants expressing SHR in the epidermis?"
#PairE: E30 vs 315.1 addressing "What are the typical differences between cortex and endodermis?"
#PairF: CASP+CIF vs pWS - for generating Venn Diagrams
#PairG: CASP vs pWS+CIF - for generating Venn Diagrams

DEtestA = exactTest(y2, pair=c("CASP","CASP+CIF"))
DEtestB = exactTest(y2, pair=c("CASP","pWS"))
DEtestC = exactTest(y2, pair=c("CASP+CIF","pWS+CIF"))
DEtestD = exactTest(y2, pair=c("pWS","pWS+CIF"))
DEtestE = exactTest(q2, pair=c("E30","315.1"))
DEtestF = exactTest(y2, pair=c("CASP+CIF","pWS"))
DEtestG = exactTest(y2, pair=c("CASP","pWS+CIF"))

#create a table of the results, with multiple testing correction
#Table of the Top Differentially Expressed Tags
resultsA = topTags(DEtestA,n=Inf)
resultsB = topTags(DEtestB,n=Inf)
resultsC = topTags(DEtestC,n=Inf)
resultsD = topTags(DEtestD,n=Inf)
resultsE = topTags(DEtestE,n=Inf)
resultsF = topTags(DEtestF,n=Inf)
resultsG = topTags(DEtestG,n=Inf)

#Export all results
resultsA = (resultsA$table)
resultsB = (resultsB$table)
resultsC = (resultsC$table)
resultsD = (resultsD$table)
resultsE = (resultsE$table)
resultsF = (resultsF$table)
resultsG = (resultsG$table)

write.table(resultsA, "CASP+CIF_vs_CASP.txt")
write.table(resultsB, "CASP_vs_pWS.txt")
write.table(resultsC, "CASP+CIF_vs_pWS+CIF.txt")
write.table(resultsD, "pWS vs pWS+CIF.txt")
write.table(resultsE, "E30 vs 315.1.txt")
write.table(resultsF, "CASP+CIF_vs_pWS.txt")
write.table(resultsG, "CASP_vs_pWS+CIF.txt")

#filteredresults
zA = data.frame(resultsA)
zB = data.frame(resultsB)
zC = data.frame(resultsC)
zD = data.frame(resultsD)
zE = data.frame(resultsE)
zF = data.frame(resultsF)
zG = data.frame(resultsG)

#filtered results for pWS+CIF vs pWS (less stringent)
filterD = subset(zD, zD$PValue<0.01)

#filtered results for all other pairs
filterA = subset(zA, zA$FDR<0.01)
filterA = subset(filterA, filterA$PValue<0.01)
filterA = subset(filterA, abs(filterA$logFC)>2)

filterB = subset(zB, zB$FDR<0.01)
filterB = subset(filterB, filterB$PValue<0.01)
filterB = subset(filterB, abs(filterB$logFC)>2)

filterC = subset(zC, zC$FDR<0.01)
filterC = subset(filterC, filterC$PValue<0.01)
filterC = subset(filterC, abs(filterC$logFC)>2)

filterE = subset(zE, zE$FDR<0.01)
filterE = subset(filterE, filterE$PValue<0.01)
filterE = subset(filterE, abs(filterE$logFC)>2)

filterF = subset(zF, zF$FDR<0.01)
filterF = subset(filterF, filterF$PValue<0.01)
filterF = subset(filterF, abs(filterF$logFC)>2)

filterG = subset(zG, zG$FDR<0.01)
filterG = subset(filterG, filterG$PValue<0.01)
filterG = subset(filterG, abs(filterG$logFC)>2)

#can check # DEGs quickly
dim(filterA)
dim(filterB)
dim(filterC)
dim(filterD)
dim(filterE)
dim(filterF)
dim(filterG)

#make new file of filtered reads.
write.table(filterA, "CASP+CIF_vs_CASP_filtered.txt")
write.table(filterB, "CASP_vs_pWS_filtered.txt")
write.table(filterC, "CASP+CIF_vs_pWS+CIF_filtered.txt")
write.table(filterD, "pWS vs pWS+CIF_filtered.txt")
write.table(filterE, "E30 vs 315.1_filtered.txt")
write.table(filterF, "CASP+CIF_vs_pWS_filtered.txt")
write.table(filterG, "CASP_vs_pWS+CIF_filtered.txt")

#how many genes in each direction?
summary(decideTestsDGE(DEtestA))
summary(decideTestsDGE(DEtestB))
summary(decideTestsDGE(DEtestC))
summary(decideTestsDGE(DEtestD))
summary(decideTestsDGE(DEtestE))
summary(decideTestsDGE(DEtestF))
summary(decideTestsDGE(DEtestG))




