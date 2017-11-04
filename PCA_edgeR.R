#Running the edgeR R package on RNA-seq data from CASP1-mCherry sorted cells from CASP1-mCherry, CASP1-mCherry+CIF2, pWERSHRGFP/CASP1-mCherry and pWERSHRGFP/CASP1-mcherry+CIF2 samples. Samples were sequenced on the HiSeq4000, 50bp SR.
#Used to install EdgeR:
source("http://bioconductor.org/biocLite.R")
biocLite('edgeR')
biocLite("sva")
library(edgeR)
library(sva)
library(colorspace)

#Set your working directory
setwd("~/Desktop/CodingDesktop/")

#Read in input counts
##input file was raw counts from CASP/pWERSHR samples and raw counts from E30 and 315 cell-types from RootMap3
count <- read.csv("EdgeR_raw_CIF_MY_e30.csv", row.names=1)

#check the input file
head(count)
summary(count)

#RootMap3.0 is missing 3 ATGIDs present in CASP/CIF/pWERSHRdataset 
#We must omit this data
count <- na.omit(count)

#Name groups of biological reps
groups = c(rep("pWS+CIF", each=3), rep("CASP", each=3), rep("COR", each=3), rep("E30", each=3))

#generate data list object
y <- DGEList(counts = count, group = groups)
#remove genes that are no expressed
A <- aveLogCPM(y)
y1 <- y[A>1,]
#normalize data
y1 <- calcNormFactors(y1)
y2 <- estimateCommonDisp(y1)
y3 <- estimateTagwiseDisp(y2)
y3$common.dispersion

#make multidimensional scaling plot
logCPM <- cpm(y3, log=TRUE, prior.count=7)
pal <-choose_palette()
colors = rep(pal(4), each=3)
plotMDS(logCPM, top = 3000, gene.selection = "common", pch=19, col= alpha(colors, 0.5), cex=3, xlim=c(-2,3), ylim=c(-1.5,2))
legend("topleft", pch=19, col=(alpha(pal(4), 0.5)),
       legend=unique(groups), bty="n")

#we also tested if there was significant variation from the sequencing of CASP/pWERSHR at a different time than RootMap3.0
#we concluded normalizing data together was sufficient to overcome major variability in gene expression from sequencing at different times

#we tested for "batch effect" using a modification of Jeff Leek's method
#http://jtleek.com/genstats/inst/doc/02_13_batch-effects.html#adjusting-for-batch-effects-with-a-linear-model

#estimate batch variable with svaseq command 
pWS_CIF = c(rep(1,each=3), rep(0,each=9))
CASP = c(rep(0,each=3), rep(1,each=3), rep(0,each=6))
COR = c(rep(0,each=6), rep(1,each=3), rep(0,each=3))
E30 = c(rep(0,each=9), rep(1,each=3))
mod1 = mod1 = cbind(pWS_CIF, CASP, COR, E30)
mod0 = rep(1, each=12)
dat0 = as.matrix(y3)
svseq = svaseq(dat0, mod1, mod0)
#sva will estimate surrogate variables (this may take a few seconds)
#sva estimates 2 surrogate variables
#do the surrogate variables correlate with batch? 
batch = c(rep(1, each=6), rep(2, each=6))
summary(lm(svseq$sv ~ batch))

#the pvalues are both far above 0.05, meaning they do not significantly contribute to the variability in gene expression



