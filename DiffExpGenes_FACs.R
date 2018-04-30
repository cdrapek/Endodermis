#RNA-seq data from CASP1-mCherry FACs from the following genotypes/treatments:
#CASP1-mCherry (wild-type) 
#CASP1-mCherry +CIF2 
#pWERSHRGFP/CASP1-mCherry
#pWERSHRGFP/CASP1-mcherry +CIF2 samples
#
#Samples were sequenced on platform HiSeq4000, 50bp SR.


#Data from RootMap 3 can be accessed from publication Li, Yamada et al 2016.
#Here, I am using it to analyze/filter the same way as the FACs data.
#The normalized/filtered data is available in Table S1.
#If you would like raw count data, please contact correspondng author. 

#Running the edgeR R package on 
biocLite('edgeR') #To install. Do not need to repeat once EdgeR is installed. 
library(edgeR)

#Set working directory and read in files
#this file countains raw counts from sequencing described above and
#raw counts from RootMap 3 (RM3) from E30 (mature endodermis) and 315(mature cortex) sorted libraries

setwd("~/Desktop/")
x <- read.csv("RawCounts_FACs.csv", header = TRUE, row.names=1)

#Make groups of samples (replicates)
groups=as.factor(c(rep("CM CIF+",3),rep("CM CIF-",3),rep("pWS CIF+",3),rep("pWS CIF-",3)))

#save gene names
genes=rownames(x)

#make edgeR object list
x <- DGEList(counts = x, group=groups, genes = genes)

#remove genes that are low/not expressed
#to do this, make cpm and lcpm
cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm))

#on the lcpm plot, look at the left tail. The peak below -5 x intersect are the non-expressed genes.

#keep only genes that are expressed
dim(x)
keep.exprs <- rowSums(cpm>1)>=7
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x) #compare to dim(x) above

#check new distribution
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm))

#normalize data after removing low expressed genes
#this adds normalization factor in the edgeR object
x <- calcNormFactors(x)
write.csv(cpm(x), "Norm_CPM.csv")

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
#tests of interests are: 
#effect of CIF2 on CASP1 cells from wild-type plants (CM.NvsCM.CIF)
#effect of CIF2 on cASP1 cells from pWERSHR plants (pWS.NvspWS.CIF)
#effect of ectopic SHR on CASP1 cells (CM.NvspWS.N) 
#effect of ectopic SHR on CASP1 cells that have been treated with CIF2 (CM.CIFvspWS.CIF)
#synergistic - removes effect CIF2 has on wild-type plants to determine if genes are regulated differently
#(pWS.CIF-pWS.N)-(CM.CIF-CM.N)

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

#Make table of all diff exp genes that also adds calculated FDR value
wtCIF = data.frame(topTags(qlf_wtCIF,n=Inf))
pWSCIF = data.frame(topTags(qlf_pWS_CIF,n=Inf))
pWS_wt = data.frame(topTags(qlf_pWS_wt,n=Inf))
pWSCIF_wtCIF = data.frame(topTags(qlf_pWS_wt_CIF,n=Inf))
eff = data.frame(topTags(qlf_eff,n=Inf))

#For many of the tests, the FDRs for all genes are insigifnificant.
#Therefore, I filtered based on pval and logFC.
#The two exceptions where there are significant FDRs are noted here.

fil_wtCIF = subset(wtCIF, wtCIF$PValue<0.01) 
fil_wtCIF = subset(fil_wtCIF, abs(fil_wtCIF$logFC)>0.3)

fil_pWSCIF = subset(pWSCIF, pWSCIF$PValue<0.01)
fil_pWSCIF = subset(fil_pWSCIF, abs(fil_pWSCIF$logFC)>0.3)

#FDR
fil_pWS_wt = subset(pWS_wt, pWS_wt$FDR<0.05) 
fil_pWS_wt = subset(fil_pWS_wt, abs(fil_pWS_wt$logFC)>0.3)

#FDR
fil_pWSCIF_wtCIF = subset(pWSCIF_wtCIF, pWSCIF_wtCIF$FDR<0.05) 
fil_pWSCIF_wtCIF = subset(fil_pWSCIF_wtCIF, abs(fil_pWSCIF_wtCIF$logFC)>0.3)

fil_eff = subset(eff, eff$PValue<0.01) 
fil_eff = subset(fil_eff, abs(fil_eff$logFC)>0.3)

#Write filtered results

write.csv(fil_wtCIF, "CM.NvsCM.CIF.csv")
write.csv(fil_pWSCIF, "pWS.NvspWS.CIF.csv")
write.csv(fil_pWS_wt, "WT.NvspWS.N.csv")
write.csv(fil_pWSCIF_wtCIF, "CM.CIFvspWS.CIF.csv")
write.csv(fil_eff, "CMeffpWS.csv")


####
#### DEGs from tests above expressed in RootMap 3
library(colorspace)

#Call in RootMap 3

FPKM <- read.csv("RM3_FPKM.csv", header=TRUE, row.names=1)
RM3 <- data.frame(FPKM)

#pWSCIF vs WTCIF Figure S5B
pWSCIF_up = subset(fil_pWSCIF_wtCIF, fil_pWSCIF_wtCIF$logFC>0.0)
pWSCIF_dn = subset(fil_pWSCIF_wtCIF, fil_pWSCIF_wtCIF$logFC<0.0)

genes_up = rownames(pWSCIF_up)
genes_dn = rownames(pWSCIF_dn)

exps_up = RM3[row.names(RM3) %in% genes_up,]
exps_dn = RM3[row.names(RM3) %in% genes_dn,]

#pick colors for heatmap
pal <- choose_palette() #I picked multiple hue divering, 30 colors

#heatmap up
coul = pal(30)
dat <- as.matrix(exps_up)
heatmap(dat, col = coul)

#heatmap down
coul = pal(30)
dat <- as.matrix(exps_dn)
heatmap(dat, col = coul)
