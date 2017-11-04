install.packages("plotly")
install.packages("RColorBrewer")
install.packages("gplots")
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(gplots)

#if plot gives an error about invalid graphics, reset by running this command
dev.off()
setwd("~/Desktop/CodingDesktop")

#Read in FPKMs. This is how we will cluster.
FPKM <- read.csv("RM3_FPKM.csv", header=TRUE, row.names=1)

#Find our DEGs between CASP/pWS+CIF and CASP/pWSnoCIF
datac <- read.csv("CASP+CIF_vs_pWS+CIF_filtered.csv", header=TRUE)
datan <- read.csv("CASP_vs_pWS_filtered.csv", header = TRUE)

#First give direction of FC for comparing
cif_up = subset(datac, datac$logFC>0)
cif_dwn = subset(datac, datac$logFC<0)

non_up = subset(datan, datan$logFC>0)
non_dwn = subset(datan, datan$logFC<0)

#then call AgID
ID_cifup = cif_up[,1]
ID_cifdwn = cif_dwn[,1]
ID_nonup = non_up[,1]
ID_nondwn = non_dwn[,1]

#compare AgID lists
int_up = intersect(ID_cifup, ID_nonup)
int_dwn = intersect(ID_cifdwn, ID_nondwn)
#Make into one AgID list for heatmap 
intersect = c(int_up, int_dwn)

#How are these AgIDs expressed in rootmap3.0? 
RootMap = data.frame(FPKM)
DEGsMap = RootMap[row.names(RootMap) %in% intersect,]
write.table (DEGsMap, "OverlapDEGs_wtRootMap-10-18-17.txt")

#Making the heat map
#color pallete from http://www.r-graph-gallery.com/215-the-heatmap-function/
#colors and labels for all the heat maps

labcol = c("Whole Root", "Epidermis/LRC", "Hair Cells", "Non-hair cells", "Columella", "Cortex(CO2)", "Cortex(315)", "Endodermis(E30)", "Endodermis(SCR)", "Phloem Pole Pericycle","Xylem","Protophloem","Protoxylem","Proto and metaphloem","Stele", "QC")
coul = colorRampPalette(brewer.pal(8, "PiYG"))(25)

DEGsMatrix = as.matrix(DEGsMap)
heatmap(DEGsMatrix, scale="row", col=coul, labCol = labcol)

#Export heatmap and clear plot
grid.newpage();

#How are the DEGs overlapping with SCRvs315 or missing in SCRvs315 expressed in RootMap?

#First, call intersecting/different genes
dataref <- read.csv("filtered_results-SCRvs315.csv", header=TRUE)
ref_up = subset(dataref, dataref$logFC>0)
ref_dwn = subset(dataref, dataref$logFC<0)

ID_refup = ref_up[,1]
ID_refdwn = ref_dwn[,1]

ref_int_up = intersect(ID_refup, int_up)
ref_int_dwn = intersect(ID_refdwn, int_dwn)

ref_intersect = c(ref_int_up, ref_int_dwn)
aberrant = setdiff(intersect, ref_intersect)

#How are the intersecting SCR/315 vs pWERSHRCASP+-CIF AgIDs expressed in rootmap3.0?
RefMap = RootMap[row.names(RootMap) %in% ref_intersect,]
RefMatrix = as.matrix(RefMap)
heatmap(RefMatrix, scale="row", col=coul, labCol = labcol)

#How are the "abberant" (missing DEGs in SCR/315) AdIDs expressed in rootmap3.0?
AbMap = RootMap[row.names(RootMap) %in% aberrant,]
AbMatrix = as.matrix(AbMap)
heatmap(AbMatrix, scale="row", col=coul, labCol = labcol)
  