setwd("~/Desktop/")
library(RColorBrewer)
library(ggplot2)
library(reshape2)

data <- read.csv("ICI_mat.csv", row.names=1, header=TRUE)
celltype <- c("Hair Cell",	"Non-hair cell",	"Meristem epidermis",	"QC",	"Columella"	," Mature Cortex"	, "Meristem Cortex",	"Mature Endodermis",	"Meristem Endodermis",	"Stele",	"Pericycle",	"Meristem Xylem",	"Xylem",	"Phloem (S32)",	"Phloem (APL)")
samples <- rownames(data)
rownames(data) <- samples #for putting things in correct order
colnames(data) <- celltype
datalong <- c(data[,1],data[,2],data[,3],data[,4],data[,5],data[,6],data[,7],data[,8],data[,9],data[,10],data[,11],data[,12],data[,13],data[,14],data[,15])
df <- expand.grid(samples = samples, celltype = celltype)
df$ICI.val <- (datalong)
grid <- ggplot(data = df, aes(x = celltype, y = samples)) + geom_tile(aes(fill = ICI.val), color = "white")
grid

#change colors
grid2 <- grid + scale_fill_gradient(low = "blue4", high = "yellow") + theme_bw()

#add labels
grid3 <- grid2 + labs(fill = "ICI value") + ylab("Sample (Mean)") + xlab("Cell Types")

#add theme
grid4 <- grid3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
grid5 <- grid4 + theme(axis.text.x=element_text(angle=90,hjust=1))
