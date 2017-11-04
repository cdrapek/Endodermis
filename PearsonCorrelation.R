setwd("~/Desktop")

data <- read.csv("Normalized.csv", header=TRUE, row.names=1)
head(data)

pearcor <- cor(data)
round(pearcor, 2)

library("corrplot")
corrplot(pearcor, type = "lower", method="color", title = "Pearson Correlation of Libraries", outline = TRUE, cl.lim=c(0,1), tl.col="black", col=colorRampPalette(c("blue","black","white"))(500))
datatable = capture.output(pearcor)
datatable
