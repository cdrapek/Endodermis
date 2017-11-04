setwd("~/Desktop/CodingDesktop")

datac <- read.csv("CASP+CIF_vs_pWS+CIF_filtered.csv", header=TRUE)
datan <- read.csv("CASP_vs_pWS_filtered.csv", header = TRUE)
dataref <- read.csv("E30 vs 315.1_filtered.csv", header=TRUE)

#which genes intersect?
#First, which genes are up vs down
cif_up = subset(datac, datac$logFC>0)
cif_dwn = subset(datac, datac$logFC<0)

non_up = subset(datan, datan$logFC>0)
non_dwn = subset(datan, datan$logFC<0)

#Then we look at gene ID in each

ID_cifup = cif_up[,1]
ID_cifdwn = cif_dwn[,1]
ID_nonup = non_up[,1]
ID_nondwn = non_dwn[,1]

#We then ask wich AGIds intersect in each list and get a total #
int_up = intersect(ID_cifup, ID_nonup)
sum_int_up = length(int_up)

int_dwn = intersect(ID_cifdwn, ID_nondwn)
sum_int_dwn = length(int_dwn)

total_int = sum_int_up + sum_int_dwn

write.table(int_up, "Overlap_DE_up-10-18-17.txt")
write.table(int_dwn, "Overlap_DE_dwn-10-18-17.txt")

#which gene IDs do not intersect?
pWSnone_unique_up = setdiff(ID_nonup, ID_cifup)
pWSnone_unique_dwn = setdiff(ID_nondwn, ID_cifdwn)
pWSCIF_unique_up = setdiff(ID_cifup, ID_nonup)
pWSCIF_unique_dwn = setdiff(ID_cifdwn, ID_nondwn)

write.table(pWSnone_unique_up, "ATGs_DE-UP_noCIFonly.txt")
write.table(pWSnone_unique_dwn, "ATGs_DE-dwn-noCIFonly.txt")
write.table(pWSCIF_unique_up, "ATGs_DE_UP_CIFonly.txt")
write.table(pWSCIF_unique_dwn, "ATGs_DE_dwn-CIFonly.txt")

library ('VennDiagram')
library(extrafont)

A = length(datac$agID)
B = length(datan$agID)

#Basic Venn Plot of overlapping DEGs as result of SHR in the epidermis
venn.plot <- draw.pairwise.venn(A, B, total_int, fill = c("darkorchid", "darkolivegreen3"), alpha = rep(0.5, 2), scaled = FALSE, fontfamily ="Arial", cat.fontfamily = "Arial", rotation.degree = 0, cex=1.5, fontface = "bold", cat.cex = 1.2, cat.pos = 180, cat.dist = 0.03, margin = 0.2)
grid.draw(venn.plot);
grid.newpage();

#Now we will ask a different question
#How many of these DEGs are also DEGs when we compare SCR and 315?
#do intersects above intersect with E30 vs 315?
ref_up = subset(dataref, dataref$logFC>0)
ref_dwn = subset(dataref, dataref$logFC<0)

ID_refup = ref_up[,1]
ID_refdwn = ref_dwn[,1]

ref_int_up = intersect(ID_refup, int_up)
sum_ref_int_up = length(ref_int_up)

ref_int_dwn = intersect(ID_refdwn, int_dwn)
ref_sum_int_dwn = length(ref_int_dwn)

DEGsOverlapRef = c(ref_int_dwn, ref_int_up)


write.table(ref_int_up, "DEG_ref_vs_DEG+-CIF_up.txt")
write.table(ref_int_dwn, "DEG_ref_vs_DEG+-CIF_dwn.txt")

#how about intersect with CIF or no CIF alone with SCR vs 315?
ref_cif_up = intersect(ID_refup, ID_cifup)
ref_cif_dwn = intersect(ID_refdwn, ID_cifdwn)

ref_non_up = intersect(ID_refup, ID_nonup)
ref_non_dwn = intersect(ID_refdwn, ID_nondwn)

CIF_overlapRef = c(ref_cif_up, ref_cif_dwn)
NON_overlapRef = c(ref_non_up, ref_non_dwn)

write.table(ref_cif_up, "DEGref_v_DEGcif_up.txt")
write.table(ref_cif_dwn, "DEGref_v_DEGcif_dwn.txt")
write.table(ref_non_up, "DEGref_v_DEGnone_up.txt")
write.table(ref_non_dwn, "DEGref_v_DEGnone_dwn.txt")

#Make 3 way Venn Diagram
A = length(datac$CASPCIFvspWERCIF)
B = length(datan$pWSnone.vs.CASPnone)
C = length(dataref$atgId)
draw.triple.venn(A, B, C, total_int, length(NON_overlapRef), length(CIF_overlapRef), length(DEGsOverlapRef), fill = c("darkorchid", "olivedrab", "cornflowerblue"), col = "transparent", fontfamily ="Arial", cat.fontfamily = "Arial", rotation.degree = 0, cex=1.5, fontface = "bold", cat.cex = 1.2, cat.pos = 180, cat.dist = 0.03, margin = 0.2)
grid.newpage();
#any go in opposite direction?
int_opp1 = intersect(ID_refup, int_dwn)
sum_opp1 = length(int_opp1)

int_opp2 = intersect(ID_refdwn, int_up)
sum_opp2 = length(int_opp2)

write.table(int_opp2, "DEGrefup_DEG+-cifdwn.txt")
write.table(int_opp2, "DEGrefdwn_DEG+-cif_up.txt")

#which gene IDs do not intersect?
DEGncif_unique_up = setdiff(int_up, ID_refup)
DEGncif_unique_dwn = setdiff(int_dwn, ID_refdwn)
DEGref_unique_up = setdiff(ID_refup, int_up)
DEGref_unique_dwn = setdiff(ID_refdwn, int_dwn)

write.table(DEGncif_unique_up, "DEG_+-CIF_missinginREF_up.txt")
write.table(DEGncif_unique_dwn, "DEG_+-CIF_missinginREF_dwn.txt")
write.table(DEGref_unique_up, "DEGinREF_missingin+-CIF_up.txt")
write.table(DEGref_unique_dwn, "DEGinREF_missingin+-CIF_dwn.txt")






