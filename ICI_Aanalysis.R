#Calculating Spec Scores from New Data
cisource <- as.matrix(read.delim("RM2_Means_MatchRM3.csv", row.names=1))
head(cisource)


#header <- c("APL", "APL", "APL","CO2","CO2", "CO2","COBL9", "COBL9","COBL9","COR", "COR","COR","E30", "E30","E30","GL2","GL2","GL2", "PET111", "PET111","PET111","S17", "S17","S17","S18", "S18","S18","S32", "S32","S32","S32", "S4","S4","S4","SCR","SCR","SCR", "WER",  "WER", "WER","WOL", "WOL","WOL","WOX5","WOX5","WOX5")
header <- c("AGL42","APL","COBL9","CORE","CORTEX","CYCD6","E30","GL2","J0571","J2661","J0121","LRC","PET111","RM1000","S4","S17","S18","S32","SCR","SUC2","WER","WOL","XYLEM_2501")

#Initial iteration could not assign 
spec_ci_RM2_MatchRM3 <- getAllSpec(cisource, header, medianfilter=0, cuts=FALSE, distshape=0)

#Calculating ICI from Unknowns
markers <- getMarkerList(spec_ci_RM2_MatchRM3,30,rownames(cisource))
data <- as.matrix(read.delim("FPKM_RM3_MatchRM2.csv", row.names=1))
head(data)
id_matrix <- getIdentity(data, spec_ci_RM2_V2, markers, TRUE, rownames(cisource))

write.csv(id_matrix[[1]], file = "Normalized_ID_Matrix.csv")
write.csv(id_matrix[[2]], file = "Raw_ICI_Score.csv")
write.csv(id_matrix[[3]], file = "Number_Markers_For_ICI_Score.csv")
write.csv(id_matrix[[4]], file = "p-value_ICI.csv")
write.csv(id_matrix[[5]], file = "Corrected_FDR.csv")
