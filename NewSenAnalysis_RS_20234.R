####libraries####
library(tidyr)
library(dplyr)
library(VennDiagram)
library(gplots)
library(remotes)
library(progeny)
library(tidyr)
library(dplyr)
library(gplots)
library(corrr)
library(Hmisc)
library(corrplot)
library(ggpubr)
library(factoextra)
library(annotables)
library(clusterProfiler)
library(enrichplot)
library(RColorBrewer)
library(KEGGREST)
library(RankProd)
library(ggplot2)
library(rstatix)
library(matrixStats)

setwd("/Volumes/Bex5TB/Sen23/NewAnalysis/Analysis2024")
AnalysisDir = ("/Volumes/Bex5TB/Sen23/NewAnalysis/Analysis2024/")

####open df####
Total_Data = read.csv("Total_filtered_dat_24Jan24.csv", row.names = 1)
Total_Data = mutate(Total_Data, Comparison = paste0(Study, "_", X))
X = Total_Data[which(!(duplicated(Total_Data$Comparison))),]

####create SENvCON####
SENvCON = filter(Total_Data, SENvCON == "Yes")
SENvCON = filter(SENvCON, Treatment == "none")
SENvCON = filter(SENvCON, Disease == "none")
SENvCON = filter(SENvCON, Control_type == "Prolif")

write.csv(SENvCON, "SENvCON.csv")
Y = SENvCON[which(!(duplicated(SENvCON$Comparison))),]

####Run IQ####
# Calculate IQR for each gene
iqr_df = SENvCON %>%
  group_by(Gene) %>%
  summarise(IQR_LogFC = IQR(LogFC))

# Merge the IQR information back to the original dataframe
SENvCON_IQR <- SENvCON %>%
  left_join(iqr_df, by = "Gene")

# Filter values outside of the IQR
SENvCON_IQR_filtered <- SENvCON_IQR %>%
  filter(LogFC >= quantile(LogFC, 0.25) - 1.5 * IQR_LogFC,
         LogFC <= quantile(LogFC, 0.75) + 1.5 * IQR_LogFC)

write.csv(SENvCON_IQR, "SENvCON_IQR.csv")
write.csv(SENvCON_IQR_filtered, "SENvCON_IQR_filtered.csv")

####zscore - NO LONGER USING####
SENvCON$z <- ave(SENvCON$LogFC, SENvCON$Sen_type, FUN=scale)
#creates dataframe with z-scores between 3.29 and -3.29
filteredSENvCON = filter(SENvCON, (z < 3.29 & z > -3.29))

setwd("/Volumes/Bex5TB/Sen23/NewAnalysis/NewAnalysis/zscore/")
write.csv(filteredSENvCON, "Zscore_SENvCON.csv")

####Venn and Heatmap####
####Venn####
#setwd("/Volumes/Bex5TB/Sen23/NewAnalysis/jan_16th")
#SENvCON_IQR_filtered = read.csv("/Volumes/Bex5TB/Sen23/NewAnalysis/jan_16th/SENvCON_IQR_filtered.csv", row.names = 1)
SENvCON_venn = filter(SENvCON_IQR_filtered, Sen_type == "DDIS" | Sen_type == "OIS" | Sen_type == "REP" |Sen_type == "Bystander")
write.csv(SENvCON_venn, "SENvCON_VENN.csv")

Get_SEN_SPECIFIC = function(Datframe, Name){
  SEN_data = filter(Datframe, Sen_type == Name)
  return(SEN_data)
}

SEN_Names = c("OIS", "DDIS", "REP", "BYS")

for(i in SEN_Names){
  assign(paste0(i, "_data"), Get_SEN_SPECIFIC(SENvCON_venn, i))
}

OIS_data = filter(SENvCON_venn, Sen_type == "OIS")
DDIS_data = filter(SENvCON_venn, Sen_type =="DDIS")
REP_data = filter(SENvCON_venn, Sen_type =="REP")
BYS_data = filter(SENvCON_venn, Sen_type =="Bystander")

#get medians for pvalues, logfc and megap
Get_medians = function(Dataframe, Name){
  Gene_list = unique(Dataframe$Gene)
  Init_dat = filter(Dataframe, Gene == Gene_list[1])
  MedianPval = median(Init_dat$Pvalues)
  MedianMegaP = median(Init_dat$MegaP)
  MedianLog = median(Init_dat$LogFC)
  Init_dat2 = as.data.frame(cbind(Gene = as.character(Gene_list[1]), MedianPval, MedianMegaP, MedianLog))
  Gene_list = Gene_list[2:length(Gene_list)]
  for(i in Gene_list){
    Added_dat = filter(Dataframe, Gene == i)
    MedianPval = median(Added_dat$Pvalues)
    MedianMegaP = median(Added_dat$MegaP)
    MedianLog = median(Added_dat$LogFC)
    Added = as.data.frame(cbind(Gene = i, MedianPval, MedianMegaP, MedianLog))
    Init_dat2 = dplyr::bind_rows(Init_dat2, Added)
  }
  #write.csv(Init_dat2, paste0(Name, "_PostSASP_medians.csv"))
  return(Init_dat2)
}

OIS_GeneSig = Get_medians(OIS_data, "OIS")
REP_GeneSig = Get_medians(REP_data, "REP")
DDIS_GeneSig = Get_medians(DDIS_data, "DDIS")
BYS_GeneSig = Get_medians(BYS_data, "Bystander")

Gene_list = c(as.character(OIS_GeneSig$Gene), as.character(DDIS_GeneSig$Gene), 
              as.character(REP_GeneSig$Gene), as.character(BYS_GeneSig$Gene))
Gene_list = unique(Gene_list)

#this adds all genes which aren't present as a sig gene to the dfs, with a pvalue =1 , mega p =0 and log =0. so not sig
Add_extra_genes = function(Datframe){
  Added_genes = Gene_list[which(!(Gene_list %in% Datframe$Gene))]
  Gene_frame = cbind(Gene = Added_genes[1], MedianPval = 1, MedianMegaP = 0, MedianLog = 0)
  Added_genes = Added_genes[2:length(Added_genes)]
  for(i in Added_genes){
    Added_frame = cbind(Gene = i, MedianPval = 1, MedianMegaP = 0, MedianLog = 0)
    Gene_frame = rbind(Gene_frame, Added_frame)
  }
  TotDat = rbind(Datframe, Gene_frame)
  TotDat$Gene = as.character(TotDat$Gene)
  TotDat = arrange(TotDat, Gene)
  TotDat = TotDat[which(TotDat$Gene != "NA"),]
  return(TotDat)
}

OIS_GeneSig = Add_extra_genes(OIS_GeneSig)
DDIS_GeneSig = Add_extra_genes(DDIS_GeneSig)
REP_GeneSig = Add_extra_genes(REP_GeneSig)
BYS_GeneSig = Add_extra_genes(BYS_GeneSig)

##IF using the MegaP values from the initial dataset
Get_MegaP = function(Datframe, Name){
  # Pos_neg = Datframe %>% mutate(Pos_Neg = if_else(.$MedianLog < 0, -1, 1))
  # Pos_neg = mutate(Pos_neg, MedianMegaP = Pos_Neg/MedianPval)
  # Pos_neg = Pos_neg[,c(1:4)]
  write.csv(Datframe, paste0(Name, "_SENvCON_medians.csv"))
  return(Datframe)
}

OIS_GeneSig = Get_MegaP(OIS_GeneSig, "OIS")
DDIS_GeneSig = Get_MegaP(DDIS_GeneSig, "DDIS")
REP_GeneSig = Get_MegaP(REP_GeneSig, "REP")
BYS_GeneSig = Get_MegaP(BYS_GeneSig, "BYS")

OIS_GeneSig$MedianMegaP = as.numeric(as.character(OIS_GeneSig$MedianMegaP))
DDIS_GeneSig$MedianMegaP = as.numeric(as.character(DDIS_GeneSig$MedianMegaP))
REP_GeneSig$MedianMegaP = as.numeric(as.character(REP_GeneSig$MedianMegaP))
BYS_GeneSig$MedianMegaP = as.numeric(as.character(BYS_GeneSig$MedianMegaP))

SigREPGenes = REP_GeneSig %>%
  mutate(Significant_MegaP = case_when(
    MedianMegaP <= -20 ~ 'Yes',
    MedianMegaP >= 20 ~ 'Yes',
    TRUE ~ 'No'
  ))
SigREPGenes = SigREPGenes[c("Gene", "Significant_MegaP")]

SigDDISGenes = DDIS_GeneSig %>%
  mutate(Significant_MegaP = case_when(
    MedianMegaP <= -20 ~ 'Yes',
    MedianMegaP >= 20 ~ 'Yes',
    TRUE ~ 'No'
  ))
SigDDISGenes = SigDDISGenes[c("Gene", "Significant_MegaP")]

SigOISGenes = OIS_GeneSig %>%
  mutate(Significant_MegaP = case_when(
    MedianMegaP <= -20 ~ 'Yes',
    MedianMegaP >= 20 ~ 'Yes',
    TRUE ~ 'No'
  ))
SigOISGenes = SigOISGenes[c("Gene", "Significant_MegaP")]

SigBYSGenes = BYS_GeneSig %>%
  mutate(Significant_MegaP = case_when(
    MedianMegaP <= -20 ~ 'Yes',
    MedianMegaP >= 20 ~ 'Yes',
    TRUE ~ 'No'
  ))
SigBYSGenes = SigBYSGenes[c("Gene", "Significant_MegaP")]

AllSenTypeComparison = cbind(SigOISGenes, SigDDISGenes$Significant_MegaP, SigREPGenes$Significant_MegaP,
                             SigBYSGenes$Significant_MegaP)

colnames(AllSenTypeComparison) = c("Gene", "OIS", "DDIS", "REP", "BYS")

write.csv(AllSenTypeComparison, "SigGenes_YesNo_Matrix_MegaP.csv")

Gene_listSIG = unique(AllSenTypeComparison$Gene)

OIS = c()
DDIS = c()
REP = c()
Bys = c()
OIS_DDIS =c()
OIS_REP = c()
OIS_Bys = c()
DDIS_REP = c()
DDIS_Bys = c()
REP_Bys = c()
O_D_R = c()
O_R_B = c()
D_R_B =c()
O_D_B = c()
O_D_R_B = c()

for(i in Gene_list){
  w = filter(AllSenTypeComparison, Gene == i)
  if(w[2] == "Yes" & w[3] == "No" & w[4] == "No" & w[5] == "No"){
    OIS = append(OIS, i)
  } else if(w[2] == "No" & w[3] == "Yes" & w[4] == "No" & w[5] == "No"){
    DDIS = append(DDIS, i)
  } else if(w[2] == "No" & w[3] == "No" & w[4] == "Yes" & w[5] == "No"){
    REP = append(REP, i)
  } else if(w[2] == "No" & w[3] == "No" & w[4] == "No" & w[5] == "Yes"){
    Bys = append(Bys, i)
  } else if(w[2] == "Yes" & w[3] == "Yes" & w[4] == "No" & w[5] == "No"){
    OIS_DDIS = append(OIS_DDIS, i)
  } else if(w[2] == "Yes" & w[3] == "No" & w[4] == "Yes" & w[5] == "No"){
    OIS_REP = append(OIS_REP, i)
  } else if(w[2] == "Yes" & w[3] == "No" & w[4] == "No" & w[5] == "Yes"){
    OIS_Bys = append(OIS_Bys, i)
  } else if(w[2] == "No" & w[3] == "Yes" & w[4] == "Yes" & w[5] == "No"){
    DDIS_REP = append(DDIS_REP, i)
  } else if(w[2] == "No" & w[3] == "Yes" & w[4] == "No" & w[5] == "Yes"){
    DDIS_Bys = append(DDIS_Bys, i)
  } else if(w[2] == "No" & w[3] == "No" & w[4] == "Yes" & w[5] == "Yes"){
    REP_Bys = append(REP_Bys, i)
  } else if(w[2] == "Yes" & w[3] == "Yes" & w[4] == "Yes" & w[5] == "No"){
    O_D_R = append(O_D_R, i)
  } else if(w[2] == "Yes" & w[3] == "Yes" & w[4] == "No" & w[5] == "Yes"){
    O_D_B = append(O_D_B, i)
  } else if(w[2] == "Yes" & w[3] == "No" & w[4] == "Yes" & w[5] == "Yes"){
    O_R_B = append(O_R_B, i)
  } else if(w[2] == "No" & w[3] == "Yes" & w[4] == "Yes" & w[5] == "Yes"){
    D_R_B = append(D_R_B, i)
  } else if(w[2] == "Yes" & w[3] == "Yes" & w[4] == "Yes" & w[5] == "Yes"){
    O_D_R_B = append(O_D_R_B, i)
  }
}

area1 = length(OIS) + length(OIS_DDIS) + length(OIS_REP) + length(OIS_Bys) + 
  length(O_D_R) + length(O_D_B) + length(O_R_B) + length(O_D_R_B)
area2 = length(DDIS) + length(OIS_DDIS) + length(DDIS_REP) + length(DDIS_Bys) + 
  length(O_D_R) + length(O_D_B) + length(D_R_B) + length(O_D_R_B)
area3 =  length(REP) + length(OIS_REP) + length(DDIS_REP) + length(REP_Bys) + 
  length(O_D_R) + length(D_R_B) + length(O_R_B) + length(O_D_R_B)
area4 =  length(Bys) + length(OIS_Bys) + length(DDIS_Bys) + length(REP_Bys) + 
  length(O_D_B) + length(O_R_B) + length(D_R_B) + length(O_D_R_B)
n12 = length(OIS_DDIS) + length(O_D_R) + length(O_D_B) + length(O_D_R_B)
n13 = length(OIS_REP) + length(O_D_R) + length(O_R_B) + length(O_D_R_B)
n14 = length(OIS_Bys) + length(O_D_B) + length(O_R_B) + length(O_D_R_B)
n23 = length(DDIS_REP) + length(O_D_R) + length(D_R_B) + length(O_D_R_B)
n24 = length(DDIS_Bys) + length(O_D_B) + length(D_R_B) + length(O_D_R_B)
n34 = length(REP_Bys) + length(O_R_B) + length(D_R_B) + length(O_D_R_B)
n123 = length(O_D_R) + length(O_D_R_B)
n124 = length(O_D_B) + length(O_D_R_B)
n134 = length(O_R_B) + length(O_D_R_B)
n234 = length(D_R_B) + length(O_D_R_B)
n1234 = length(O_D_R_B)

grid.newpage()
v = draw.quad.venn(area1=area1,
                   area2=area2,
                   area3=area3,
                   area4=area4,
                   n12=n12,
                   n13=n13,
                   n14=n14,
                   n23=n23,
                   n24=n24,
                   n34=n34,
                   n123=n123,
                   n124=n124,
                   n134=n134,
                   n234=n234,
                   n1234=n1234, 
                   category = c("OIS","DDIS","REP","BYS"),
                   col= "Black", fill=c("Red","Blue","Yellow","Green"),
                   cex = 1.5,
                   cat.cex = 2,
                   margin = 0.05)





####Heatmap####
Get_Top_genes = function(Datframe, Name){
  Top_genes = Datframe[which(Datframe$Gene %in% O_D_R_B),]
  Top_genes = arrange(Top_genes, Gene)
  colnames(Top_genes) = c("Gene", "Pvalues", "MegaP", "LogFC")
  Top_genes$Pvalues = as.numeric(as.character(Top_genes$Pvalues))
  Top_genes$MegaP = as.numeric(as.character(Top_genes$MegaP))
  Top_genes$LogFC = as.numeric(as.character(Top_genes$LogFC))
  Top_genes = Top_genes %>% pivot_longer(`Pvalues`:`LogFC`, names_to = "Test", values_to = Name)
  return(Top_genes)
}

Top_OIS = Get_Top_genes(OIS_GeneSig, "OIS")
Top_DDIS = Get_Top_genes(DDIS_GeneSig, "DDIS")
Top_REP = Get_Top_genes(REP_GeneSig, "REP")
Top_BYS = Get_Top_genes(BYS_GeneSig, "BYS")
Top_Genes = cbind(Top_OIS, DDIS = Top_DDIS$DDIS, REP = Top_REP$REP, BYS = Top_BYS$BYS)
#Top_Genes_sen = cbind(Top_OIS, DDIS = Top_DDIS$DDIS, REP = Top_REP$REP)
#For heatmap
setwd("/Volumes/Bex5TB/Sen23/NewAnalysis/Analysis2024/Venn_Heatmap")

write.csv(Top_Genes, "SigGenes_AllSentypes_SENvCON.csv")
Top_Genes2 = Top_Genes
#Top_Genes_sen2 = Top_Genes_sen
#write.csv(Top_Genes_sen, "SigGenes_ODR_SENvCON.csv")

#Now get complete table
Top_Genes = mutate(Top_Genes, Group = "O_D_R_B")
#Top_Genes_sen = mutate(Top_Genes_sen, Group = "O_D_R")
Sig_list = lst(OIS, DDIS, REP, Bys, OIS_DDIS, OIS_REP, OIS_Bys, DDIS_REP, DDIS_Bys, REP_Bys, 
               O_R_B, O_D_R, D_R_B, O_D_B)

for(i in 1:length(Sig_list)){
  Gene_list1 = Sig_list[[i]]
  Namy = names(Sig_list[i])
  #Gene_list1 = as.data.frame(Gene_list1)
  Gene_list1 = as.character(Gene_list1)
  #Gene_list1 = unlist(Gene_list1)
  print(Sig_list[i])
  print(names(Sig_list[i]))
  #Gene_list1 = as.character(Gene_list1)
  Next_genes = OIS_GeneSig[which(OIS_GeneSig$Gene %in% Gene_list1),]
  Next_genes = arrange(Next_genes, Gene)
  colnames(Next_genes) = c("Gene", "Pvalues", "MegaP", "LogFC")
  Next_genes$Pvalues = as.numeric(as.character(Next_genes$Pvalues))
  Next_genes$MegaP = as.numeric(as.character(Next_genes$MegaP))
  Next_genes$LogFC = as.numeric(as.character(Next_genes$LogFC))
  Next_OIS = Next_genes %>% pivot_longer(`Pvalues`:`LogFC`, names_to = "Test", values_to = "OIS")
  Next_OIS$Gene = as.character(Next_OIS$Gene)
  Next_OIS = arrange(Next_OIS, Gene)
  
  Next_genes = DDIS_GeneSig[which(DDIS_GeneSig$Gene %in% Gene_list1),]
  Next_genes = arrange(Next_genes, Gene)
  colnames(Next_genes) = c("Gene", "Pvalues", "MegaP", "LogFC")
  Next_genes$Pvalues = as.numeric(as.character(Next_genes$Pvalues))
  Next_genes$MegaP = as.numeric(as.character(Next_genes$MegaP))
  Next_genes$LogFC = as.numeric(as.character(Next_genes$LogFC))
  Next_DDIS = Next_genes %>% pivot_longer(`Pvalues`:`LogFC`, names_to = "Test", values_to = "DDIS")
  Next_DDIS$Gene = as.character(Next_DDIS$Gene)
  Next_DDIS = arrange(Next_DDIS, Gene)
  
  Next_genes = REP_GeneSig[which(REP_GeneSig$Gene %in% Gene_list1),]
  Next_genes = arrange(Next_genes, Gene)
  colnames(Next_genes) = c("Gene", "Pvalues", "MegaP", "LogFC")
  Next_genes$Pvalues = as.numeric(as.character(Next_genes$Pvalues))
  Next_genes$MegaP = as.numeric(as.character(Next_genes$MegaP))
  Next_genes$LogFC = as.numeric(as.character(Next_genes$LogFC))
  Next_REP = Next_genes %>% pivot_longer(`Pvalues`:`LogFC`, names_to = "Test", values_to = "REP")
  Next_REP$Gene = as.character(Next_REP$Gene)
  Next_REP = arrange(Next_REP, Gene)
  
  Next_genes = BYS_GeneSig[which(BYS_GeneSig$Gene %in% Gene_list1),]
  Next_genes = arrange(Next_genes, Gene)
  colnames(Next_genes) = c("Gene", "Pvalues", "MegaP", "LogFC")
  Next_genes$Pvalues = as.numeric(as.character(Next_genes$Pvalues))
  Next_genes$MegaP = as.numeric(as.character(Next_genes$MegaP))
  Next_genes$LogFC = as.numeric(as.character(Next_genes$LogFC))
  Next_BYS = Next_genes %>% pivot_longer(`Pvalues`:`LogFC`, names_to = "Test", values_to = "BYS")
  Next_BYS$Gene = as.character(Next_BYS$Gene)
  Next_BYS = arrange(Next_BYS, Gene)
  
  #Next_Genes = cbind(Next_OIS, DDIS = Next_DDIS$DDIS, REP = Next_REP$REP, BYS = Next_BYS$BYS, Group = "O_D_B")
  Next_Genes = cbind(Next_OIS, DDIS = Next_DDIS$DDIS, REP = Next_REP$REP, BYS = Next_BYS$BYS, Group = Namy)
  Top_Genes = rbind(Top_Genes, Next_Genes)
}

write.csv(Top_Genes, "SigGenes_SENvCON_filteredDirec.csv")

#Top_Genes2 = read.csv("SENvCON/SigGenes_AllSentypes_SENvCON_filteredDirec.csv", row.names = 1)
Top_Genes3 = filter(Top_Genes, Test == "MegaP")
rownames(Top_Genes3) = Top_Genes3$Gene
Top_Genes3 = Top_Genes3[, c(3:6)]
Top_Genes3 = data.matrix(Top_Genes3, rownames.force = TRUE)
column_order = c("OIS", "REP", "BYS", "DDIS")
Top_Genes3 = Top_Genes3[, column_order]
#Top_Genes3 = as.numeric(Top_Genes3)

#Top_Genes_sen3 = filter(Top_Genes_sen, Test == "MegaP")
#rownames(Top_Genes_sen3) = Top_Genes_sen3$Gene
#Top_Genes_sen3 = Top_Genes_sen3[, c(3:5)]
#Top_Genes_sen3 = data.matrix(Top_Genes_sen3, rownames.force = TRUE)
#column_order = c("OIS", "REP", "DDIS")
#Top_Genes_sen3 = Top_Genes_sen3[, column_order]

colors <- colorRampPalette(c("blue", "white", "red")) (9)
breaks <- c(-100000, -10000, -1000, -100, -20, 20, 100, 1000, 10000, 100000)

pdftotalname2 = "SigGenes_AllSentypes_heatmap.pdf"
pdf(pdftotalname2)
par(oma=c(4,0,0,0))
heatmap.2(Top_Genes, dendrogram = "none", Rowv = TRUE, Colv = FALSE, key = FALSE, col = colors, na.color = "black", 
          density.info = "none", scale = "none", trace = "none", cexCol = 1, cexRow = 0.7, breaks = breaks, margins =c(5,10),
          hclustfun = function(x) hclust(x, method = 'centroid'), distfun = function(x) dist(x, method = 'euclidean'))
dev.off()
print("f")
pdftotalname2 = "SigGenes_AllSentypes_heatmap.png"
png(pdftotalname2)
par(oma=c(4,0,0,0))
heatmap.2(Top_Genes3, dendrogram = "none", Rowv = TRUE, Colv = FALSE, key = FALSE, col = colors, na.color = "black", 
          density.info = "none", scale = "none", trace = "none", cexCol = 1, cexRow = 0.7, breaks = breaks, 
          hclustfun = function(x) hclust(x, method = 'centroid'), distfun = function(x) dist(x, method = 'euclidean'))
dev.off()

pdftotalname2 = "SigGenes_AllSentypes_heatmap.svg"
svg(pdftotalname2)
par(oma=c(4,0,0,0))
heatmap.2(Top_Genes3, dendrogram = "none", Rowv = TRUE, Colv = FALSE, key = FALSE, col = colors, na.color = "black", 
          density.info = "none", scale = "none", trace = "none", cexCol = 1, cexRow = 0.7, breaks = breaks, 
          hclustfun = function(x) hclust(x, method = 'centroid'), distfun = function(x) dist(x, method = 'euclidean'))
dev.off()


####gsea####
#Linux
#term2gene = read.gmt("/media/c0068011/Bex5TB/SenAnalysisFiles/h.all.v2023.2.Hs.symbols.gmt")
#mac
term2gene = read.gmt("/Volumes/Bex5TB/SenAnalysisFiles/h.all.v2023.2.Hs.symbols.gmt")

term2gene$term = gsub("HALLMARK_", "", term2gene$term)
grch38 = grch38

#for each run of the GSEA, need to change the Group, file names and the rows selected on row 445

Top_Genes3 = filter(Top_Genes, Group == "REP")
Top_Genes3 = filter(Top_Genes3, Test == "LogFC")
Top_Genes3 = mutate(Top_Genes3, Median = rowMedians(as.matrix(Top_Genes3[,c(5)])))

genelist = pull(Top_Genes3, Median)
names(genelist) = pull(Top_Genes3, Gene)
genelist = sort(genelist, decreasing = TRUE)

gsea_all = GSEA(geneList = genelist,
                exponent = 1,
                #nPerm = 10000,
                minGSSize = 5,
                maxGSSize = 500,
                pvalueCutoff = 1,
                pAdjustMethod = "none",
                TERM2GENE = term2gene,
                by = "fgsea")
gsea_all_df = as.data.frame(gsea_all)

gsea_0.1 = GSEA(geneList = genelist,
                exponent = 1,
                #nPerm = 10000,
                minGSSize = 5,
                maxGSSize = 500,
                pvalueCutoff = 0.1,
                pAdjustMethod = "none",
                TERM2GENE = term2gene,
                by = "fgsea")
gsea_0.1_df = as.data.frame(gsea_0.1)

gsea_sig = GSEA(geneList = genelist,
                exponent = 1,
                #nPerm = 10000,
                minGSSize = 5,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "none",
                TERM2GENE = term2gene,
                by = "fgsea")
gsea_sig_df = as.data.frame(gsea_sig)

write.csv(gsea_all_df, paste0("Lung", "_all_gsea.csv"))
write.csv(gsea_sig_df, paste0("Lung", "_sig_gsea.csv"))
write.csv(gsea_0.1_df, paste0("Lung", "_0.1_gsea.csv"))

if(nrow(gsea_all_df) > 0){
  # Dotplot
  g1 = dotplot(gsea_all, color = "pvalue",  showCategory=nrow(gsea_all), split=".sign") + facet_grid(.~.sign)+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(text = element_text(size = 16),
          plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("Lung", "_all_dotplot.pdf"), g1. width = 10, height =12)
  # Ridge plot
  g3 = ridgeplot(gsea_all, fill = "pvalue") + labs(x = "enrichment distribution")+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("Lung", "_all_ridgeplot.pdf"), g3, width = 10, height =12)
}

if(nrow(gsea_sig_df) > 0){
  ## Dotplot
  g1 = dotplot(gsea_sig, color = "pvalue",  showCategory=nrow(gsea_sig), split=".sign") + facet_grid(.~.sign)+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(text = element_text(size = 16),
          plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("Lung", "_sig_dotplot.pdf"), g1, width = 10, height =12)
  ## Ridge plot
  g3 = ridgeplot(gsea_sig, fill = "pvalue") + labs(x = "enrichment distribution")+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("Lung","_sig_ridgeplot.pdf"), g3, width = 10, height =12)
}

if(nrow(gsea_0.1_df) > 0){
  ## Dotplot
  g1 = dotplot(gsea_0.1, color = "pvalue",  showCategory=nrow(gsea_0.1), split=".sign") + facet_grid(.~.sign)+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(text = element_text(size = 16),
          plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("Lung", "_sig0.1_dotplot.pdf"), g1, width = 10, height =12)
  ## Ridge plot
  g3 = ridgeplot(gsea_0.1, fill = "pvalue") + labs(x = "enrichment distribution")+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("Lung","_sig0.1_ridgeplot.pdf"), g3, width = 10, height =12)
}

####unknown currently#####
##Clustering analysis
SENvCON$Study = as.factor(SENvCON$Study)
SENvCON$Study = as.numeric(SENvCON$Study)
unique(SENvCON$Sen_type)
Pivoted_data_MegaP = mutate(SENvCON, Study_Comp = paste0(as.factor(Study), "_", Sen_type, Organ, Cell_line, numeric_time))
Pivoted_data_MegaP = mutate(Pivoted_data_MegaP, Study_Sentype = paste0(Study, "_", Sen_type), .before = 2)

Pivoted_data_MegaP = Pivoted_data_MegaP[, c(1,2,47,50)]

Pivoted_data_MegaP$Study_Comp = gsub("OIS", "O", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("REP", "R", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("DDIS", "D", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("Bystander", "B", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("NIS", "N", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("RNIS", "Nr", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("NBIS", "Nb", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("RiboMature", "Rm", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("dNTP", "Np", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("CR", "Cr", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("OSKM", "Ok", Pivoted_data_MegaP$Study_Comp)

Pivoted_data_MegaP$Study_Comp = gsub("Lung", "L", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("Skin", "S", Pivoted_data_MegaP$Study_Comp)

Pivoted_data_MegaP$Study_Comp = gsub("IMR", "I", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("BJ", "B", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("WI38", "W", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("LF1", "L", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("HCA2", "H", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("FL2", "F", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("Tig3", "T", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("HFF", "Hf", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("MRC", "M", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("CAF", "C", Pivoted_data_MegaP$Study_Comp)
Pivoted_data_MegaP$Study_Comp = gsub("HDF", "Hd", Pivoted_data_MegaP$Study_Comp)

#NOW get rid of studies that have the same comparisons outside those (of interest) in Study_Comp
Pivoted_data_MegaP = mutate(Pivoted_data_MegaP, Comp_Gene = paste0(Study_Comp, "_", Gene))
Pivoted_data_MegaP = Pivoted_data_MegaP[which(!(duplicated(Pivoted_data_MegaP$Comp_Gene))),]
Pivoted_data_MegaP2 = Pivoted_data_MegaP[,c(1,2,3,4)]

P_P = pivot_wider(Pivoted_data_MegaP2, names_from = Gene, values_from = MegaP)
write.csv(P_P, "Pivoted_genes_SENvCON_unfiltered_withSentype.csv")
#Comps_Pivoted = Pivoted_data_MegaP$Study_Comp
#rownames(P_P) = Pivoted_data_MegaP$Study_Comp

Pivoted_data_Matrix = P_P[, c(3:ncol(P_P))]
Pivoted_data_Matrix = data.matrix(Pivoted_data_Matrix)
head(Pivoted_data_Matrix[,c(1:10)])

m = apply(Pivoted_data_Matrix,2,mean)
s= apply(Pivoted_data_Matrix,2,sd)
Pivoted_data_Matrix = scale(Pivoted_data_Matrix,m,s)

distance = dist(Pivoted_data_Matrix)
#Cluster with complete linkage
hcc = hclust(distance)
png(file="Complete_linkage_cluster_SENvCON.png")
plot(hcc, labels=P_P$Study_Sentype, xlab = "Complete linkage clustering", cex = 0.3, hang = -1) #, hang=-1
dev.off()
svg(file="Complete_linkage_cluster_All_detail_SENvCON.svg")
plot(hcc, labels=P_P$Study_Comp, xlab = "Complete linkage clustering", cex = 0.275, hang = -1) #, hang=-1
dev.off()
#Cluster with average linkage
hca = hclust(distance, method = "average")
svg(file="Average_linkage_cluster_SENvCON.svg")
plot(hca, labels=P_P$Study_Sentype, xlab = "Average linkage clustering", cex = 0.275, hang =-1)
dev.off()
member_c = cutree(hcc,60)
member_a = cutree(hca,60)
R = as.data.frame.matrix(table(member_c, member_a))
#Cluster means
Cluster_means = aggregate(Pivoted_data_Matrix, list(member_c), mean)

library(cluster)
plot(silhouette(cutree(hcc,3),distance))

Pivoted_data_Matrix2 = Pivoted_data_Matrix[ , colSums(is.na(Pivoted_data_Matrix)) == 0]

wss = (nrow(Pivoted_data_Matrix2)-1)*sum(apply(Pivoted_data_Matrix,2,var))
for(i in 2:50){ wss[i] = sum(kmeans(Pivoted_data_Matrix2, centers=i)$withinss)}
svg(file="Scree_plot_clusters.svg")
plot(1:50, wss, type = "b", xlab = "Number of clusters", ylab = "Within group SS")
dev.off()

#K means clustering
kc = kmeans(Pivoted_data_Matrix2, 6)
Cluster_list = kc[["cluster"]]
P_P = cbind(Cluster = Cluster_list, P_P)
Cluster_table = P_P[, c(1:3)]
write.csv(Cluster_table, "Table_clusters_6_kmeans.csv")
plot(TP53,TGFB1, Pivoted_data_MegaP, col = kc$cluster)
####Temporal gene graphs####
setwd("/Volumes/Bex5TB/Sen23/NewAnalysis/Analysis2024")
SENvCON = read.csv("SENvCON_VENN.csv", row.names =1)

#create df with Sen_types of interest
sen3 = filter(SENvCON, Sen_type == "DDIS" | Sen_type == "OIS" | Sen_type == "REP")

#make NA = 0 when sen_type = REP
sen3 <- sen3 %>%
  mutate(numeric_time = if_else(Sen_type == 'REP' & is.na(numeric_time), 0, numeric_time))

#Creating numeric_time2 column so every line has a timegroup
sen3 = sen3 %>%
  mutate(numeric_time2 = case_when(
    Sen_type %in% c("DDIS", "OIS") ~ case_when(
      numeric_time <= 4 ~ "0-4",
      numeric_time > 4 & numeric_time <= 7 ~ "5-7",
      numeric_time > 7 & numeric_time <= 11 ~ "8-11",
      numeric_time > 11 & numeric_time <= 14 ~ "12-14",
      numeric_time > 14 ~ "15+",
      is.na(numeric_time) ~ "NA",
      TRUE ~ NA_character_
    ),
    Sen_type == "REP" ~ case_when(
      numeric_time < 41 ~ "0-40",
      numeric_time >= 41 ~ "41+",
      TRUE ~ NA_character_
    ),
    TRUE ~ NA_character_
  ))

unique(sen3$numeric_time2)

TimelineDO = filter(sen3, Sen_type == "DDIS" | Sen_type == "OIS")
TimelineDO = filter(TimelineDO, numeric_time != "NA")
TimelineDO = mutate(TimelineDO, SEN = "DDIS/OIS")
unique(TimelineDO$numeric_time2)
#create REP df
TimelineR = filter(sen3, Sen_type == "REP")
TimelineR = mutate(TimelineR, SEN = "REP")
unique(TimelineR$numeric_time2)
#bind both dfs together
sen3 = rbind(TimelineDO, TimelineR)
#removing 15+ times (done by JW in original analysis)
#RS 14Dec23 11:00 > with zscore hopefully 15+ shouldn't cause major outliers skewing outputs
#RS 14Dec23 11:40 > As only DDIS has 15+, remove from df
sen3 = filter(sen3, numeric_time2 != "15+")

# Create a vector with unique values of numeric_time2
Time_list = c("0-4", "5-7", "8-11", "12-14", "0-40", "41+")
# Order the levels of numeric_time2 based on Time_list
sen3$numeric_time2 = factor(sen3$numeric_time2, levels = Time_list)
unique(sen3$numeric_time2)

Create_pathway_graphs = function(pathway_data, expression_data){
  for(i in pathway_data){
    Gene_data = filter(expression_data, Gene == i)
    Gene_name = paste0(AnalysisDir, "/", i)
    print(i)
    Gene_graph = ggplot(Gene_data, aes(x= factor(numeric_time2), y=LogFC, fill=Sen_type))+
      stat_boxplot(geom ='errorbar') +
      geom_boxplot() + #outlier.shape = NA
      #geom_jitter(width=0.1,alpha=0.2) +
      stat_compare_means(method = "wilcox.test", label = "p.signif", size =4,
                         symnum.args = list(cutpoint = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                                            symbols = c("****", "***", "**", "*", "ns")),
                         position = position_nudge(y = 0.1)) +  # Add Wilcoxon test results
      labs(x = "Treatment", y = "Log FC", fill = "Senescence", title = i)+
      #ylim(Min_Numb, Max_Numb)+
      #facet_wrap(~REP) +
      facet_grid(cols = vars(SEN), scales = "free_x", space = "free") +
      scale_x_discrete("Days after induction", breaks=factor(Time_list), drop=TRUE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), panel.background = element_blank(), 
            strip.background = element_blank(), strip.text = element_blank(), 
            axis.text.y = element_text(size = 14), axis.title=element_text(size=14),
            legend.text=element_text(size=13), legend.title=element_text(size=14),
            axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))
    ggsave(paste0(Gene_name, "_Timeline.png"))
    ggsave(paste0(Gene_name, "_Timeline.pdf"))
    #ggsave(paste0(Gene_name, "_Timeline.svg"))#for then wilcox.test label = "p.signif"
    #ggsave(paste0(Gene_name, "_Timeline_pvalue.png")) #for when wilcox.test label = "p"
  }
}

Gene_list = c("TP53", "IL6","IL1B", "TGFB1", "TGFBR1", "COL1A1", "PDGFA", "ACTA2", "CXCL8", "HES1", "NOTCH1", "MAPK14",
              "CDKN2A", "ATM", "CHEK1", "MDM2", "GADD45A", "CDKN1A", "ATR", "CHEK2", "CDC25A", "HEY1")

AnalysisDir = "/Volumes/Bex5TB/Sen23/NewAnalysis/Analysis2024/Temporal"
Create_pathway_graphs(Gene_list, sen3)

####Temporal KDs####

KDs = filter(Total_Data, SENvCON == "Yes")
KDs = filter(KDs, Control_type == "Prolif")
KDs = filter(KDs, Sen_type == "OIS" | Sen_type == "DDIS")

iqr_df = KDs %>%
  group_by(Gene) %>%
  summarise(IQR_LogFC = IQR(LogFC))

# Merge the IQR information back to the original dataframe
IQR <- KDs %>%
  left_join(iqr_df, by = "Gene")

# Filter values outside of the IQR
IQR_filtered <- IQR %>%
  filter(LogFC >= quantile(LogFC, 0.25) - 1.5 * IQR_LogFC,
         LogFC <= quantile(LogFC, 0.75) + 1.5 * IQR_LogFC)

KDs = filter(IQR_filtered, Gene_down == "p53" | Gene_down == "RELA" | Gene_down == "none")
KDs = filter(KDs, numeric_time !=28)
KDs = filter(KDs, Control_gene_down == "none" | Control_gene_down == "p53" | Control_gene_down == "RELA")
KDs$Gene_down = gsub("none", "None", KDs$Gene_down)
KDs$Control_gene_down = gsub("none", "None", KDs$Control_gene_down)
KDs$Gene_down = gsub("p53", "TP53", KDs$Gene_down)
KDs = filter(KDs, Control_gene_down == "None")

setwd("/Volumes/Bex5TB/Sen23/NewAnalysis/Analysis2024/KD")
#write.csv(KDs, "Filtered_KD_df.csv")

Sen_list = c("DDIS", "OIS")

KDs$Sen_type = factor(KDs$Sen_type, levels = Sen_list)
unique(KDs$Sen_type)

Create_pathway_graphs_p53KDs = function(pathway_data, expression_data){
  for(i in pathway_data){
    Gene_data = filter(expression_data, Gene == i)
    Gene_name = paste0(AnalysisDir, "/", i)
    print(i)
    Gene_graph = ggplot(Gene_data, aes(x= factor(Sen_type), y=LogFC, fill=Gene_down))+
      stat_boxplot(geom ='errorbar') +
      geom_boxplot() + #outlier.shape = NA +
      stat_compare_means(method = "wilcox.test", label = "p.signif", size =4,
                         symnum.args = list(cutpoint = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                                            symbols = c("****", "***", "**", "*", "ns")),
                         position = position_nudge(y = 0.1)) +  # Add Wilcoxon test results
      labs(x = "Senescence", y = "Log FC", fill = "Gene Inhibition", title = i)+
      scale_x_discrete("Senescence", breaks=factor(Sen_list), drop=TRUE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), panel.background = element_blank(), 
            strip.background = element_blank(), strip.text = element_blank(), 
            axis.text.y = element_text(size = 14), axis.title=element_text(size=14),
            legend.text=element_text(size=13), legend.title=element_text(size=14),
            axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))
    ggsave(paste0(Gene_name, "_p53KD.png"))
    ggsave(paste0(Gene_name, "_p53KD.pdf"))
  }
}

Create_pathway_graphs_RELAKDs = function(pathway_data, expression_data){
  for(i in pathway_data){
    Gene_data = filter(expression_data, Gene == i)
    Gene_name = paste0(AnalysisDir, "/", i)
    print(i)
    Gene_graph = ggplot(Gene_data, aes(x= factor(Sen_type), y=LogFC, fill=Gene_down))+
      stat_boxplot(geom ='errorbar') +
      geom_boxplot() + #outlier.shape = NA +
      stat_compare_means(method = "wilcox.test", label = "p.signif", size =4,
                         symnum.args = list(cutpoint = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
                                            symbols = c("****", "***", "**", "*", "ns")),
                         position = position_nudge(y = 0.1)) +  # Add Wilcoxon test results
      labs(x = "Senescence", y = "Log FC", fill = "Gene Inhibition", title = i)+
      scale_x_discrete("Senescence", breaks=factor(Sen_list), drop=TRUE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), panel.background = element_blank(), 
            strip.background = element_blank(), strip.text = element_blank(), 
            axis.text.y = element_text(size = 14), axis.title=element_text(size=14),
            legend.text=element_text(size=13), legend.title=element_text(size=14),
            axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))
    ggsave(paste0(Gene_name, "_RELAKD.png"))
    ggsave(paste0(Gene_name, "_RELAKD.pdf"))
  }
}

Gene_list = c("ATM", "CHEK1", "TP53", "MDM2", "GADD45A", "GADD45B", "CDKN1A","CDKN2A", "MAPK14", "IL1B", "IL6", "CXCL8", "TGFB1")

p53 = filter(KDs, Gene_down == "TP53" | Gene_down == "None")
rela = filter(KDs, Gene_down == "RELA" | Gene_down == "None")

AnalysisDir = ("/Volumes/Bex5TB/Sen23/NewAnalysis/Analysis2024/KD")

Create_pathway_graphs_p53KDs(Gene_list, p53)
Create_pathway_graphs_RELAKDs(Gene_list, rela)

#write.csv(p53, "Filtered_p53KD.csv")
#write.csv(rela, "Filtered_relaKD.csv")

####Cell line comparisons####
Cell_type = filter(sen3, Sen_type != "REP")
Cell_type = filter(Cell_type, (numeric_time > 4 & numeric_time < 12))
unique(Cell_type$Cell_line)

Cell_list = c("BJ" , "HCA2", "HDF", "HFF", "IMR", "MRC5", "Tig3", "WI38")
Cell_type$Cell_line = factor(Cell_type$Cell_line, levels = Cell_list)
unique(Cell_type$Cell_line)
Sen_list = c("DDIS", "OIS")
Cell_type$Sen_type = factor(Cell_type$Sen_type, levels = Sen_list)

Create_pathway_graphs_cell = function(pathway_data, expression_data){
  for(i in pathway_data){
    Gene_data = filter(expression_data, Gene == i)
    Gene_name = paste0(AnalysisDir, "/", i)
    print(i)
    Gene_graph = ggplot(Gene_data, aes(x= factor(Sen_type), y=LogFC, fill=factor(Cell_line, levels =Cell_list)))+
      stat_boxplot(geom ='errorbar') +
      geom_boxplot() +
      labs(x = "", y = "Log FC", fill = "Cell line", title = i)+
      #ylim(Min_Numb, Max_Numb)+
      #facet_wrap(~REP) +
      #facet_grid(cols = vars(REP), scales = "free_x", space = "free") +
      scale_x_discrete("", breaks=factor(Sen_list), drop=TRUE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), panel.background = element_blank(), 
            strip.background = element_blank(), strip.text = element_blank(), 
            axis.text.y = element_text(size = 14), axis.title=element_text(size=14),
            legend.text=element_text(size=13), legend.title=element_text(size=14),
            axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))
    ggsave(paste0(Gene_name, "_IQR_CellLine.png"))
  }
}

Gene_list = c("CXCL8", "IL6","IL1B", "TGFB1", "TGFBR1", "COL1A1", "PDGFA", "ACTA2", "RB1", "CDKN1A", "CDKN2A")
AnalysisDir =  ("/Volumes/Bex5TB/Sen23/NewAnalysis/Analysis2024/Cell")
Create_pathway_graphs_cell(Gene_list, Cell_type)

#do with no IQR
Cell_type = filter(Total_Data, SENvCON == "Yes")
Cell_type = filter(Cell_type, Treatment == "none")
Cell_type = filter(Cell_type, Disease == "none")
Cell_type = filter(Cell_type, Control_type == "Prolif")
Cell_type = filter(Cell_type, Sen_type == "DDIS" | Sen_type == "OIS")
Cell_type = filter(Cell_type, (numeric_time > 4 & numeric_time < 12))
unique(Cell_type$Cell_line)

Cell_list = c("BJ" , "HCA2", "HDF", "HFF", "IMR", "MRC5", "Tig3", "WI38")
Cell_type$Cell_line = factor(Cell_type$Cell_line, levels = Cell_list)
unique(Cell_type$Cell_line)
Sen_list = c("DDIS", "OIS")
Cell_type$Sen_type = factor(Cell_type$Sen_type, levels = Sen_list)

Create_pathway_graphs_cell = function(pathway_data, expression_data){
  for(i in pathway_data){
    Gene_data = filter(expression_data, Gene == i)
    Gene_name = paste0(AnalysisDir, "/", i)
    print(i)
    Gene_graph = ggplot(Gene_data, aes(x= factor(Sen_type), y=LogFC, fill=factor(Cell_line, levels =Cell_list)))+
      stat_boxplot(geom ='errorbar') +
      geom_boxplot() +
      labs(x = "", y = "Log FC", fill = "Cell line", title = i)+
      geom_text(aes(label = Cell_line), position = position_dodge(width =0.8), vjust = -0.5)+
      #facet_wrap(~REP) +
      #facet_grid(cols = vars(REP), scales = "free_x", space = "free") +
      scale_x_discrete("", breaks=factor(Sen_list), drop=TRUE) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), panel.background = element_blank(), 
            strip.background = element_blank(), strip.text = element_blank(), 
            axis.text.y = element_text(size = 14), axis.title=element_text(size=14),
            legend.text=element_text(size=13), legend.title=element_text(size=14),
            axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))
    ggsave(paste0(Gene_name, "_noIQR_CellLine.png"))
  }
}

Create_pathway_graphs_cell(Gene_list, Cell_type)

####TO DO hes1, tgfb1 0-4, 5-11####
SENvCON = read.csv("SENvCON_IQR_filtered.csv", row.names =1)
# SENvCON2 = mutate(SENvCON, Details = paste0(as.numeric(factor(Study)), "_",Sen_type, "_", Cell_line, "_", numeric_time))
# SENvCON2 = mutate(SENvCON2, StudyNo = paste0(as.character(as.numeric(factor(Study)))))

TGFB = filter(SENvCON2, Gene == "TGFB1")
TGFB = filter(SENvCON2, numeric_time < 5)

TGFB = filter(TGFB, Control_specifics != "Q0")
unique(TGFB$Study)

TGFB2 = SENvCON2[which(SENvCON2$StudyNo %in% TGFB$StudyNo),]
TGFB = filter(TGFB2, numeric_time > 4 & numeric_time < 12)
TGFB = filter(TGFB, Control_specifics != "Q0")


TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP >= 10000, "UP", ""))
TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 10000 & .$MegaP >= 1000, "2", Sig_nif))
TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 1000 & .$MegaP >= 100, "3", Sig_nif))
TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 100 & .$MegaP >= 20, "4", Sig_nif))
TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 20 & .$MegaP > -20, "NS", Sig_nif))
TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -20 & .$MegaP > -100, "6", Sig_nif))
TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -100 & .$MegaP > -1000, "7", Sig_nif))
TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -1000 & .$MegaP > -10000, "8", Sig_nif))
TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -10000, "DOWN", Sig_nif))
Signif_levs = c("UP", "2", "3", "4", "NS", "6", "7", "8", "DOWN")
TGFB$Sig_nif = factor(TGFB$Sig_nif, levels = Signif_levs)

TGFB = TGFB[which(TGFB$LogFC != min(TGFB$LogFC)),]
TGFB = TGFB[which(TGFB$Sig_nif != "NS"),]
TGFB = mutate(TGFB, Details = paste0(Sen_type, "_", Cell_line, "_", numeric_time))
Sen_list = factor(TGFB$Details)

colors <- colorRampPalette(c("blue", "white", "red")) (9)
colors = c("#FF0000","#FF3F3F","#FF7F7F","#FFBFBF","#FFFFFF","#BFBFFF","#7F7FFF","#3F3FFF","#0000FF")

Gene_list = c("TGFB1", "ACTA2", "HES1", "HEY1", "PDGFA", "CTGF", "TGFBR1", "COL1A1")

getwd()

Total_Data = mutate(Total_Data, StudyNo = paste0(as.character(as.numeric(factor(Study)))))
Total_Data = mutate(Total_Data, Details = paste0(as.numeric(factor(Study)), "_",Sen_type, "_", Cell_line, "_", numeric_time))

SENvCON = filter(Total_Data, SENvCON == "Yes")
SENvCON = filter(SENvCON, Treatment == "none")
SENvCON = filter(SENvCON, Disease == "none")
SENvCON = filter(SENvCON, Control_type == "Prolif")

SENvCON = mutate(SENvCON, Comparison = paste0(Study, "_", X))

Total_Data = Total_Data[which(!(is.na(Total_Data$Control_time))),]
Total_Data$Control_time = as.numeric(as.character(Total_Data$Control_time))

unique(Total_Data$StudyNo)


setwd("/media/njw262/MASSIVE/Transcriptomic_Analysis/Combined_analysis/NEW/analysis/SENvCON/GENE/Indiv_stud/New")

for(i in Gene_list){
  TGFB = filter(SENvCON, Gene == i)
  TGFB = filter(TGFB, numeric_time < 5)
  #TGFB = filter(TGFB, LogFC != 0)
  
  TGFB = filter(TGFB, Control_specifics != "Q0")
  
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP >= 10000, "UP", ""))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 10000 & .$MegaP >= 1000, "2", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 1000 & .$MegaP >= 100, "3", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 100 & .$MegaP >= 20, "4", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 20 & .$MegaP > -20, "NS", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -20 & .$MegaP > -100, "6", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -100 & .$MegaP > -1000, "7", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -1000 & .$MegaP > -10000, "8", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -10000, "DOWN", Sig_nif))
  Signif_levs = c("UP", "2", "3", "4", "NS", "6", "7", "8", "DOWN")
  TGFB$Sig_nif = factor(TGFB$Sig_nif, levels = Signif_levs)
  
  #TGFB = TGFB[which(TGFB$Sig_nif != "NS"),]
  Sen_list = factor(TGFB$Details)
  TGFB = TGFB %>% mutate(LogFC = if_else(.$LogFC > 5, 5, LogFC))
  TGFB = TGFB %>% mutate(LogFC = if_else(.$LogFC < -5, -5, LogFC))
  
  signifs1 = which(Signif_levs %in% TGFB$Sig_nif)
  colors2 = colors[signifs1]
  ggplot(TGFB, aes(x= factor(Comparison), y=LogFC, fill=Sig_nif))+
    geom_col(width = 0.9, colour="black", size =0.05) + #outlier.shape = NA
    scale_fill_manual(values= colors2)+
    #geom_jitter(width=0.1,alpha=0.2) +
    labs(x = "", y = "Log FC", fill = "Gene KD/KO", title = i)+
    #ylim(NA, 10)+
    #facet_wrap(~REP) +
    #facet_grid(cols = vars(REP), scales = "free_x", space = "free") +
    scale_x_discrete("0-4 days vs Prolif", labels=TGFB$Details, drop=TRUE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 11), panel.background = element_blank(), 
          strip.background = element_blank(), strip.text = element_blank(), 
          axis.text.y = element_text(size = 14), axis.title=element_text(size=14),
          legend.text=element_blank(), legend.title=element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none", axis.line.y = element_line(colour = "black"),
          axis.line.x = element_blank(), plot.title = element_text(hjust = 0.5))
  ggsave(paste0(i, "_All_studs_T1.svg"))
  
  TGFB2 = SENvCON[which(SENvCON$StudyNo %in% TGFB$StudyNo),]
  
  TGFB = filter(TGFB2, Gene == i)
  TGFB = filter(TGFB, numeric_time > 4 & numeric_time < 12)
  #TGFB = filter(TGFB, LogFC != 0)
  
  TGFB = filter(TGFB, Control_specifics != "Q0")
  
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP >= 10000, "UP", ""))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 10000 & .$MegaP >= 1000, "2", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 1000 & .$MegaP >= 100, "3", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 100 & .$MegaP >= 20, "4", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 20 & .$MegaP > -20, "NS", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -20 & .$MegaP > -100, "6", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -100 & .$MegaP > -1000, "7", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -1000 & .$MegaP > -10000, "8", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -10000, "DOWN", Sig_nif))
  Signif_levs = c("UP", "2", "3", "4", "NS", "6", "7", "8", "DOWN")
  TGFB$Sig_nif = factor(TGFB$Sig_nif, levels = Signif_levs)
  
  TGFB = TGFB[which(TGFB$Sig_nif != "NS"),]
  Sen_list = factor(TGFB$Details)
  
  TGFB = TGFB %>% mutate(LogFC = if_else(.$LogFC > 5, 5, LogFC))
  TGFB = TGFB %>% mutate(LogFC = if_else(.$LogFC < -5, -5, LogFC))
  
  signifs1 = which(Signif_levs %in% TGFB$Sig_nif)
  colors2 = colors[signifs1]
  ggplot(TGFB, aes(x= factor(Comparison), y=LogFC, fill=Sig_nif))+
    geom_col(width = 0.9, colour="black", size =0.05) + #outlier.shape = NA
    scale_fill_manual(values= colors2)+
    #geom_jitter(width=0.1,alpha=0.2) +
    labs(x = "", y = "Log FC", fill = "Gene KD/KO", title = i)+
    #ylim(Min_Numb, Max_Numb)+
    #facet_wrap(~REP) +
    #facet_grid(cols = vars(REP), scales = "free_x", space = "free") +
    scale_x_discrete("5-11 days vs Prolif", labels=TGFB$Details, drop=TRUE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, size = 14), panel.background = element_blank(), 
          strip.background = element_blank(), strip.text = element_blank(), 
          axis.text.y = element_text(size = 14), axis.title=element_text(size=14),
          legend.text=element_blank(), legend.title=element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none", axis.line.y = element_line(colour = "black"),
          axis.line.x = element_blank(), plot.title = element_text(hjust = 0.5))
  ggsave(paste0(i, "_All_studs_T2.svg"))
  
}



Tim = filter(Time_comp, StudyNo == "83")
Tim = Tim[which(!(duplicated(Tim$Comparison))),]

Time_comp = filter(Total_Data, Control_time < 5 & Control_time > 0)
Time_comp = filter(Time_comp, SENvCON == "No")
Time_comp = filter(Time_comp, Treatment == "none")
Check = Time_comp[which(!(duplicated(Time_comp$X))),]

Time_comp = filter(Time_comp, Sen_type != "REP")
Time_comp = filter(Time_comp, numeric_time > 4 & numeric_time < 12)
Time_comp = mutate(Time_comp, CompT = paste(numeric_time, "_", Control_time))
Time_comp = filter(Time_comp, Control_time >= 1)

Time_comp = mutate(Time_comp, Details = paste0(StudyNo, "_",Sen_type, "_", Cell_line, "_", numeric_time, "_", Control_time))
#Time_comp = mutate(Time_comp, StudyNo = paste0(as.character(as.numeric(factor(Study)))))
Time_comp = mutate(Time_comp, Comparison = paste0(StudyNo, "_", X))

for(i in Gene_list){
  TGFB = filter(Time_comp, Gene == i)
  
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP >= 10000, "UP", ""))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 10000 & .$MegaP >= 1000, "2", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 1000 & .$MegaP >= 100, "3", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 100 & .$MegaP >= 20, "4", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP < 20 & .$MegaP > -20, "NS", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -20 & .$MegaP > -100, "6", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -100 & .$MegaP > -1000, "7", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -1000 & .$MegaP > -10000, "8", Sig_nif))
  TGFB = TGFB %>% mutate(Sig_nif = if_else(.$MegaP <= -10000, "DOWN", Sig_nif))
  Signif_levs = c("UP", "2", "3", "4", "NS", "6", "7", "8", "DOWN")
  TGFB$Sig_nif = factor(TGFB$Sig_nif, levels = Signif_levs)
  
  Sen_list = factor(TGFB$Details)
  TGFB = TGFB %>% mutate(LogFC = if_else(.$LogFC > 5, 5, LogFC))
  TGFB = TGFB %>% mutate(LogFC = if_else(.$LogFC < -5, -5, LogFC))
  
  signifs1 = which(Signif_levs %in% TGFB$Sig_nif)
  colors2 = colors[signifs1]
  ggplot(TGFB, aes(x= factor(Comparison), y=LogFC, fill=Sig_nif))+
    geom_col(width = 0.9, colour="black", size =0.05) + #outlier.shape = NA
    scale_fill_manual(values= colors2)+
    #geom_jitter(width=0.1,alpha=0.2) +
    labs(x = "", y = "Log FC", fill = "Gene KD/KO", title = i)+
    #ylim(NA, 10)+
    #facet_wrap(~REP) +
    #facet_grid(cols = vars(REP), scales = "free_x", space = "free") +
    scale_x_discrete("5-11 days vs 0-4 days", labels=TGFB$Details, drop=TRUE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 14), panel.background = element_blank(), 
          strip.background = element_blank(), strip.text = element_blank(), 
          axis.text.y = element_text(size = 14), axis.title=element_text(size=14),
          legend.text=element_blank(), legend.title=element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none", axis.line.y = element_line(colour = "black"),
          axis.line.x = element_blank(), plot.title = element_text(hjust = 0.5))
  ggsave(paste0(i, "_All_studs_T2_T1.svg"))
}
