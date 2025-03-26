#install.packages("corrr")
#install.packages("Hmisc")
#install.packages("corrplot")
#install.packages("ggpubr")
#install.packages("factoextra")
#BiocManager::install("RankProd")
#For gmp library, Need MP C library. This solved it:
#sudo apt-get install libgmp-dev
#install.packages("gmp")
#For Rmpfr library, Need C library. This solved it:
#sudo apt-get install libmpfr-dev
#install.packages("Rmpfr")
#install.packages("Rmpfr", dependencies=TRUE, repos='http://cran.rstudio.com/')
#BiocManager::install("progeny")


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

#RS
study_list = list.dirs("/Volumes/Bex5TB/Sen23/AwiaitngResponseBeforeESET", full.names = FALSE, recursive = FALSE)
print(study_list)


##UPDATE 090722 JW

Workdir = setwd("/Volumes/Bex5TB/Sen23/AwiaitngResponseBeforeESET")
Workdir
LogFC_Folder = paste0(Workdir, "/LogFC") 
MegaP_Folder = paste0(Workdir, "/MegaP") 
Pvalues_Folder = paste0(Workdir, "/Pvalues") 

fnames <- list.files(LogFC_Folder)
studynames = c()
for(i in fnames){
  study = gsub("_Comparison_limma_log.csv", "", i)
  #names(bob) = gsub("_[^_]+$", "", names(bob))
  studynames = append(studynames, study)
}
an = "Victorelli"
LOG = read.csv(paste0(LogFC_Folder, "/", an, "_Comparison_limma_log.csv"))
LOG = pivot_longer(LOG, `X5_8S_rRNA`:`ZZZ3`, names_to = "Gene", values_to = "LogFC")
Pvalues = read.csv(paste0(Pvalues_Folder, "/", an, "_Comparison_limma_Pvalues.csv"))
Pvalues = pivot_longer(Pvalues, `X5_8S_rRNA`:`ZZZ3`, names_to = "Gene", values_to = "Pvalues")
MegaP = read.csv(paste0(MegaP_Folder, "/", an, "_Comparison_limma_MegP.csv"))
MegaP = pivot_longer(MegaP, `X5_8S_rRNA`:`ZZZ3`, names_to = "Gene", values_to = "MegaP")
Total = cbind(LOG, Pvalues = Pvalues$Pvalues, MegaP = MegaP$MegaP)
Total = Total[which(Total$LogFC != "NA"),]

studynames = studynames[2:length(studynames)]

for(i in studynames){
  LOG = read.csv(paste0(LogFC_Folder, "/", i, "_Comparison_limma_log.csv"))
  LOG = pivot_longer(LOG, `X5_8S_rRNA`:`ZZZ3`, names_to = "Gene", values_to = "LogFC")
  Pvalues = read.csv(paste0(Pvalues_Folder, "/", i, "_Comparison_limma_Pvalues.csv"))
  Pvalues = pivot_longer(Pvalues, `X5_8S_rRNA`:`ZZZ3`, names_to = "Gene", values_to = "Pvalues")
  MegaP = read.csv(paste0(MegaP_Folder, "/", i, "_Comparison_limma_MegP.csv"))
  MegaP = pivot_longer(MegaP, `X5_8S_rRNA`:`ZZZ3`, names_to = "Gene", values_to = "MegaP")
  Added = cbind(LOG, Pvalues = Pvalues$Pvalues, MegaP = MegaP$MegaP)
  Added = Added[which(Added$LogFC != "NA"),]
  Total = rbind(Total, Added)
}


#RS: If there is unexpected columns, remove them. e.g. I'm unsure why these below columns appeared
removecolumns = c("V1", "X5S_rRNA", "hsa.mir.1253", "hsa.mir.423", "hsa.mir.8069.1", "snoZ196")
Total = select(Total, -removecolumns)
setwd("/Volumes/Bex5TB/Sen23/NewAnalysis")
getwd()
write.csv(Total, "Victorelli_combined.csv")

numeric_time <- Total$Time_after_induction
numeric_time <- gsub("d", "", numeric_time)
numeric_time <- as.numeric(numeric_time)
is.numeric(numeric_time)
Total2 = cbind(numeric_time = numeric_time, Total)
#RS:when making StudyName_Comparison_info.csv files, for one of them I put 'Y' insted of 'Yes' for SENvCON so have updated it here
#Total3 = Total2 %>%
#  mutate(SENvCON = ifelse(SENvCON == 'Y', 'Yes', SENvCON))
write.csv(Total2, "ToAdd_20Dec23.csv")

print(colnames(Total2))
Previous_Filtered_Dat = read.csv("/Volumes/Bex5TB/Sen23/NewAnalysis/Total_filtered_dat_19Dec23.csv", row.names =1)
#Previous_Filtered_Dat = Previous_Filtered_Dat[,1:(ncol(Previous_Filtered_Dat)-1)]
#to order columns so they match as identical
Total2 = Total2 %>% select(order(names(.)))
Previous_Filtered_Dat = Previous_Filtered_Dat %>% select(order(names(.)))
identical(colnames(Previous_Filtered_Dat), colnames(Total2))

NewData = rbind(Previous_Filtered_Dat, Total2)
#Previous_Filtered_Dat$Immortal_mechanism = gsub("TERT", "hTERT", Previous_Filtered_Dat$Immortal_mechanism)
write.csv(NewData, "Total_filtered_dat_20Dec23.csv")
print("f")
Previous_Filtered_Dat = Previous_Filtered_Dat[which(Previous_Filtered_Dat$Study != "Yang"),]

unique(Previous_Filtered_Dat$Study)

Etop = filter(Previous_Filtered_Dat, Gene_down == "Tp53")
Etop = filter(Etop, Organ != "Skin")

#######################################################################################################################

m = read.csv("/media/njw262/MASSIVE/Transcriptomic_Analysis/Combined_analysis/vizioli_2_Comparison_limma_Pvalues.csv")
n = m[, 1:50]

##LP once all the datasheets have been created I copied them all to the same folder called Combined_info
#One for LogFC, P_values, and Mega P

#Make sure the folder only has the correct docs in
setwd("/media/njw262/MASSIVE/Transcriptomic_Analysis/Combined_analysis/Pvalues")
fnames <- list.files()
csv <- lapply(fnames, read.csv)
result <- do.call(rbind, csv)
result <- as.data.frame(result)

setwd("/media/njw262/MASSIVE/Transcriptomic_Analysis/Combined_analysis/MegaP")
fnames <- list.files()
csv2 <- lapply(fnames, read.csv)
result2 <- do.call(rbind, csv2)
result2 <- as.data.frame(result2)

result3 = rbind(result, result2)

# setwd("/media/njw262/MASSIVE/Transcriptomic_Analysis/Combined_analysis/RNAseq2/MegaP")
# fnames <- list.files()
# csv3 <- lapply(fnames, read.csv)
# result4 <- do.call(rbind, csv3)
# result4 <- as.data.frame(result4)
# 
# result5 = rbind(result3, result4)

write.csv(result, "/media/njw262/MASSIVE/Transcriptomic_Analysis/Combined_analysis/Combined_Pvalues_Comparison_info.csv")

MegaP = read.csv("/media/njw262/MASSIVE/Transcriptomic_Analysis/Combined_analysis/Initial_combined_MegaP_Comparison_info.csv", row.names = 1)
Pvalues = pivot_longer(result, `X5_8S_rRNA`:`ZZZ3`, names_to = "Gene", values_to = "Pvalues")

write.csv(Pvalues, "/media/njw262/MASSIVE/Transcriptomic_Analysis/Combined_analysis/Pivoted_Pvalues_Comparison_info.csv")

levels(MegaP$Study)

setwd("/home/njw262/Documents/Combined_info/Analysis")
write.csv(result, "Allgenes_NormExp.csv")
#start from import
setwd("/home/njw262/Documents/Combined_info/Analysis")
result = read.csv("Allgenes_MegaP.csv", row.names = 1)

#Combined_data_info = result[, 1:35]
#write.csv(Combined_data_info, "Combined_data_info.csv")



####################################################################################################################################
#IF the data show that some docs have different colnames and cannot be rbinded
#Then this code is used to identify which csv/studies are the problem
Addlist =c()
ColnamesData = colnames(csv2[[1]])

head(ColnamesData)
head(ColnamesData3)

for(i in csv2){
  comparison = colnames(i)
  if(identical(comparison, ColnamesData) == TRUE){
    Addlist = append(Addlist, "One")
  } else{
    Addlist = append(Addlist, "Two")
  }
}
print(Addlist)

#Then once you have identified the study you can use the following code to identify the problems within the study
ColnamesData2 = colnames(csv[[1]])
ColnamesData3 = colnames(csv2[[15]])

Addlist2 = c()
for(i in seq_along(ColnamesData3)){
  if(identical(ColnamesData3[i], ColnamesData[i]) == TRUE){
    Addlist2 = append(Addlist2, "One")
  } else{
    Addlist2 = append(Addlist2, "Two")
  }
}
Addlist2_mat = cbind(ColnamesData3, Addlist2, ColnamesData)

print(Addlist2)
identical(colnames(csv[[1]]), colnames(csv[[3]]))
checkdata = cbind(colnames(csv2[[1]]), colnames(csv2[[15]]))
georg = colnames(csv2[[15]])
georg = georg[1:20]
print(georg)
georg = read.csv("/media/njw262/MASSIVE/Transcriptomic_Analysis/RNAseq_Analysis3/georgilis_1_data/analysis/georgilis_1_Comparison_limma_log.csv")
georg = colnames(georg)
###################################################################################################################################
##LP the rest is all about analysing it. 

#The following filters the dataset according to the conditions of interest

#For only senescent cells that have not undergone treatment
SENvsCON <- filter(NewData, SENvCON == "Yes", Treatment == "none", Control_condition == "Prolif", Disease == "none")
Y <- seq(1:nrow(SENvsCON))

#head(colnames(result), 20)

#You can also do further filtering, and naming according to need
SENvsCON2 <- mutate(SENvsCON, Fullname = paste(Sen_type, Time_after_induction, Cell_line, Study, Y, sep = "_"))

numeric_time <- SENvsCON2$Time_after_induction
print(numeric_time)
print(SENvsCON$Sen_type)
numeric_time <- gsub("d", "", numeric_time)
numeric_time <- as.numeric(numeric_time)
is.numeric(numeric_time)
SENvsCON3 <- cbind(SENvsCON2, numeric_time)
SENvsCON3 <- SENvsCON3 %>% arrange(numeric_time)

#Now we create filters for everything of interest
#All cells all stimuli
All = SENvsCON3[c(1:49, 55,56, 62:67),]
NotREP = filter(SENvsCON3, Sen_type != "REP")
Day0_4 = filter(SENvsCON3, Sen_type != "REP", numeric_time <= 4)
Day5_7 = filter(SENvsCON3, Sen_type != "REP", (numeric_time >= 5 & numeric_time <=7))
Day10_12 = filter(SENvsCON3, Sen_type != "REP", (numeric_time >= 10 & numeric_time <= 12))
Day14_15 = filter(SENvsCON3, Sen_type != "REP", (numeric_time >= 14 & numeric_time <= 15))
Day20 = filter(SENvsCON3, Sen_type != "REP", numeric_time == 20)

##OIS only
OIS = filter(SENvsCON3, Sen_type == "OIS")
OIS_Day0_4 = filter(SENvsCON3, Sen_type == "OIS", numeric_time <= 4)
OIS_Day5_7 = filter(SENvsCON3, Sen_type == "OIS", (numeric_time >= 5 & numeric_time <=7))
OIS_Day10_12 = filter(SENvsCON3, Sen_type == "OIS", (numeric_time >= 10 & numeric_time <= 12))
#DDIS only
DDIS = filter(SENvsCON3, Sen_type == "DDIS")
#REP only
REP = filter(SENvsCON3, Sen_type == "REP")
#Celllines only
IMR = filter(SENvsCON3, Cell_line == "IMR")
MRC = filter(SENvsCON3, Cell_line == "MRC")
HFF = filter(SENvsCON3, Cell_line == "HFF")

IMR_OIS = filter(SENvsCON3, Cell_line == "IMR", Sen_type == "OIS")
IMR_DDIS = filter(SENvsCON3, Cell_line == "IMR", Sen_type == "DDIS")
IMR_REP = filter(SENvsCON3, Cell_line == "IMR", Sen_type == "REP")

Filter_list = list(All, NotREP, Day0_4, Day10_12, Day20, Day5_7, OIS, OIS_Day10_12, 
                   OIS_Day0_4, OIS_Day5_7, DDIS, REP, IMR, IMR_REP, IMR_DDIS, IMR_OIS, MRC, HFF)
names(Filter_list) = c("All", "NotREP", "Day0_4", "Day10_12", "Day20", "Day5_7", "OIS", "OIS_Day10_12", 
                       "OIS_Day0_4", "OIS_Day5_7", "DDIS", "REP", "IMR", "IMR_REP", "IMR_DDIS", "IMR_OIS",
                       "MRC", "HFF")

#Now remove the columns without genes and transpose to get genes in rows
#df = dataframe, Desc_cols = number of columns with descriptive data before the genes, 
#Added_cols = number of descriptive columns after the genes, NA_number = 1 + How many NAs are allowed
######
######CHECK THIS THOROUGHLY

RemoveNAs <- function(df, Desc_cols, Added_cols, NA_number) {
  rownames(df) <- df$Fullname
  df <- df[, Desc_cols:(length(df) - Added_cols)]
  df2 <- t(df)
  names <- rownames(df2)
  df2 <- cbind(df2, names)
  M <- as.data.frame(rowSums(is.na(df2)))
  colnames(M) <- c("NaCol")
  df3 <- as.data.frame(df2)
  df3 <- cbind(df3, NaCol = M)
  NA_divis <- (ncol(df3) - Added_cols)/ NA_number
  df3 <- filter(df3, NaCol < (ncol(df3) - Added_cols)/NA_divis)
  rownames(df3) <- df3$names
  df3 <- df3[, 1:(length(df3) -2)]
  df4 <- as.matrix(df3)
  df4 <- data.matrix(df4, rownames.force = T)
  df4 <- apply(df4, MARGIN = 2, as.numeric)
  rownames(df4) <- rownames(df3)
  apply(df4, MARGIN = 2, class)
  return(df4)
}
