library(dplyr)
library(tidyverse)
library(GEOquery)

#To download data:
#SRAexplorer: https://sra-explorer.info/#
#On GEO, get SRP number to put onto SRA explorer. If SRP not on GEO, check ebi.ac.uk/ena (type GSE number in and find study)
#Download files under 'Bash script for downloading FastQ files'
#In terminal, set directory. Then copy and paste from under the bash into terminal
#setting directory to external hard drive via mac: cd /Volumes/Bex5TB/Sen23/GSE214409

#to get metadata:
#Change GSE and remember directory
gse = getGEO(GEO = 'GSE210020', GSEMatrix = TRUE)
gse
metadata = pData(phenoData(gse[[1]]))
setwd("/Volumes/Bex5TB/Sen23/GSE210020")
write.csv(metadata, "GSE210020_metadata.csv")

#If doesn't work for whatever reason, on GEO go to SRA Run Selector
#Select samples wanted (if all select all)
#Download Metadata using button above. This downloads file to Downloads folder. Move to correct directory

a = read.csv("/Users/rebekahscanlon/Downloads/SraRunTable.txt")
write.csv(a, "GSE102537_metadata.csv")