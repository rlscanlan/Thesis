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
library(usethis)
library(devtools)
library(ggbiplot)
library(pcaMethods)
library(fgsea)
library(ggridges)
library(ggrepel)


#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("pcaMethods")

####df####
setwd("/media/c0068011/Bex5TB/Analysis2024")
total_data = read.csv("/media/c0068011/Bex5TB/SenOmic_Project/Working_Database/Total_filtered_dat_24Jan24.csv", row.names =1)
X = total_data[which(!(duplicated(total_data$Comparison))),]

####df filtering####
SENvCON = filter(total_data, SENvCON == "Yes")
SENvCON = filter(SENvCON, Treatment == "none")
SENvCON = filter(SENvCON, Disease == "none")
SENvCON = filter(SENvCON, Control_type == "Prolif")
Y = SENvCON[which(!(duplicated(SENvCON$Comparison))),]

#this counts how many unique variables are in the Sen_type column of Y
n_distinct(Y$Sen_type)
#this counts how many of each variable appears in each column
count_table <- apply(Y, 2, table)
print(count_table)

FourSen_df = filter(SENvCON, Sen_type == "DDIS" | Sen_type == "OIS" | Sen_type == "REP" |Sen_type == "Bystander")
K = FourSen_df[which(!(duplicated(FourSen_df$Comparison))),]

####IQR on total_data####
# Calculate IQR for each gene
total_iqr_df = total_data %>%
  group_by(Gene) %>%
  summarise(IQR_LogFC = IQR(LogFC))

# Merge the IQR information back to the original dataframe
total_IQR <- total_data %>%
  left_join(total_iqr_df, by = "Gene")

# Filter values outside of the IQR
total_IQR_filtered <- total_IQR %>%
  filter(LogFC >= quantile(LogFC, 0.25) - 1.5 * IQR_LogFC,
         LogFC <= quantile(LogFC, 0.75) + 1.5 * IQR_LogFC)

write.csv(total_IQR, "total_IQR.csv")
write.csv(total_IQR_filtered, "total_IQR_filtered.csv")

Xi = total_IQR_filtered[which(!(duplicated(total_IQR_filtered$Comparison))),]
####IQR on SENvCON####
# Calculate IQR for each gene
SENvCON_iqr_df = SENvCON %>%
  group_by(Gene) %>%
  summarise(IQR_LogFC = IQR(LogFC))

# Merge the IQR information back to the original dataframe
SENvCON_IQR <- SENvCON %>%
  left_join(SENvCON_iqr_df, by = "Gene")

# Filter values outside of the IQR
SENvCON_IQR_filtered <- SENvCON_IQR %>%
  filter(LogFC >= quantile(LogFC, 0.25) - 1.5 * IQR_LogFC,
         LogFC <= quantile(LogFC, 0.75) + 1.5 * IQR_LogFC)

write.csv(SENvCON_IQR, "SENvCON_IQR.csv")
write.csv(SENvCON_IQR_filtered, "SENvCON_IQR_filtered.csv")

Yi = SENvCON_IQR_filtered[which(!(duplicated(SENvCON_IQR_filtered$Comparison))),]

####IQR on FourSen_df####
# Calculate IQR for each gene
four_iqr_df = FourSen_df %>%
  group_by(Gene) %>%
  summarise(IQR_LogFC = IQR(LogFC))

# Merge the IQR information back to the original dataframe
four_IQR <- FourSen_df %>%
  left_join(four_iqr_df, by = "Gene")

# Filter values outside of the IQR
four_IQR_filtered <- four_IQR %>%
  filter(LogFC >= quantile(LogFC, 0.25) - 1.5 * IQR_LogFC,
         LogFC <= quantile(LogFC, 0.75) + 1.5 * IQR_LogFC)

write.csv(four_IQR, "FourSen_IQR.csv")
write.csv(four_IQR_filtered, "FourSen_IQR_filtered.csv")

Ki = four_IQR_filtered[which(!(duplicated(four_IQR_filtered$Comparison))),]

####filter senvcon for four sens####
FourSen_df_IQR = filter(SENvCON_IQR_filtered, Sen_type == "DDIS" | Sen_type == "OIS" | Sen_type == "REP" |Sen_type == "Bystander")

#Using this method is better for my analysis than IQR just the FourSen_df imo
####PCA of whole database without IQR####
####using prcomp and na.omit####
a = cbind(total_data[,c(17,22,49)])
a = a %>%
  pivot_wider(names_from = Comparison, values_from = LogFC)

a1 =  a[-c(1)]
a2 = as.data.frame(t(a1))
colnames(a2) = a$Gene

#Sen_type, Comparison
b = cbind(X[,c(4, 36, 49)])
b = b %>%
  mutate(
    Sen_type = replace(Sen_type, is.na(Sen_type), 'None')
  )
#keep b and b1 separate in case anything goes wrong, it is more accessible to resolve
b1 = b
# Set the row names of g1 based on the Comparison column and ensure it's treated as characters
b1$Comparison = as.character(b1$Comparison)
# Set the row names to the Comparison column
rownames(b1) = b1$Comparison 
#Remove the Comparison column from g1
b1$Comparison = NULL

# Make the row names a column in both data frames
a2$RowName <- rownames(a2)
b1$RowName <- rownames(b1)
#Merge the dfs based on the RowName column. doing g1 first menas those columns will be at front of f df
c = merge(b1, a2, by = "RowName", all = TRUE) 
# Set the RowName column as row names in the df
rownames(c) = c$RowName
# Remove the RowName column
c$RowName = NULL

#run prcomp with na.omit
Sen.pca = prcomp(na.omit(c[c(3:27088)], center = TRUE, scale. = TRUE))
#prcomp with na.omit does not work as it removed all the data as each line has at least one missing value

####PCA Method 1####
#using prcomp when na is computed to equal 0

#columns wanted: Gene, LogFC, Comparison
h = cbind(total_data[,c(17,22,49)])
h = h %>%
  pivot_wider(names_from = Comparison, values_from = LogFC)

h1 =  h[-c(1)]
h2 = as.data.frame(t(h1))
colnames(h2) = h$Gene

#Cell_line, numeric_time, Organ, Sen_type, Comparison
g = cbind(X[,c(4, 24, 25, 36, 49)])
g = g %>%
  mutate(
    Sen_type = replace(Sen_type, is.na(Sen_type), 'None'),
    numeric_time = replace(numeric_time, is.na(numeric_time), 'None')
  )
#keep g and g1 separate in case anythign goes wrong, it is more accessible to resolve
g1 = g
# Set the row names of g1 based on the Comparison column
# Ensure it's treated as characters
g1$Comparison = as.character(g1$Comparison)
# Set the row names to the Comparison column
rownames(g1) = g1$Comparison 
#Remove the Comparison column from g1
g1$Comparison = NULL

# Make the row names a column in both data frames
h2$RowName <- rownames(h2)
g1$RowName <- rownames(g1)
#Merge the dfs based on the RowName column. doing g1 first menas those columns will be at front of f df
f = merge(g1, h2, by = "RowName", all = TRUE)  # 'all = TRUE' does a full join
# Set the RowName column as row names in the df
rownames(f) = f$RowName
# Remove the RowName column
f$RowName = NULL

f2 = f
f2[is.na(f2)] = 0

##prcomp when na = 0 ##
#can use f2 because f2 is f but with na = 0
Sen.pca = prcomp(f2[c(5:27090)], center = TRUE, scale. = TRUE)
summary(Sen.pca)
str(Sen.pca)

#with comparison labels. for the group
#Using current p1 you can change the grouping to: Sen_type, cell_line, organ and numeric_time
#with points only
ggbiplot(Sen.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups = f2$Sen_type)+
  labs(title = "PCA plot total_data using prcomp and NA == 0")

####PCA Method 2####
##using NIPALS from pcaMethods##
resNipals5 = pca(f, method="nipals", center=FALSE, nPcs=5)
resNipals15 = pca(f, method="nipals", center=FALSE, nPcs=15)
resNipals = pca(f, method="nipals")

# Check the PCA object
resNipals
resNipals5
resNipals15

summary(resNipals5)
plot(resNipals15)

scores5 = as.data.frame(resNipals5@scores)
scores15 = as.data.frame(resNipals15@scores)
scores2 = as.data.frame(resNipals@scores)

# Extract loadings (contributions of variables to PCs)
loadings5 = as.data.frame(resNipals5@loadings)
loadings15 = as.data.frame(resNipals15@loadings)
loadings2 = as.data.frame(resNipals@loadings)

#add varibales for labelling back into scores
pcaplot5 = cbind(scores5, g)
pcaplot15 = cbind(scores15, g)
pcaplot2 = cbind(scores2, g)

ggplot(pcaplot2, aes(x = PC1, y = PC2, color = Sen_type)) +
  geom_point(size = 3) +
  #labs(title = "PCA plot total_data using NIPALS, nPC = 5", x = "PC1", y = "PC2") +
  #labs(title = "PCA plot total_data using NIPALS, nPC = 15", x = "PC1", y = "PC2") +
  labs(title = "PCA plot total_data using NIPALS, nPC = 2", x = "PC1", y = "PC2") +
  theme_minimal()

####PCA Method 3####
#compute all na values to equal 0
f2 = f
f2[is.na(f2)] = 0

#check computation was successful
any_na = any(is.na(f2))
print(any_na)

resNipals5_zero = pca(f2, method="nipals", center=FALSE, nPcs=5)
resNipals15_zero = pca(f2, method="nipals", center=FALSE, nPcs=15)
resNipals_zero = pca(f2, method="nipals")

# Check the PCA object
resNipals_zero
resNipals5_zero
resNipals15_zero

summary(resNipals5_zero)
plot(resNipals_zero)

scores5 = as.data.frame(resNipals5_zero@scores)
scores15 = as.data.frame(resNipals15_zero@scores)
scores2 = as.data.frame(resNipals_zero@scores)

# Extract loadings (contributions of variables to PCs)
loadings5 = as.data.frame(resNipals5_zero@loadings)
loadings15 = as.data.frame(resNipals15_zero@loadings)
loadings2 = as.data.frame(resNipals_zero@loadings)

#add varibales for labelling back into scores
pcaplot5 = cbind(scores5, g)
pcaplot15 = cbind(scores15, g)
pcaplot2 = cbind(scores2, g)

ggplot(pcaplot5, aes(x = PC1, y = PC2, color = Sen_type)) +
  geom_point(size = 3) +
  labs(title = "PCA plot total_data using NIPALS and NA == 0, nPC = 5", x = "PC1", y = "PC2") +
  #labs(title = "PCA plot total_data using NIPALS and NA == 0, nPC = 15", x = "PC1", y = "PC2") +
  #labs(title = "PCA plot total_data using NIPALS and NA == 0, nPC = 2", x = "PC1", y = "PC2") +
  theme_minimal()

###PCA Method 4####
#on total data non-IQR
a = cbind(total_data[,c(17,22,49)])
a = a %>%
  pivot_wider(names_from = Comparison, values_from = LogFC)

a1 =  a[-c(1)]
a2 = as.data.frame(t(a1))
colnames(a2) = a$Gene

#Sen_type, Comparison
b = cbind(X[,c(4, 36, 49)])
b = b %>%
  mutate(
    Sen_type = replace(Sen_type, is.na(Sen_type), 'None')
  )
#keep b and b1 separate in case anything goes wrong, it is more accessible to resolve
b1 = b
# Set the row names of g1 based on the Comparison column
# Ensure it's treated as characters
b1$Comparison = as.character(b1$Comparison)
# Set the row names to the Comparison column
rownames(b1) = b1$Comparison 
#Remove the Comparison column from g1
b1$Comparison = NULL

# Make the row names a column in both data frames
a2$RowName <- rownames(a2)
b1$RowName <- rownames(b1)
#Merge the dfs based on the RowName column. doing g1 first menas those columns will be at front of f df
c = merge(b1, a2, by = "RowName", all = TRUE)  # 'all = TRUE' does a full join
# Set the RowName column as row names in the df
rownames(c) = c$RowName
# Remove the RowName column
c$RowName = NULL

#rem colswhen there are more than 50% NA values in a column####
Threshold = 0.5  # Set threshold as 50%
new_df = c[, colMeans(is.na(c)) <= threshold]

cols_to_keep = colMeans(is.na(c)) <= threshold

# Create separate data frames
df_kept = c[, cols_to_keep]  
df_removed = c[, !cols_to_keep] 

resNipals5 = pca(df_kept, method="nipals", center=FALSE, nPcs=5)
resNipals15 = pca(df_kept, method="nipals", center=FALSE, nPcs=15)
resNipals = pca(df_kept, method="nipals")

# Check the PCA object
resNipals
resNipals5
resNipals15

summary(resNipals5)
plot(resNipals)

scores5 = as.data.frame(resNipals5@scores)
scores15 = as.data.frame(resNipals15@scores)
scores2 = as.data.frame(resNipals@scores)

#df of just variables 
b2 = b[, 1:2]

#add varibales for labelling back into scores
pcaplot5 = cbind(scores5, b2)
pcaplot15 = cbind(scores15, b2)
pcaplot2 = cbind(scores2, b2)

ggplot(pcaplot5, aes(x = PC1, y = PC2, color = Sen_type)) +
  geom_point(size = 3) +
  #labs(title = "PCA plot total NOIQR using NIPALS, columns removed when NA>50%, nPC = 5", x = "PC1", y = "PC2") +
  #labs(title = "PCA plot total NO IQR using NIPALS, columns removed when NA>50%, nPC = 15", x = "PC1", y = "PC2") +
  labs(title = "PCA plot total NO IQR using NIPALS, columns removed when NA>50%, nPC = 2", x = "PC1", y = "PC2") +
  theme_minimal()

####PCA of whole database with IQR####
# using total_IQR_filtered 
a = cbind(total_IQR_filtered[,c(17,22,49)])
a = a %>%
  pivot_wider(names_from = Comparison, values_from = LogFC)

a1 =  a[-c(1)]
a2 = as.data.frame(t(a1))
colnames(a2) = a$Gene

#Cell_line, numeric_time, Organ, Sen_type, Comparison
b = cbind(Xi[,c(4, 24, 25, 36, 49)])
b = b %>%
  mutate(
    Sen_type = replace(Sen_type, is.na(Sen_type), 'None'),
    numeric_time = replace(numeric_time, is.na(numeric_time), 'None')
  )
#keep b and b1 separate in case anythign goes wrong, it is more accessible to resolve
b1 = b
# Set the row names of g1 based on the Comparison column
# Ensure it's treated as characters
b1$Comparison = as.character(b1$Comparison)
# Set the row names to the Comparison column
rownames(b1) = b1$Comparison 
#Remove the Comparison column from g1
b1$Comparison = NULL

# Make the row names a column in both data frames
a2$RowName <- rownames(a2)
b1$RowName <- rownames(b1)
#Merge the dfs based on the RowName column. doing g1 first menas those columns will be at front of f df
c = merge(b1, a2, by = "RowName", all = TRUE)  # 'all = TRUE' does a full join
# Set the RowName column as row names in the df
rownames(c) = c$RowName
# Remove the RowName column
c$RowName = NULL

#below is PCA Method 2 (lines 331-362)
resNipals5 = pca(c, method="nipals", center=FALSE, nPcs=5)
resNipals15 = pca(c, method="nipals", center=FALSE, nPcs=15)
resNipals = pca(c, method="nipals")

# Check the PCA object
resNipals
resNipals5
resNipals15

summary(resNipals5)
plot(resNipals)

scores5 = as.data.frame(resNipals5@scores)
scores15 = as.data.frame(resNipals15@scores)
scores2 = as.data.frame(resNipals@scores)

# Extract loadings (contributions of variables to PCs)
loadings5 = as.data.frame(resNipals5@loadings)
loadings15 = as.data.frame(resNipals15@loadings)
loadings2 = as.data.frame(resNipals@loadings)

#add varibales for labelling back into scores
pcaplot5 = cbind(scores5, b)
pcaplot15 = cbind(scores15, b)
pcaplot2 = cbind(scores2, b)

ggplot(pcaplot15, aes(x = PC1, y = PC2, color = Sen_type)) +
  geom_point(size = 3) +
  #labs(title = "PCA plot total_data using NIPALS, IQR, nPC = 5", x = "PC1", y = "PC2") +
  labs(title = "PCA plot total_data using NIPALS, IQR, nPC = 15", x = "PC1", y = "PC2") +
  #labs(title = "PCA plot total_data using NIPALS, IQR, nPC = 2", x = "PC1", y = "PC2") +
  theme_minimal()

##NIPALS but make all NAs 0 before using method - PCA Method 3 (lines 365-402)##
c2 = c
c2[is.na(c2)] = 0

any_na = any(is.na(c2))
print(any_na)

resNipals5_zero = pca(c2, method="nipals", center=FALSE, nPcs=5)
resNipals15_zero = pca(c2, method="nipals", center=FALSE, nPcs=15)
resNipals_zero = pca(c2, method="nipals")

# Check the PCA object
resNipals_zero
resNipals5_zero
resNipals15_zero

summary(resNipals5_zero)
plot(resNipals15_zero)

scores5 = as.data.frame(resNipals5_zero@scores)
scores15 = as.data.frame(resNipals15_zero@scores)
scores2 = as.data.frame(resNipals_zero@scores)

# Extract loadings (contributions of variables to PCs)
loadings5 = as.data.frame(resNipals5_zero@loadings)
loadings15 = as.data.frame(resNipals15_zero@loadings)
loadings2 = as.data.frame(resNipals_zero@loadings)

#add varibales for labelling back into scores
pcaplot5 = cbind(scores5, b)
pcaplot15 = cbind(scores15, b)
pcaplot2 = cbind(scores2, b)

ggplot(pcaplot15, aes(x = PC1, y = PC2, color = Sen_type)) +
  geom_point(size = 3) +
  #labs(title = "PCA plot total_data using NIPALS, IQR, and NA == 0, nPC = 5", x = "PC1", y = "PC2") +
  labs(title = "PCA plot total_data using NIPALS, IQR, and NA == 0, nPC = 15", x = "PC1", y = "PC2") +
  #labs(title = "PCA plot total_data using NIPALS, IQR, and NA == 0, nPC = 2", x = "PC1", y = "PC2") +
  theme_minimal()

##prcomp when na = 0 ##
#can use c2 because c2 is c but with na = 0
#Sen.pca = prcomp(c2[c(5:26381)], center = TRUE, scale. = TRUE)
#summary(Sen.pca)
#str(Sen.pca)

#with comparison labels. for the group
#Using current p1 you can change the grouping to: Sen_type, cell_line, organ and numeric_time
#with points only
#ggbiplot(Sen.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups = c2$Sen_type)+
#  labs(title = "PCA plot total_data using prcomp, IQR, and NA == 0")


####PCA of SenOmic filtered to paper filters####
#-labelled by sen_type, cell_line, organ, sex, numeric_time
unique(FourSen_df_IQR$Cell_line)
#add Sex column to which DF
DF = FourSen_df_IQR %>%
  mutate(Sex = case_when(
    Cell_line %in% c("BJ", "HFF", "MRC", "MRC5", "Tig3", "HDF", "HFL1") ~ "Male",
    Cell_line %in% c("HCA2", "IMR", "LF1", "WI38") ~ "Female",
    Cell_line %in% c("CAF", "FL2") ~ "N/A"
  ))
write.csv(DF, "FourSen_IQR_Sex.csv")
DF = read.csv("FourSen_IQR_Sex.csv", row.names =1)

DFi = DF[which(!(duplicated(DF$Comparison))),]

sum(is.na(l2))

l = cbind(DF[,c(17,22,49)])
l = l %>%
  pivot_wider(names_from = Comparison, values_from = LogFC)

l1 =  l[-c(1)]
l2 = as.data.frame(t(l1))
colnames(l2) = l$Gene

#Cell_line, numeric_time, Organ, Sen_type, Comparison, Sex
m = cbind(DFi[,c(4, 24, 25, 36, 49, 51)])
m = m %>%
  mutate(
    Sen_type = replace(Sen_type, is.na(Sen_type), 'None'),
    numeric_time = replace(numeric_time, is.na(numeric_time), 'None')
  )
#keep b and b1 separate in case anythign goes wrong, it is more accessible to resolve
m1 = m
# Set the row names of g1 based on the Comparison column
# Ensure it's treated as characters
m1$Comparison = as.character(m1$Comparison)
# Set the row names to the Comparison column
rownames(m1) = m1$Comparison 
#Remove the Comparison column from g1
m1$Comparison = NULL

# Make the row names a column in both data frames
l2$RowName <- rownames(l2)
m1$RowName <- rownames(m1)
#Merge the dfs based on the RowName column. doing g1 first menas those columns will be at front of f df
n = merge(m1, l2, by = "RowName", all = TRUE)  # 'all = TRUE' does a full join
# Set the RowName column as row names in the df
rownames(n) = n$RowName
# Remove the RowName column
n$RowName = NULL


resNipals5 = pca(n, method="nipals", center=FALSE, nPcs=5)
resNipals15 = pca(n, method="nipals", center=FALSE, nPcs=15)
resNipals = pca(n, method="nipals")

# Check the PCA object
resNipals
resNipals5
resNipals15

summary(resNipals5)
plot(resNipals)

scores5 = as.data.frame(resNipals5@scores)
scores15 = as.data.frame(resNipals15@scores)
scores2 = as.data.frame(resNipals@scores)

# Extract loadings (contributions of variables to PCs) if plotting
#loadings5 = as.data.frame(resNipals5@loadings)
#loadings15 = as.data.frame(resNipals15@loadings)
#loadings2 = as.data.frame(resNipals@loadings)

#add varibales for labelling back into scores
pcaplot5 = cbind(scores5, m)
pcaplot15 = cbind(scores15, m)
pcaplot2 = cbind(scores2, m)

ggplot(pcaplot15, aes(x = PC1, y = PC2, color = Sen_type)) +
  geom_point(size = 3) +
  #labs(title = "PCA plot foursen using NIPALS, nPC = 5", x = "PC1", y = "PC2") +
  labs(title = "PCA plot foursen using NIPALS, nPC = 15", x = "PC1", y = "PC2") +
  #labs(title = "PCA plot foursen using NIPALS, nPC = 2", x = "PC1", y = "PC2") +
  theme_minimal()
  #geom_text_repel(aes(label = rownames(pcaplot2)), size = 3, max.overlaps = 500)#max.overlaps was to see where bystander was on the plot
  #coord_cartesian(xlim = c(-0.01, 0.01), ylim = c(-0.01, 0.01))  # Set x and y axis limits wihtout removing data outside axis limits
  #xlim(-100, 100) +  # Set x-axis limits
  #ylim(-500, 500)      # Set y-axis limits

  
##NIPALS but make all NAs 0 before using method##
f2 = f
f2[is.na(f2)] = 0

any_na = any(is.na(f2))
print(any_na)

resNipals5_zero = pca(f2, method="nipals", center=FALSE, nPcs=5)
resNipals15_zero = pca(f2, method="nipals", center=FALSE, nPcs=15)
resNipals_zero = pca(f2, method="nipals")

# Check the PCA object
resNipals_zero
resNipals5_zero
resNipals15_zero

summary(resNipals5_zero)
plot(resNipals_zero)

scores5 = as.data.frame(resNipals5_zero@scores)
scores15 = as.data.frame(resNipals15_zero@scores)
scores2 = as.data.frame(resNipals_zero@scores)

# Extract loadings (contributions of variables to PCs)
loadings5 = as.data.frame(resNipals5_zero@loadings)
loadings15 = as.data.frame(resNipals15_zero@loadings)
loadings2 = as.data.frame(resNipals_zero@loadings)

#add varibales for labelling back into scores
pcaplot5 = cbind(scores5, g)
pcaplot15 = cbind(scores15, g)
pcaplot2 = cbind(scores2, g)

ggplot(pcaplot5, aes(x = PC1, y = PC2, color = Sen_type)) +
  geom_point(size = 3) +
  labs(title = "PCA plot total_data using NIPALS and NA == 0, nPC = 5", x = "PC1", y = "PC2") +
  #labs(title = "PCA plot total_data using NIPALS and NA == 0, nPC = 15", x = "PC1", y = "PC2") +
  #labs(title = "PCA plot total_data using NIPALS and NA == 0, nPC = 2", x = "PC1", y = "PC2") +
  theme_minimal()


##prcomp when na = 0 ##
#can use f2 because f2 is f but with na = 0
Sen.pca = prcomp(f2[c(5:27090)], center = TRUE, scale. = TRUE)
summary(Sen.pca)
str(Sen.pca)

#with comparison labels. for the group
#Using current p1 you can change the grouping to: Sen_type, cell_line, organ and numeric_time
#with points only
ggbiplot(Sen.pca, obs.scale = 1, var.scale = 1, var.axes = FALSE, groups = f2$Sen_type)+
  labs(title = "PCA plot total_data using prcomp and NA == 0")




o = n
o[is.na(o)] = 0








####cell type and organ analysis####
####Venn of organ specific sig genes####

Get_Organ_SPECIFIC = function(Datframe, Name){
  Organ_data = filter(Datframe, Organ == Name)
  return(Organ_data)
}

Organ_Names = c("Lung", "Skin")

for(i in Organ_Names){
  assign(paste0(i, "_data"), Get_Organ_SPECIFIC(DF, i))
}

Lung_data = filter(DF, Organ == "Lung")
Skin_data = filter(DF, Organ == "Skin")


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

Lung_GeneSig = Get_medians(Lung_data, "Lung")
Skin_GeneSig = Get_medians(Skin_data, "Skin")


Gene_list = c(as.character(Lung_GeneSig$Gene), as.character(Skin_GeneSig$Gene))
Gene_list = unique(Gene_list)
print("e")
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

Lung_GeneSig = Add_extra_genes(Lung_GeneSig)
Skin_GeneSig = Add_extra_genes(Skin_GeneSig)


##IF using the MegaP values from the initial dataset
Get_MegaP = function(Datframe, Name){
  # Pos_neg = Datframe %>% mutate(Pos_Neg = if_else(.$MedianLog < 0, -1, 1))
  # Pos_neg = mutate(Pos_neg, MedianMegaP = Pos_Neg/MedianPval)
  # Pos_neg = Pos_neg[,c(1:4)]
  write.csv(Datframe, paste0(Name, "_SENvCON_medians.csv"))
  return(Datframe)
}

Lung_GeneSig = Get_MegaP(Lung_GeneSig, "Lung")
Skin_GeneSig = Get_MegaP(Skin_GeneSig, "Skin")


Lung_GeneSig$MedianMegaP = as.numeric(as.character(Lung_GeneSig$MedianMegaP))
Skin_GeneSig$MedianMegaP = as.numeric(as.character(Skin_GeneSig$MedianMegaP))

SigLungGenes = Lung_GeneSig %>%
  mutate(Significant_MegaP = case_when(
    MedianMegaP <= -20 ~ 'Yes',
    MedianMegaP >= 20 ~ 'Yes',
    TRUE ~ 'No'
  ))
SigLungGenes = SigLungGenes[c("Gene", "Significant_MegaP")]

SigSkinGenes = Skin_GeneSig %>%
  mutate(Significant_MegaP = case_when(
    MedianMegaP <= -20 ~ 'Yes',
    MedianMegaP >= 20 ~ 'Yes',
    TRUE ~ 'No'
  ))
SigSkinGenes = SigSkinGenes[c("Gene", "Significant_MegaP")]

#will need checking when running
OrganTypeComparison = cbind(SigLungGenes, SigSkinGenes$Significant_MegaP)
#will need updatign when above
colnames(OrganTypeComparison) = c("Gene", "Lung", "Skin")

write.csv(OrganTypeComparison, "SigGenes_Organ_YesNo_Matrix_MegaP.csv")

Gene_listSIG = unique(OrganTypeComparison$Gene)

Lung = c()
Skin = c()
Lung_Skin = c()

for(i in Gene_list){
  w = filter(OrganTypeComparison, Gene == i)
  if(w[2] == "Yes" & w[3] == "No"){
    Lung = append(Lung, i)
  } else if(w[2] == "No" & w[3] == "Yes"){
    Skin = append(Skin, i)
  } else if(w[2] == "Yes" & w[3] == "Yes"){
    Lung_Skin = append(Lung_Skin, i)}
}


area1 = length(Lung) + length(Lung_Skin)
area2 = length(Skin) + length(Lung_Skin)
area3 = length(Lung_Skin)


grid.newpage()
v = draw.pairwise.venn(area1 = area1,
                       area2 = area2,
                       cross.area = area3,
                       category = c("Lung","Skin"),
                       col= "Black",
                       fill=c("#1434A4", "#E30B5C"),
                       cex = 2.5,
                       cat.cex = 2,
                       margin = 0.1,
                       cat.pos = c(-90, 90),
                       cat.dist = 0.05)
#save file
print("w")
#### GSEA of organ specific sig genes ####
term2gene = read.gmt("/media/c0068011/Bex5TB/SenOmic_Project/SenAnalysisFiles/h.all.v2023.2.Hs.symbols.gmt")
#mac
#term2gene = read.gmt("/Volumes/Bex5TB/SenOmic_Project/SenAnalysisFiles/h.all.v2023.2.Hs.symbols.gmt")

term2gene$term = gsub("HALLMARK_", "", term2gene$term)
grch38 = grch38
#below dfs needed for below analysis
#Top_Genes = read.csv("SigGenes_SENvCON_filteredDirec.csv", row.names =1)
#for each run of the GSEA, need to change the Group, file names and the rows selected on row 445

Top_Genes5 = filter(Top_Genes, Group == "Lung" | Group == "Skin")
Top_Genes5 = filter(Top_Genes5, Test == "LogFC")
Top_Genes5 = Top_Genes5 %>%
  filter(Group == "Lung_Skin") %>%
  distinct() %>%
  bind_rows(Top_Genes5 %>% filter(Group != "Lung_Skin"))
Top_Genes5 = mutate(Top_Genes5, Median = rowMedians(as.matrix(Top_Genes5[,c(3:4)])))

write.csv(Top_Genes5, "Top_Genes5_Organ.csv")

genelist = pull(Top_Genes5, Median)
names(genelist) = pull(Top_Genes5, Gene)
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

write.csv(gsea_all_df, paste0("Skin", "_all_gsea.csv"))
write.csv(gsea_sig_df, paste0("Skin", "_sig_gsea.csv"))
write.csv(gsea_0.1_df, paste0("Skin", "_0.1_gsea.csv"))

#remember to change file name depending on input!!
if(nrow(gsea_all_df) > 0){
  # Dotplot
  g1 = dotplot(gsea_all, color = "pvalue",  showCategory=nrow(gsea_all), split=".sign") + facet_grid(.~.sign)+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(text = element_text(size = 16),
          plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("Skin", "_all_dotplot.pdf"), g1, width = 10, height = 12)
  # Ridge plot
  g3 = ridgeplot(gsea_all, fill = "pvalue") + labs(x = "enrichment distribution")+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("Skin", "_all_ridgeplot.pdf"), width = 10, height = 12)
}

if(nrow(gsea_sig_df) > 0){
  ## Dotplot
  g1 = dotplot(gsea_sig, color = "pvalue",  showCategory=nrow(gsea_sig), split=".sign") + facet_grid(.~.sign)+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(text = element_text(size = 16),
          plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("Lung", "_sig_dotplot.pdf"), g1, width = 10, height = 12)
  ## Ridge plot
  g3 = ridgeplot(gsea_sig, fill = "pvalue") + labs(x = "enrichment distribution")+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("Lung","_sig_ridgeplot.pdf"), g3, width = 10, height = 12)
}

if(nrow(gsea_0.1_df) > 0){
  ## Dotplot
  g1 = dotplot(gsea_0.1, color = "pvalue",  showCategory=nrow(gsea_0.1), split=".sign") + facet_grid(.~.sign)+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(text = element_text(size = 16),
          plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("Lung", "_sig0.1_dotplot.pdf"), g1, width = 10, height = 12)
  ## Ridge plot
  g3 = ridgeplot(gsea_0.1, fill = "pvalue") + labs(x = "enrichment distribution")+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("Lung","_sig0.1_ridgeplot.pdf"), g3, width = 10, height = 12)
}

#### heatmap top 25 up and top 25 down for each organ group (S, L)####
#this only gets the top genes for the shared genes the 1294 in the centre
Get_Top_genes = function(Datframe, Name){
  Top_genes = Datframe[which(Datframe$Gene %in% Lung_Skin),]
  Top_genes = arrange(Top_genes, Gene)
  colnames(Top_genes) = c("Gene", "Pvalues", "MegaP", "LogFC")
  Top_genes$Pvalues = as.numeric(as.character(Top_genes$Pvalues))
  Top_genes$MegaP = as.numeric(as.character(Top_genes$MegaP))
  Top_genes$LogFC = as.numeric(as.character(Top_genes$LogFC))
  Top_genes = Top_genes %>% pivot_longer(`Pvalues`:`LogFC`, names_to = "Test", values_to = Name)
  return(Top_genes)
}

Top_Lung = Get_Top_genes(Lung_GeneSig, "Lung")
Top_Skin = Get_Top_genes(Skin_GeneSig, "Skin")
Top_Genes = cbind(Top_Lung, Skin = Top_Skin$Skin)
write.csv(Top_Genes, "SigGenes_Organs.csv")
Top_Genes2 = Top_Genes

Top_Genes = mutate(Top_Genes, Group = "Lung_Skin")

Sig_list = lst(Lung, Skin, Lung_Skin)

for(i in 1:length(Sig_list)){
  Gene_list1 = Sig_list[[i]]
  Namy = names(Sig_list[i])
  #Gene_list1 = as.data.frame(Gene_list1)
  Gene_list1 = as.character(Gene_list1)
  #Gene_list1 = unlist(Gene_list1)
  print(Sig_list[i])
  print(names(Sig_list[i]))
  #Gene_list1 = as.character(Gene_list1)
  Next_genes = Lung_GeneSig[which(Lung_GeneSig$Gene %in% Gene_list1),]
  Next_genes = arrange(Next_genes, Gene)
  colnames(Next_genes) = c("Gene", "Pvalues", "MegaP", "LogFC")
  Next_genes$Pvalues = as.numeric(as.character(Next_genes$Pvalues))
  Next_genes$MegaP = as.numeric(as.character(Next_genes$MegaP))
  Next_genes$LogFC = as.numeric(as.character(Next_genes$LogFC))
  Next_Lung = Next_genes %>% pivot_longer(`Pvalues`:`LogFC`, names_to = "Test", values_to = "Lung")
  Next_Lung$Gene = as.character(Next_Lung$Gene)
  Next_Lung = arrange(Next_Lung, Gene)
  
  Next_genes = Skin_GeneSig[which(Skin_GeneSig$Gene %in% Gene_list1),]
  Next_genes = arrange(Next_genes, Gene)
  colnames(Next_genes) = c("Gene", "Pvalues", "MegaP", "LogFC")
  Next_genes$Pvalues = as.numeric(as.character(Next_genes$Pvalues))
  Next_genes$MegaP = as.numeric(as.character(Next_genes$MegaP))
  Next_genes$LogFC = as.numeric(as.character(Next_genes$LogFC))
  Next_Skin = Next_genes %>% pivot_longer(`Pvalues`:`LogFC`, names_to = "Test", values_to = "Skin")
  Next_Skin$Gene = as.character(Next_Skin$Gene)
  Next_Skin = arrange(Next_Skin, Gene)
  
  
  #Next_Genes = cbind(Next_OIS, DDIS = Next_DDIS$DDIS, REP = Next_REP$REP, BYS = Next_BYS$BYS, Group = "O_D_B")
  Next_Genes = cbind(Next_Lung, Skin = Next_Skin$Skin, Group = Namy)
  Top_Genes = rbind(Top_Genes, Next_Genes)
}

write.csv(Top_Genes, "SigGenes_SENvCON_filteredDirec.csv")

Top_Genes3 = filter(Top_Genes, Test == "MegaP")

# Remove all duplicated rows based on all columns
Top_Genes4 = Top_Genes3 %>%
  filter(Group == "Lung_Skin") %>%
  distinct() %>%
  bind_rows(Top_Genes3 %>% filter(Group != "Lung_Skin"))

rownames(Top_Genes4) = Top_Genes4$Gene
Top_Genes4 = Top_Genes4[, c(3:4)]
Top_Genes4 = data.matrix(Top_Genes4, rownames.force = TRUE)
#column_order = c("Lung", "Skin")
#Top_Genes4 = Top_Genes4[, column_order]

#get the top 50 upregulated and top 50 downregulated genes to plot on heatmap (1294 is too many to heatmap!!)
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

colors <- colorRampPalette(c("blue", "white", "red")) (9)
breaks <- c(-100000, -10000, -1000, -100, -20, 20, 100, 1000, 10000, 100000)

#pdf version
pdftotalname2 = "SigGenes_Organ_heatmap.pdf"
pdf(pdftotalname2)
par(oma=c(4,0,0,0))
heatmap.2(Top_Genes4, dendrogram = "none", Rowv = TRUE, Colv = FALSE, key = FALSE, col = colors, na.color = "black", 
          density.info = "none", scale = "none", trace = "none", cexCol = 1, cexRow = 0.7, breaks = breaks, margins =c(5,10),
          hclustfun = function(x) hclust(x, method = 'centroid'), distfun = function(x) dist(x, method = 'euclidean'))
dev.off()

#png version
pdftotalname2 = "SigGenes_Organ_heatmap.png"
png(pdftotalname2)
par(oma=c(4,0,0,0))
heatmap.2(Top_Genes4, dendrogram = "none", Rowv = TRUE, Colv = FALSE, key = FALSE, col = colors, na.color = "black", 
          density.info = "none", scale = "none", trace = "none", cexCol = 1, cexRow = 0.7, breaks = breaks, 
          hclustfun = function(x) hclust(x, method = 'centroid'), distfun = function(x) dist(x, method = 'euclidean'))
dev.off()
#svg version
pdftotalname2 = "SigGenes_Organ_heatmap.svg"
svg(pdftotalname2)
par(oma=c(4,0,0,0))
heatmap.2(Top_Genes4, dendrogram = "none", Rowv = TRUE, Colv = FALSE, key = FALSE, col = colors, na.color = "black", 
          density.info = "none", scale = "none", trace = "none", cexCol = 1, cexRow = 0.7, breaks = breaks, 
          hclustfun = function(x) hclust(x, method = 'centroid'), distfun = function(x) dist(x, method = 'euclidean'))
dev.off()



####cell line####
genes = unique(Top_Genes$Gene)
gene_skin = filter(Top_Genes, Group == "Skin" | Group == "Lung_Skin")
gene_skin = unique(gene_skin$Gene)
gene_lung = filter(Top_Genes, Group == "Lung" | Group == "Lung_Skin")
gene_lung = unique(gene_lung$Gene)
FourSen28Genes = c("FGF9", "EREG", "DUSP6", "CLDN1", "WLS", "BMP2", "STC1", "PLAT", "SMIM3",
                   "CNN2", "CITED4", "RIPOR3", "TMEM97", "IDH2", "MARCKS", "SKP2", "OLFML3", "TUBB",
                   "CXCL12", "RRM2", "PTN", "GAS6", "KIF20A", "DIAPH3", "ALKAL1",
                   "CCNB2", "NSD2", "CDC20")

# Find the overlap between sen28 and skin v lung

overlap = intersect(gene_skin, FourSen28Genes)
print(overlap)
The28 = DF %>% filter(Gene %in% FourSen28Genes)

#### sex venn ####
#remove FL2 and CAF from df -- doing this because cannot discern sex
Sex_df = filter(DF, Cell_line != "CAF" & Cell_line != "FL2")

setwd("/media/c0068011/Bex5TB/Analysis2024/sex")

Get_Sex_SPECIFIC = function(Datframe, Name){
  Sex_data = filter(Datframe, Sex == Name)
  return(Sex_data)
}

Sex_Names = c("Female", "Male")

for(i in Sex_Names){
  assign(paste0(i, "_data"), Get_Sex_SPECIFIC(Sex_df, i))
}

Female_data = filter(Sex_df, Sex == "Female")
Male_data = filter(Sex_df, Sex == "Male")


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

Female_GeneSig = Get_medians(Female_data, "Female")
Male_GeneSig = Get_medians(Male_data, "Male")

Gene_list = c(as.character(Female_GeneSig$Gene), as.character(Male_GeneSig$Gene))
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

Female_GeneSig = Add_extra_genes(Female_GeneSig)
Male_GeneSig = Add_extra_genes(Male_GeneSig)


##IF using the MegaP values from the initial dataset
Get_MegaP = function(Datframe, Name){
  # Pos_neg = Datframe %>% mutate(Pos_Neg = if_else(.$MedianLog < 0, -1, 1))
  # Pos_neg = mutate(Pos_neg, MedianMegaP = Pos_Neg/MedianPval)
  # Pos_neg = Pos_neg[,c(1:4)]
  write.csv(Datframe, paste0(Name, "_SENvCON_medians.csv"))
  return(Datframe)
}

Female_GeneSig = Get_MegaP(Female_GeneSig, "Female")
Male_GeneSig = Get_MegaP(Male_GeneSig, "Male")


Female_GeneSig$MedianMegaP = as.numeric(as.character(Female_GeneSig$MedianMegaP))
Male_GeneSig$MedianMegaP = as.numeric(as.character(Male_GeneSig$MedianMegaP))

SigFemaleGenes = Female_GeneSig %>%
  mutate(Significant_MegaP = case_when(
    MedianMegaP <= -20 ~ 'Yes',
    MedianMegaP >= 20 ~ 'Yes',
    TRUE ~ 'No'
  ))
SigFemaleGenes = SigFemaleGenes[c("Gene", "Significant_MegaP")]

SigMaleGenes = Male_GeneSig %>%
  mutate(Significant_MegaP = case_when(
    MedianMegaP <= -20 ~ 'Yes',
    MedianMegaP >= 20 ~ 'Yes',
    TRUE ~ 'No'
  ))
SigMaleGenes = SigMaleGenes[c("Gene", "Significant_MegaP")]

#will need checking when running
SexTypeComparison = cbind(SigFemaleGenes, SigMaleGenes$Significant_MegaP)
#will need updatign when above
colnames(SexTypeComparison) = c("Gene", "Female", "Male")

write.csv(SexTypeComparison, "SigGenes_Sex_YesNo_Matrix_MegaP.csv")

Gene_listSIG = unique(SexTypeComparison$Gene)

Female = c()
Male = c()
F_M = c()

for(i in Gene_list){
  w = filter(SexTypeComparison, Gene == i)
  if(w[2] == "Yes" & w[3] == "No"){
    Female = append(Female, i)
  } else if(w[2] == "No" & w[3] == "Yes"){
    Male = append(Male, i)
  } else if(w[2] == "Yes" & w[3] == "Yes"){
    F_M = append(F_M, i)}
}


area1 = length(Female) + length(F_M)
area2 = length(Male) + length(F_M)
area3 = length(F_M)


grid.newpage()
v = draw.pairwise.venn(area1 = area1,
                       area2 = area2,
                       cross.area = area3,
                       category = c("Female","Male"),
                       col= "Black",
                       fill=c("#1434A4", "#E30B5C"),
                       cex = 2.5,
                       cat.cex = 2,
                       margin = 0.1,
                       cat.pos = c(-90, 90),
                       cat.dist = 0.05)


####sex heatmap ####
Get_Top_genes = function(Datframe, Name){
  Top_genes = Datframe[which(Datframe$Gene %in% F_M),]
  Top_genes = arrange(Top_genes, Gene)
  colnames(Top_genes) = c("Gene", "Pvalues", "MegaP", "LogFC")
  Top_genes$Pvalues = as.numeric(as.character(Top_genes$Pvalues))
  Top_genes$MegaP = as.numeric(as.character(Top_genes$MegaP))
  Top_genes$LogFC = as.numeric(as.character(Top_genes$LogFC))
  Top_genes = Top_genes %>% pivot_longer(`Pvalues`:`LogFC`, names_to = "Test", values_to = Name)
  return(Top_genes)
}

Top_Female = Get_Top_genes(Female_GeneSig, "Female")
Top_Male = Get_Top_genes(Male_GeneSig, "Male")
Top_Genes = cbind(Top_Female, Male = Top_Male$Male)
write.csv(Top_Genes, "SigGenes_bothsex.csv")
Top_Genes2 = Top_Genes

Top_Genes = mutate(Top_Genes, Group = "F_M")

Sig_list = lst(Female, Male, F_M)

for(i in 1:length(Sig_list)){
  Gene_list1 = Sig_list[[i]]
  Namy = names(Sig_list[i])
  #Gene_list1 = as.data.frame(Gene_list1)
  Gene_list1 = as.character(Gene_list1)
  #Gene_list1 = unlist(Gene_list1)
  print(Sig_list[i])
  print(names(Sig_list[i]))
  #Gene_list1 = as.character(Gene_list1)
  Next_genes = Female_GeneSig[which(Female_GeneSig$Gene %in% Gene_list1),]
  Next_genes = arrange(Next_genes, Gene)
  colnames(Next_genes) = c("Gene", "Pvalues", "MegaP", "LogFC")
  Next_genes$Pvalues = as.numeric(as.character(Next_genes$Pvalues))
  Next_genes$MegaP = as.numeric(as.character(Next_genes$MegaP))
  Next_genes$LogFC = as.numeric(as.character(Next_genes$LogFC))
  Next_Female = Next_genes %>% pivot_longer(`Pvalues`:`LogFC`, names_to = "Test", values_to = "Female")
  Next_Female$Gene = as.character(Next_Female$Gene)
  Next_Female = arrange(Next_Female, Gene)
  
  Next_genes = Male_GeneSig[which(Male_GeneSig$Gene %in% Gene_list1),]
  Next_genes = arrange(Next_genes, Gene)
  colnames(Next_genes) = c("Gene", "Pvalues", "MegaP", "LogFC")
  Next_genes$Pvalues = as.numeric(as.character(Next_genes$Pvalues))
  Next_genes$MegaP = as.numeric(as.character(Next_genes$MegaP))
  Next_genes$LogFC = as.numeric(as.character(Next_genes$LogFC))
  Next_Male = Next_genes %>% pivot_longer(`Pvalues`:`LogFC`, names_to = "Test", values_to = "Male")
  Next_Male$Gene = as.character(Next_Male$Gene)
  Next_Male = arrange(Next_Male, Gene)
  
  
  #Next_Genes = cbind(Next_OIS, DDIS = Next_DDIS$DDIS, REP = Next_REP$REP, BYS = Next_BYS$BYS, Group = "O_D_B")
  Next_Genes = cbind(Next_Female, Male = Next_Male$Male, Group = Namy)
  Top_Genes = rbind(Top_Genes, Next_Genes)
}

write.csv(Top_Genes, "SigGenes_sex_SENvCON_filteredDirec.csv")

Top_Genes3 = filter(Top_Genes, Test == "MegaP")

# Remove all duplicated rows based on all columns
Top_Genes4 = Top_Genes3 %>%
  distinct()
rownames(Top_Genes4) = Top_Genes4$Gene
#keeps the megap value for female and male
Top_Genes4 = Top_Genes4[, c(3:4)]
Top_Genes4 = data.matrix(Top_Genes4, rownames.force = TRUE)
#column_order = c("Lung", "Skin")
#Top_Genes4 = Top_Genes4[, column_order]

#get the top 50 upregulated and top 50 downregulated genes to plot on heatmap (1294 is too many to heatmap!!)
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

colors <- colorRampPalette(c("blue", "white", "red")) (9)
breaks <- c(-100000, -10000, -1000, -100, -20, 20, 100, 1000, 10000, 100000)

#pdf version
pdftotalname2 = "SigGenes_Organ_heatmap.pdf"
pdf(pdftotalname2)
par(oma=c(4,0,0,0))
heatmap.2(Top_Genes4, dendrogram = "none", Rowv = TRUE, Colv = FALSE, key = FALSE, col = colors, na.color = "black", 
          density.info = "none", scale = "none", trace = "none", cexCol = 1, cexRow = 0.7, breaks = breaks, margins =c(5,10),
          hclustfun = function(x) hclust(x, method = 'centroid'), distfun = function(x) dist(x, method = 'euclidean'))
dev.off()

#png version
pdftotalname2 = "SigGenes_Organ_heatmap.png"
png(pdftotalname2)
par(oma=c(4,0,0,0))
heatmap.2(Top_Genes4, dendrogram = "none", Rowv = TRUE, Colv = FALSE, key = FALSE, col = colors, na.color = "black", 
          density.info = "none", scale = "none", trace = "none", cexCol = 1, cexRow = 0.7, breaks = breaks, 
          hclustfun = function(x) hclust(x, method = 'centroid'), distfun = function(x) dist(x, method = 'euclidean'))
dev.off()
#svg version
pdftotalname2 = "SigGenes_Organ_heatmap.svg"
svg(pdftotalname2)
par(oma=c(4,0,0,0))
heatmap.2(Top_Genes4, dendrogram = "none", Rowv = TRUE, Colv = FALSE, key = FALSE, col = colors, na.color = "black", 
          density.info = "none", scale = "none", trace = "none", cexCol = 1, cexRow = 0.7, breaks = breaks, 
          hclustfun = function(x) hclust(x, method = 'centroid'), distfun = function(x) dist(x, method = 'euclidean'))
dev.off()




####sex gsea ####
term2gene = read.gmt("/media/c0068011/Bex5TB/SenOmic_Project/SenAnalysisFiles/h.all.v2023.2.Hs.symbols.gmt")
#mac
#term2gene = read.gmt("/Volumes/Bex5TB/SenOmic_Project/SenAnalysisFiles/h.all.v2023.2.Hs.symbols.gmt")

term2gene$term = gsub("HALLMARK_", "", term2gene$term)
grch38 = grch38
#below dfs needed for below analysis
#Top_Genes = read.csv("SigGenes_SENvCON_filteredDirec.csv", row.names =1)
#for each run of the GSEA, need to change the Group, file names and the rows selected on row 445

Top_Genes5 = Top_Genes %>%
  distinct()

#Top_Genes6 = filter(Top_Genes5, Group == "F_M")
#Top_Genes6 = filter(Top_Genes5, Group == "F_M" | Group == "Female")
Top_Genes6 = filter(Top_Genes5, Group == "F_M" | Group == "Male")

Top_Genes6 = filter(Top_Genes6, Test == "LogFC")
#Top_Genes6 = Top_Genes6[-1,]

#this is a df of median logfc for all gene sig in the chosen group
Top_Genes7 = mutate(Top_Genes6, Median = rowMedians(as.matrix(Top_Genes6[,c(3:4)])))
#change this file name depending on what is filtered for above
write.csv(Top_Genes7, "Top_Genes7_Sex_femaleONLY.csv")

genelist = pull(Top_Genes7, Median)
names(genelist) = pull(Top_Genes7, Gene)
genelist = sort(genelist, decreasing = TRUE)

gsea_all = GSEA(geneList = genelist,
                exponent = 1,
                minGSSize = 5,
                maxGSSize = 500,
                pvalueCutoff = 1,
                pAdjustMethod = "none",
                TERM2GENE = term2gene,
                by = "fgsea")

gsea_all_df = as.data.frame(gsea_all)

gsea_0.1 = GSEA(geneList = genelist,
                exponent = 1,
                minGSSize = 5,
                maxGSSize = 500,
                pvalueCutoff = 0.1,
                pAdjustMethod = "none",
                TERM2GENE = term2gene,
                by = "fgsea")
gsea_0.1_df = as.data.frame(gsea_0.1)

gsea_sig = GSEA(geneList = genelist,
                exponent = 1,
                minGSSize = 5,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "none",
                TERM2GENE = term2gene,
                by = "fgsea")
gsea_sig_df = as.data.frame(gsea_sig)

write.csv(gsea_all_df, paste0("FemaleONLY", "_all_gsea.csv"))
write.csv(gsea_sig_df, paste0("FemaleONLY", "_sig_gsea.csv"))
write.csv(gsea_0.1_df, paste0("FemaleONLY", "_0.1_gsea.csv"))

#remember to change file name depending on input!!
if(nrow(gsea_all_df) > 0){
  # Dotplot
  g1 = dotplot(gsea_all, color = "pvalue",  showCategory=nrow(gsea_all), split=".sign") + facet_grid(.~.sign)+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(text = element_text(size = 16),
          plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("FemaleONLY", "_all_dotplot.pdf"), g1, width = 10, height = 12)
  # Ridge plot
  g3 = ridgeplot(gsea_all, fill = "pvalue") + labs(x = "enrichment distribution")+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("FemaleONLY", "_all_ridgeplot.pdf"), width = 10, height = 12)
}

if(nrow(gsea_sig_df) > 0){
  ## Dotplot
  g1 = dotplot(gsea_sig, color = "pvalue",  showCategory=nrow(gsea_sig), split=".sign") + facet_grid(.~.sign)+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(text = element_text(size = 16),
          plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("FemaleONLY", "_sig_dotplot.pdf"), g1, width = 10, height = 12)
  ## Ridge plot
  g3 = ridgeplot(gsea_sig, fill = "pvalue") + labs(x = "enrichment distribution")+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("FemaleONLY","_sig_ridgeplot.pdf"), g3, width = 10, height = 12)
}

if(nrow(gsea_0.1_df) > 0){
  ## Dotplot
  g1 = dotplot(gsea_0.1, color = "pvalue",  showCategory=nrow(gsea_0.1), split=".sign") + facet_grid(.~.sign)+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(text = element_text(size = 16),
          plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("FemaleONLY", "_sig0.1_dotplot.pdf"), g1, width = 10, height = 12)
  ## Ridge plot
  g3 = ridgeplot(gsea_0.1, fill = "pvalue") + labs(x = "enrichment distribution")+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("FemaleONLY","_sig0.1_ridgeplot.pdf"), g3, width = 10, height = 12)
}

####primary v secondary analysis####
SENvCON_IQR_filtered = read.csv("SENvCON_IQR_filtered.csv", row.names=1)
unique(SENvCON_IQR_filtered$Sen_type)
#first want to create a primary sen database and a secondary sen database from the SENvCON IQR filtered DF
#define which sentypes to include in each one
#for this analysis to keep it all similar in thesis, primary = ddis, ois, rep
Primary = filter(SENvCON_IQR_filtered, Sen_type == "DDIS" | Sen_type == "OIS" | Sen_type == "REP")
Secondary = filter(SENvCON_IQR_filtered, Sen_type == "Bystander")
#NIS and RNIS not included asseocndary types of sen as these are ectopic overexpression of notch to induce


####temporal venns and heatmaps####
setwd("/media/c0068011/Bex5TB/Analysis2024")
SENvCON_IQR_filtered = read.csv("SENvCON_IQR_filtered.csv")

Yi = SENvCON_IQR_filtered[which(!(duplicated(SENvCON_IQR_filtered$Comparison))),]

Sen2 = filter(SENvCON_IQR_filtered, Sen_type == "DDIS" | Sen_type == "OIS")
Sen2A = Sen2
Sen2 = Sen2 %>%
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
    TRUE ~ NA_character_
  ))

unique(Sen2$numeric_time2)

Sen2i = Sen2[which(!(duplicated(Sen2$Comparison))),]

T1 = filter (Sen2, numeric_time2 == "0-4")
T1i = T1[which(!(duplicated(T1$Comparison))),]
T2 = filter (Sen2, numeric_time2 == "5-7")
T2i = T2[which(!(duplicated(T2$Comparison))),]
T3 = filter (Sen2, numeric_time2 == "8-11")
T3i = T3[which(!(duplicated(T3$Comparison))),]
T4 = filter (Sen2, numeric_time2 == "12-14")
T4i = T4[which(!(duplicated(T4$Comparison))),]
T5 = filter (Sen2, numeric_time2 == "15+")
T5i = T5[which(!(duplicated(T5$Comparison))),]
TNA = filter (Sen2, numeric_time2 == "NA")
TNAi = TNA[which(!(duplicated(TNA$Comparison))),]

counts <- table(TNAi$Sen_type)

# Print the counts
print(counts)

unique(T5$Study)

####t1####
setwd("/media/c0068011/Bex5TB/Analysis2024/chp4")

#Get_Time_SPECIFIC = function(Datframe, Name){
#  Time_data = filter(Datframe, numeric_time2 == Name)
#  return(Time_data)
#}

#Sen_Names = c("DDIS", "OIS")
#T1, days 0-4
#for(i in Sen_Names){
#  assign(paste0(i, "_Time_data"), Get_Time_SPECIFIC(T1, i))
#}
DDIS_Time_data_T1 = filter(T1, Sen_type == "DDIS")
OIS_Time_data_T1 = filter(T1, Sen_type == "OIS")

#T2, days 5-7
#for(i in Sen_Names){
#  assign(paste0(i, "_Time_data"), Get_Time_SPECIFIC(T2, i))
#}
DDIS_Time_data_T2 = filter(T2, Sen_type == "DDIS")
OIS_Time_data_T2 = filter(T2, Sen_type == "OIS")

#T3, days 8-11
#for(i in Sen_Names){
#  assign(paste0(i, "_Time_data"), Get_Time_SPECIFIC(T3, i))
#}
DDIS_Time_data_T3 = filter(T3, Sen_type == "DDIS")
OIS_Time_data_T3 = filter(T3, Sen_type == "OIS")

#T4, days 12-14
#(i in Sen_Names){
#  assign(paste0(i, "_Time_data"), Get_Time_SPECIFIC(T4, i))
#}
DDIS_Time_data_T4 = filter(T4, Sen_type == "DDIS")
OIS_Time_data_T4 = filter(T4, Sen_type == "OIS")

#T5, 15+
#this only ddis so....?
#for(i in Sen_Names){
#  assign(paste0(i, "_Time_data"), Get_Time_SPECIFIC(T5, i))
#}
#DDIS_Time_data = filter(T5, Sen_type == "DDIS")
#OIS_Time_data = filter(T5, Sen_type == "OIS")

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
  return(Init_dat2)
}

#needs doing for each time df
DDIS_GeneSig = Get_medians(DDIS_Time_data_T2, "DDIS")
OIS_GeneSig = Get_medians(OIS_Time_data_T2, "OIS")

Gene_list = c(as.character(DDIS_GeneSig$Gene), as.character(OIS_GeneSig$Gene))
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

DDIS_GeneSig = Add_extra_genes(DDIS_GeneSig)
OIS_GeneSig = Add_extra_genes(OIS_GeneSig)
print("e")
##CHANGE FILE NAME EACH TIME DF CHANGES
##IF using the MegaP values from the initial dataset
Get_MegaP = function(Datframe, Name){
  # Pos_neg = Datframe %>% mutate(Pos_Neg = if_else(.$MedianLog < 0, -1, 1))
  # Pos_neg = mutate(Pos_neg, MedianMegaP = Pos_Neg/MedianPval)
  # Pos_neg = Pos_neg[,c(1:4)]
  write.csv(Datframe, paste0(Name, "TIME2redo_SENvCON_medians.csv"))
  return(Datframe)
}

DDIS_GeneSig = Get_MegaP(DDIS_GeneSig, "DDIS")
OIS_GeneSig = Get_MegaP(OIS_GeneSig, "OIS")

DDIS_T1 = read.csv("/media/c0068011/Bex5TB/Analysis2024/chp4/DDISTIME1_SENvCON_medians.csv", row.names =1)
DDIS_T2 = read.csv("/media/c0068011/Bex5TB/Analysis2024/chp4/DDISTIME2_SENvCON_medians.csv", row.names =1)
DDIS_T3 = read.csv("/media/c0068011/Bex5TB/Analysis2024/chp4/DDISTIME3_SENvCON_medians.csv", row.names =1)
DDIS_T4 = read.csv("/media/c0068011/Bex5TB/Analysis2024/chp4/DDISTIME4_SENvCON_medians.csv", row.names =1)
OIS_T1 = read.csv("/media/c0068011/Bex5TB/Analysis2024/chp4/OISTIME1_SENvCON_medians.csv", row.names =1)
OIS_T2 = read.csv("/media/c0068011/Bex5TB/Analysis2024/chp4/OISTIME2_SENvCON_medians.csv", row.names =1)
OIS_T3 = read.csv("/media/c0068011/Bex5TB/Analysis2024/chp4/OISTIME3_SENvCON_medians.csv", row.names =1)
OIS_T4 = read.csv("/media/c0068011/Bex5TB/Analysis2024/chp4/OISTIME4_SENvCON_medians.csv", row.names =1)

# List of data frame names
df_names <- c("DDIS_T1", "DDIS_T2", "DDIS_T3", "DDIS_T4", "OIS_T1", "OIS_T2", "OIS_T3", "OIS_T4")

# Loop through each data frame name
for (df_name in df_names) {
  df <- get(df_name)
  df$MedianMegaP <- as.numeric(as.character(df$MedianMegaP))
  assign(df_name, df)
}

#use this when running from beginning, but above does it when running the files loaded in
#DDIS_GeneSig$MedianMegaP = as.numeric(as.character(DDIS_GeneSig$MedianMegaP))
#OIS_GeneSig$MedianMegaP = as.numeric(as.character(OIS_GeneSig$MedianMegaP))

#cange T1 to T2 etc when running
SigDDISGenes = DDIS_T4 %>%
  mutate(Significant_MegaP = case_when(
    MedianMegaP <= -20 ~ 'Yes',
    MedianMegaP >= 20 ~ 'Yes',
    TRUE ~ 'No'
  ))
SigDDISGenes= SigDDISGenes[c("Gene", "Significant_MegaP")]

SigOISGenes = OIS_T4 %>%
  mutate(Significant_MegaP = case_when(
    MedianMegaP <= -20 ~ 'Yes',
    MedianMegaP >= 20 ~ 'Yes',
    TRUE ~ 'No'
  ))
SigOISGenes = SigOISGenes[c("Gene", "Significant_MegaP")]

#will need checking when running
TimeTypeComparison_T1 = cbind(SigDDISGenes, SigOISGenes$Significant_MegaP)
TimeTypeComparison_T2 = cbind(SigDDISGenes, SigOISGenes$Significant_MegaP)
TimeTypeComparison_T3 = cbind(SigDDISGenes, SigOISGenes$Significant_MegaP)
TimeTypeComparison_T4 = cbind(SigDDISGenes, SigOISGenes$Significant_MegaP)
#will need updatign when above
colnames(TimeTypeComparison_T4) = c("Gene", "DDIS", "OIS")
#change file write name each run through
write.csv(TimeTypeComparison_T2, "SigGenes_Time2REDO_YesNo_Matrix_MegaP.csv")

#change T for each timepoint
Gene_listSIG = unique(TimeTypeComparison_T4$Gene)

DDIS = c()
OIS = c()
DDIS_OIS = c()
#change input df_Tx
for(i in Gene_listSIG){
  w = filter(TimeTypeComparison_T4, Gene == i)
  if(w[2] == "Yes" & w[3] == "No"){
    DDIS = append(DDIS, i)
  } else if(w[2] == "No" & w[3] == "Yes"){
    OIS = append(OIS, i)
  } else if(w[2] == "Yes" & w[3] == "Yes"){
    DDIS_OIS = append(DDIS_OIS, i)}
}

area1 = length(DDIS) + length(DDIS_OIS)
area2 = length(OIS) + length(DDIS_OIS)
area3 = length(DDIS_OIS)

t1_D_O = DDIS_OIS
#write.csv(t1_D_O, "Time1_DDISandOISsig.csv")
t2_D_O = DDIS_OIS
#write.csv(t2_D_O, "Time2_DDISandOISsig.csv")
t3_D_O = DDIS_OIS
#write.csv(t3_D_O, "Time3_DDISandOISsig.csv")
t4_D_O = DDIS_OIS
#write.csv(t4_D_O, "Time4_DDISandOISsig.csv")

grid.newpage()
v = draw.pairwise.venn(area1 = area1,
                       area2 = area2,
                       cross.area = area3,
                       category = c("DDIS","OIS"),
                       col= "Black",
                       fill=c("#03ac13", "#FFA500"),
                       cex = 1.5,
                       cat.cex = 1.5,
                       margin = 0.1,
                       cat.pos = c(-90, 90),
                       cat.dist = 0.1)
#only do grid.text if want a printed title
#grid.text("Significantly expressed genes 12-14 days
#    post senescence induction", y = unit(0.8, "npc"), gp = gpar(fontsize = 16, fontface = "bold"))
#save file

####have df for each of the sig genes to both ddis and ois####
#compare and find overalps

common_genes <- Reduce(intersect, list(t1_D_O, t2_D_O, t3_D_O, t4_D_O))
print(common_genes)  # Genes shared by all lists
#72 genes across all timepoints
#heatmap of the 72 genes
#merge all 4 databases for the 72 genes

Get_Top_genes = function(Datframe, Name){
  Top_genes = Datframe[which(Datframe$Gene %in% common_genes),]
  Top_genes = arrange(Top_genes, Gene)
  colnames(Top_genes) = c("Gene", "Pvalues", "MegaP", "LogFC")
  Top_genes$Pvalues = as.numeric(as.character(Top_genes$Pvalues))
  Top_genes$MegaP = as.numeric(as.character(Top_genes$MegaP))
  Top_genes$LogFC = as.numeric(as.character(Top_genes$LogFC))
  Top_genes = Top_genes %>% pivot_longer(`Pvalues`:`LogFC`, names_to = "Test", values_to = Name)
  return(Top_genes)
}

Top_DDIS = Get_Top_genes(DDIS_T4, "DDIS")
Top_OIS = Get_Top_genes(OIS_T4, "OIS")
Top_Genes = cbind(Top_DDIS, "Day 12-14 OIS" = Top_OIS$OIS)
colnames(Top_Genes)[colnames(Top_Genes) == "DDIS"] <- "Day 12-14 DDIS"
write.csv(Top_Genes, "TOP72Genes_t4.csv")

Top_GenesT1 = Top_Genes
Top_GenesT2 = Top_Genes
Top_GenesT3 = Top_Genes
Top_GenesT4 = Top_Genes

top72 = Top_GenesT1 %>%
  full_join(Top_GenesT2, by = c("Gene", "Test")) %>%
  full_join(Top_GenesT3, by = c("Gene", "Test")) %>%
  full_join(Top_GenesT4, by = c("Gene", "Test"))

top72_megap = filter(top72, Test == "MegaP")
rownames(top72_megap) = top72_megap$Gene
top72_log = filter(top72, Test == "LogFC")
rownames(top72_log) = top72_log$Gene

top72_megap = top72_megap[, c(3:10)]
top72_megap = data.matrix(top72_megap, rownames.force = TRUE)
top72_megap <- top72_megap[, c(1,3,5,7,2,4,6,8)]

top72_log = top72_log[, c(3:10)]
top72_log = data.matrix(top72_log, rownames.force = TRUE)
top72_log <- top72_log[, c(1,3,5,7,2,4,6,8)]


colours <- colorRampPalette(c("blue", "white", "red")) (11)
breaks <- c(-1000000, -100000, -10000, -1000, -100, -20, 20, 100, 1000, 10000, 100000, 1000000)

#pdf version, needs breaks for megap but not logfc
#key =false for megap, true for logfc
pdftotalname2 = "SigGenes_temporal72_heatmap_megap.pdf"
pdf(pdftotalname2, height = 15, width = 10)
par(oma=c(4,0,0,0))
heatmap.2(top72_megap, dendrogram = "none", Rowv = TRUE, Colv = FALSE, key = FALSE, col = colours, na.color = "black", 
          density.info = "none", scale = "none", trace = "none", cexCol = 1, cexRow = 0.75,
          breaks = breaks,
          hclustfun = function(x) hclust(x, method = 'centroid'), distfun = function(x) dist(x, method = 'euclidean'))
dev.off()

#svg version
#svgtotalname2 = "SigGenes_temporal72_heatmap_log.svg"
#svg(svgtotalname2)
#par(oma=c(4,0,0,0))
#heatmap.2(top72_log, dendrogram = "none", Rowv = TRUE, Colv = FALSE, key = TRUE, col = colors, na.color = "black", 
#          density.info = "none", scale = "none", trace = "none", cexCol = 1, cexRow = 0.7,
#          #breaks = breaks, 
#          hclustfun = function(x) hclust(x, method = 'centroid'), distfun = function(x) dist(x, method = 'euclidean'))
#dev.off()

###gsea of temporal timepoints####
#linux
term2gene = read.gmt("/media/c0068011/Bex5TB/SenOmic_Project/SenAnalysisFiles/h.all.v2023.2.Hs.symbols.gmt")
#mac
#term2gene = read.gmt("/Volumes/Bex5TB/SenOmic_Project/SenAnalysisFiles/h.all.v2023.2.Hs.symbols.gmt")

term2gene$term = gsub("HALLMARK_", "", term2gene$term)
grch38 = grch38
#below dfs needed for below analysis
#Top_Genes = read.csv("SigGenes_SENvCON_filteredDirec.csv", row.names =1)

T1Genelist = t1_D_O
T2Genelist = t2_D_O
T3Genelist = t3_D_O
T4Genelist = t4_D_O

# Rename the columns
colnames(OIS_T2) <- paste(colnames(OIS_T2), "OIS", sep = "_")
colnames(DDIS_T2) <- paste(colnames(DDIS_T2), "DDIS", sep = "_")
# Combine the data frames
T1_common = cbind(DDIS_T1, OIS_T1)
T2_common = cbind(DDIS_T2, OIS_T2)
T3_common = cbind(DDIS_T3, OIS_T3)
T4_common = cbind(DDIS_T4, OIS_T4)
#make gene name the row name
rownames(T1_common) = T1_common$Gene_DDIS

T1_common = T1_common %>%
  filter(Gene_DDIS %in% T1Genelist)
T2_common = T2_common %>%
  filter(Gene_DDIS %in% T2Genelist)
T3_common = T3_common %>%
  filter(Gene_DDIS %in% T3Genelist)
T4_common = T4_common %>%
  filter(Gene_DDIS %in% T4Genelist)

#select columns for logfc
T1_common_log = T1_common[, c(1,4,8)]
T2_common_log = T2_common[, c(1,4,8)]
T3_common_log = T3_common[, c(1,4,8)]
T4_common_log = T4_common[, c(1,4,8)]

#this is a df of median logfc for all gene sig in the chosen group. note  it will median a
#neg and pos so directions may change
T1_common_log = mutate(T1_common_log, Median = rowMedians(as.matrix(T1_common_log[,c(2:3)])))
T2_common_log = mutate(T2_common_log, Median = rowMedians(as.matrix(T2_common_log[,c(2:3)])))
T3_common_log = mutate(T3_common_log, Median = rowMedians(as.matrix(T3_common_log[,c(2:3)])))
T4_common_log = mutate(T4_common_log, Median = rowMedians(as.matrix(T4_common_log[,c(2:3)])))
#change this file name depending on what is filtered for above
write.csv(T4_common_log, "T4_commongenes_log_medians.csv")

#do below for each timepoint, will need to change T1 to T2 for eg where it appears
genelist = pull(T4_common_log, Median)
names(genelist) = pull(T4_common_log, Gene_DDIS)
genelist = sort(genelist, decreasing = TRUE)

gsea_all = GSEA(geneList = genelist,
                exponent = 1,
                minGSSize = 5,
                maxGSSize = 500,
                pvalueCutoff = 1,
                pAdjustMethod = "none",
                TERM2GENE = term2gene,
                by = "fgsea")

gsea_all_df = as.data.frame(gsea_all)

gsea_0.1 = GSEA(geneList = genelist,
                exponent = 1,
                minGSSize = 5,
                maxGSSize = 500,
                pvalueCutoff = 0.1,
                pAdjustMethod = "none",
                TERM2GENE = term2gene,
                by = "fgsea")
gsea_0.1_df = as.data.frame(gsea_0.1)

gsea_sig = GSEA(geneList = genelist,
                exponent = 1,
                minGSSize = 5,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                pAdjustMethod = "none",
                TERM2GENE = term2gene,
                by = "fgsea")
gsea_sig_df = as.data.frame(gsea_sig)

write.csv(gsea_all_df, paste0("T4_log", "_all_gsea.csv"))
write.csv(gsea_sig_df, paste0("T4_log", "_sig_gsea.csv"))
write.csv(gsea_0.1_df, paste0("T4_log", "_0.1_gsea.csv"))

#remember to change file name depending on input!!
if(nrow(gsea_all_df) > 0){
  # Dotplot
  g1 = dotplot(gsea_all, color = "pvalue",  showCategory=nrow(gsea_all), split=".sign") + facet_grid(.~.sign)+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(text = element_text(size = 16),
          plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("T4_Log", "_all_dotplot.pdf"), g1, width = 10, height = 12)
  # Ridge plot
  g3 = ridgeplot(gsea_all, fill = "pvalue") + labs(x = "enrichment distribution")+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("T4_Log", "_all_ridgeplot.pdf"), width = 10, height = 12)
}

if(nrow(gsea_sig_df) > 0){
  ## Dotplot
  g1 = dotplot(gsea_sig, color = "pvalue",  showCategory=nrow(gsea_sig), split=".sign") + facet_grid(.~.sign)+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(text = element_text(size = 16),
          plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("T4_Log", "_sig_dotplot.pdf"), g1, width = 10, height = 12)
  ## Ridge plot
  g3 = ridgeplot(gsea_sig, fill = "pvalue") + labs(x = "enrichment distribution")+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("T4_Log","_sig_ridgeplot.pdf"), g3, width = 10, height = 12)
}

if(nrow(gsea_0.1_df) > 0){
  ## Dotplot
  g1 = dotplot(gsea_0.1, color = "pvalue",  showCategory=nrow(gsea_0.1), split=".sign") + facet_grid(.~.sign)+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(text = element_text(size = 16),
          plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("T4_Log", "_sig0.1_dotplot.pdf"), g1, width = 10, height = 12)
  ## Ridge plot
  g3 = ridgeplot(gsea_0.1, fill = "pvalue") + labs(x = "enrichment distribution")+
    scale_fill_gradient(low = "blue", high = "red")+
    theme(plot.margin = unit(c(1,1,1,1), "cm"))
  ggsave(paste0("T4_Log","_sig0.1_ridgeplot.pdf"), g3, width = 10, height = 12)
}

