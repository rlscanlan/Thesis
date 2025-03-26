library(ggridges)
library(Biobase)
library(tximport)
library(limma)
library(edgeR)
library(dplyr)
library(topGO)
library(org.Hs.eg.db)
library(VennDiagram)
library(KEGGREST) 
library(gplots)
library(RColorBrewer)
library(tidyverse)
library(biomaRt)
library(DESeq2)
library(annotables)
library(pheatmap)
library(ggplot2)
library(gtools)
library(qusage)
library(devtools)
library(clusterProfiler)

study_list = list.dirs("/Volumes/Bex5TB/Sen23/AwiaitngResponseBeforeESET", full.names = FALSE, recursive = FALSE)
print(study_list)

####saveEset####
#This function is used within the main function
#Writes all the outputs to files
SaveEset <- function(k, p4, res_annot2, StudyName) {
  SaveDir = paste("/Volumes/Bex5TB/Sen23/AwiaitngResponseBeforeESET/", StudyName, "/analysis", sep = "")
  kFile = paste(SaveDir, "/k.csv", sep = "")
  p4File = paste(SaveDir, "/p4.csv", sep = "")
  resFile = paste(SaveDir, "/res_annot2.csv", sep = "")
  kFile2 = paste(SaveDir, "/", StudyName, "_k.csv", sep = "")
  p4File2 = paste(SaveDir, "/", StudyName, "_p4.csv", sep = "")
  resFile2 = paste(SaveDir, "/", StudyName, "_res_annot2.csv", sep = "")
  write.csv(k, kFile)
  write.csv(p4, p4File)
  write.csv(res_annot2, resFile)
  write.csv(k, kFile2)
  write.csv(p4, p4File2)
  write.csv(res_annot2, resFile2)
  return(kFile)
}

####Save_data_forDESeq####
#This function creates a DESeq data R file in case we want to run DESeq2 instead of limma in future
#It is included in the function below so doesn't need to be run separately
Save_data_for_DESeq = function(mytxi, p3){
  deseq_data = DESeqDataSetFromTximport(txi = mytxi,
                                        colData = p3,
                                        design = ~ condition)
  
  deseq_data = DESeq(deseq_data)
  deseq_file = paste0(anal_address, StudyName, "_deseq_data.RData")
  save(deseq_data, file = deseq_file)
  
  count_table_raw = counts(deseq_data, normalized = FALSE)
  colnames(count_table_raw) = p3$condition
  count_table_normalized = counts(deseq_data, normalized = TRUE)
  colnames(count_table_normalized) = p3$condition
  
  gene_names = as.data.frame(rownames(count_table_raw))
  colnames(gene_names) = "ensgene"
  gene_names = left_join(gene_names, grch38, by = "ensgene")
  gene_names = as.data.frame(gene_names[!duplicated(gene_names$ensgene),])
  count_table_raw = cbind(count_table_raw, gene_names)
  count_table_normalized = cbind(count_table_normalized, gene_names)
  
  raw_count_name = paste0(anal_address, StudyName, "_count_table_raw.csv")
  write.csv(as.data.frame(count_table_raw), raw_count_name)
  norm_count_name = paste0(anal_address, StudyName, "_count_table_normalized.csv")
  write.csv(as.data.frame(count_table_normalized), norm_count_name)
  
  vst_data = DESeq2::vst(deseq_data, blind = TRUE)
  plot_data = plotPCA(vst_data, intgroup="condition", returnData=TRUE)
  percentVar = round(100 * attr(plot_data, "percentVar"))
  PCA_p = ggplot(plot_data, aes(PC1, PC2, color=condition, label = "label")) +
    geom_point(size=1.5) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
  
  ggsavename = paste0(anal_address, StudyName, "_Deseq_PCA.png")
  ggsave(ggsavename, PCA_p)
}

####create_and_save_eset####

Create_and_save_eset = function(){
  # wd_address = paste0("/media/njw262/MASSIVE/Transcriptomic_Analysis/RNAseq_Analysis4/", StudyName, "/Salmon/quant/")
  # setwd(wd_address)
  #mydir = paste0("/Volumes/Bex5TB/Sen23/test/", StudyName, "/salmon/")
  address = paste0("/Volumes/Bex5TB/Sen23/AwiaitngResponseBeforeESET/", StudyName, "/", StudyName, "_data_info.csv")
  #address = paste0("/Volumes/Bex5TB/Sen23/test/", StudyName, "/reference_folder/", StudyName, "_data_info.csv")
  p <- read.csv(address)
  pop <- rep("TSI", nrow(p))
  centre <- rep("UNIGE", nrow(p))
  p3 <- as.data.frame(cbind(p, pop, centre))
  #p3 = Add_replicate_row(p2)
  p3 <- mutate(p3, Fullname = paste(condition, replicate, sep = "_"))
  run <- paste0(p3$Run, "_quant")
  condition <- p3$condition
  rownames(p3) <- p3$Run
  mydir = paste0("/Volumes/Bex5TB/Sen23/AwiaitngResponseBeforeESET/", StudyName, "/quants")
  #Locate the quant.sf files and create an array of expression for expression eset:
  myfiles <- file.path(mydir, run, "quant.sf")
  names(myfiles) <- p3$Run
  mygene = read.csv("/Volumes/Bex5TB/SenAnalysisFiles/tx2gene_ensembl_v87.csv")
  #mygene <- read_tsv("/media/njw262/DATA/Ageing_review/gencode_v32.gene.map.tsv", col_names = TRUE)
  # anal_address = paste0("/media/njw262/MASSIVE/Transcriptomic_Analysis/RNAseq_Analysis3/",
  #                       StudyName, "/analysis/")
  #NB lengthScaledTPM should only be needed if comparing transcripts across genes
  print(file.exists(myfiles))
  mytxi <- tximport(myfiles, type="salmon", tx2gene=mygene, ignoreTxVersion = TRUE, countsFromAbundance = "scaledTPM")
  #Use this one for DESeq (without TPM normalisation)
  mytxi2 <- tximport(myfiles, type="salmon", tx2gene=mygene, ignoreTxVersion = TRUE)
  
  #SAVE THE DESEQ DATA
  Save_data_for_DESeq(mytxi2, p3)
  
  #Create Eset
  dge <- DGEList(mytxi$counts)
  # filtering
  keep <- filterByExpr(dge)
  dge <- dge[keep, ]
  dge <- calcNormFactors(dge)
  #turn into alterable matrix and prepare col and rownames
  k <- as.data.frame(dge$counts)
  res = as.data.frame(rownames(k))
  colnames(res) = "ensgene"
  res_annot = left_join(res, grch38, by = "ensgene")
  res_annot2 <- res_annot[!duplicated(res_annot$ensgene),]
  res_annot2 <- res_annot2[!duplicated(res_annot2$symbol),]
  res_annot2 = res_annot2[which(res_annot2$ensgene != "NA" & res_annot2$symbol != "NA"),]
  rownames(res_annot2) = res_annot2$ensgene
  k <- k[which(rownames(k) %in% rownames(res_annot2) ==T),]
  raw_k_name = paste0(anal_address, StudyName, "_k_raw.csv")
  write.csv(k, raw_k_name)
  #Then normalise the data with cpm or voom depending on the variance
  totals = c(colSums(k))
  Red_line = max(totals)/ min(totals)
  print(Red_line)
  if(Red_line <= 3){
    dge = cpm(k, log=TRUE, prior.count = 3)
    cpm_name = paste0(anal_address, StudyName, "_cpm_density_plot.png")
    png(cpm_name)
    #png(cpm_name, width = 20, height = 20, units = "in", res = 300)
    #plotDensities(dge, legend = TRUE, col =sample_colours)
    plotDensities(dge, legend = FALSE)
    dev.off()
  } else{
    vdes = as.matrix(p3$replicate)
    dge = voom(dge, vdes, plot = T)
    voom_name = paste0(anal_address, StudyName, "_voom_density_plot.png")
    png(voom_name)
    #png(voom_name, width = 20, height = 20, units = "in", res = 300)
    #plotDensities(dge, legend = TRUE, col =sample_colours)
    plotDensities(dge, legend = FALSE)
    dev.off()
  }
  #Remake k with normalised values
  k <- as.data.frame(dge)
  k <- k[which(rownames(k) %in% rownames(res_annot2) ==T),]
  p4 = p3
  rownames(p4) = p4$Fullname
  krun = colnames(k)
  knames = rownames(p4)
  p4run = as.character(p4$Run)
  #To test you can randomise the vector
  #p5run = sample(p4run)
  true_vector = c()
  #See if all the run names match
  for(i in seq_along(krun)){
    if(krun[i] == p4run[i]){
      true_vector = append(true_vector, TRUE)
    } else {
      true_vector = append(true_vector, FALSE)
    }
  }
  #This tells you if all the values in true_vector are true
  m = all(true_vector)
  if(m == TRUE){
    colnames(k) = knames
  } else {
    print("You have got your run names mixed up!!!!!!")
  }
  
  grch2 <- grch38[!duplicated(grch38$ensgene),]
  grch2 <- grch2[!duplicated(grch2$symbol),]
  grch2 = grch2[which(grch2$ensgene != "NA" & grch2$symbol != "NA"),]
  
  k2 = cbind(k, symbol = res_annot2$symbol)
  rownames(k2) = k2$symbol
  k2 = k2[, 1:length(k)]
  other_genes = grch2[which(!(grch2$symbol %in% rownames(k2))),]
  NA_vec = c(rep(NA, nrow(other_genes)))
  NA_vec = as.data.frame(NA_vec)
  NA_matrix = NA_vec[rep(names(NA_vec), ncol(k2))]
  rownames(NA_matrix) = other_genes$symbol
  colnames(NA_matrix) = colnames(k2)
  total_genes = rbind(k2, NA_matrix)
  total_genes2 = arrange(total_genes, rownames(total_genes))
  total_genes3 = t(total_genes2)
  p5 = mutate(p4, Fullername = paste(study, condition, replicate, sep = "_"))
  combined_data = cbind(p5, total_genes3)
  combined_data2 = combined_data %>% pivot_longer(`5_8S_rRNA`:`ZZZ3`, names_to = "Gene", values_to = "norm_express")
  combined_data_name = paste0(anal_address, StudyName, "_normalised_expression.csv")
  write.csv(combined_data2, combined_data_name)
  Filtered_genes = combined_data2[which(combined_data2$Gene %in% rownames(k2)) ,]
  filtered_data_name = paste0(anal_address, StudyName, "_normalised_expression_NAfiltered.csv")
  write.csv(Filtered_genes, filtered_data_name)
  #Create an expression set
  eset <- ExpressionSet(assayData = as.matrix(k),
                        phenoData = AnnotatedDataFrame(p4),
                        featureData = AnnotatedDataFrame(res_annot2))
  SaveEset(k, p4, res_annot2, StudyName)
  
  MDSname = paste0(anal_address, StudyName, "_limma_MDS.png")
  png(MDSname)
  plotMDS(eset, labels = pData(eset)[, "label"], gene.selection = "common", cex = 0.6)
  dev.off()
  return(eset)
}

#this runs fro just esets
#I'm using it instead of the overall function for everything to remove problem smaples whose curves do not match
#then moving to 'removing irregular curves' tab
for(i in study_list){
  StudyName <- i
  print(StudyName)
  anal_address = paste0("/Volumes/Bex5TB/Sen23/AwiaitngResponseBeforeESET/",
                        StudyName, "/analysis/")
  eset = Create_and_save_eset()
  }


####get_Design####
Get_Design = function(eset){
  Design = model.matrix(~0 + label, data = pData(eset))
  return(Design)
}

####makeConstrastsFromString####
makeContrastsFromString <- function(s){
  eval(parse(text = paste("makeContrasts(", s, ")")))
}

####term2gene####
term2gene = read.gmt("/Volumes/Bex5TB/SenAnalysisFiles/h.all.v2023.2.Hs.symbols.gmt")
term2gene$term = gsub("HALLMARK_", "", term2gene$term)
grch38 = grch38

####GetPandlogfromLimma####
#The below function generates the P values and Log values using the Eset
#It also runs GSEA

GetPandlogfromLimma = function(eset){
  #First we make the linear model
  comp_name = paste0(anal_address, StudyName, "_comparison_info.csv")
  comp_info = read.csv(comp_name)
  comparisons = as.character(comp_info$X)
  print(comparisons)
  comp_info2 = mutate(comp_info,  test = gsub("_.*", "", comparisons))
  comp_info2 = mutate(comp_info2,  control = gsub(".*_", "", comparisons))
  comp_info3 = mutate(comp_info2,  comparison = paste0(X, " = label", test, " - label", control, ",\n"))
  
  cm = toString(comp_info3$comparison)
  cm2 = gsub("\\\n,", "\\\n", cm)
  cm3 = paste(cm2, "levels = Design")
  print(cm3)
  cm <- makeContrastsFromString(cm3)
  
  fit <- lmFit(eset, Design)
  fit2 <- contrasts.fit(fit, contrasts = cm)
  #NB don't use trend=TRUE with voom (just limma trend)
  k_raw = read.csv(paste0(anal_address, StudyName, "_k_raw.csv"), row.names = 1)
  totals = c(colSums(k_raw))
  Red_line = max(totals)/ min(totals)
  if(Red_line <= 3){
    fit2 <- eBayes(fit2, trend = TRUE)
  } else{
    fit2 <- eBayes(fit2)
  }
  results <- decideTests(fit2)
  results <- as.data.frame(results)
  rownames(results) <- fData(eset)$symbol
  write.csv(results, paste0(anal_address, StudyName, "_limma_results.csv"))
  stats <- topTable(fit2, number = nrow(fit2), sort.by = "none")
  # stats <- as.data.frame(stats)
  # write.csv(stats, paste0(anal_address, StudyName, "_limma_stats.csv"))
  
  if(nrow(comp_info) == 1){
    colnom = as.character(comparisons)
    colnum = which(colnames(stats)=="logFC" )
    colnames(stats)[colnum] = colnom
  }
  write.csv(stats, paste0(anal_address, StudyName,  "_limma_stats.csv"))
  LogFC = stats[, which(colnames(stats) %in% comparisons)]
  if(nrow(comp_info) == 1){
    LogFC = as.data.frame(LogFC)
    colnames(LogFC) = comparisons
  }
  rownames(LogFC) = stats$Gene
  write.csv(LogFC, paste0(anal_address, StudyName, "_limma_LogFC.csv"))
  
  P_values <- fit2[["p.value"]]
  rownames(P_values) <- fData(eset)$symbol
  fdata_eset <- fData(eset)
  P_values <- cbind(P_values, external_gene_name = as.character(fData(eset)$ensgene), entrezgene = as.character(fData(eset)$entrez))
  P_values <- as.data.frame(P_values)
  write.csv(P_values, paste0(anal_address, StudyName, "_limma_Pvalues.csv"))
  #Get significant genes
  Sig_genes <- filter(stats, adj.P.Val < 0.05)
  write.csv(Sig_genes, paste0(anal_address, StudyName, "_limma_Sig_genes.csv"))
  #View data
  Histname = paste0(anal_address, StudyName, "_limma_Pvalue_hist.png")
  png(Histname)
  hist(stats[, "P.Value"])
  dev.off()
  volcanoName = paste0(anal_address, StudyName, "_limma_volcano.png")
  png(volcanoName)
  volcanoplot(fit2, highlight = 5, names = fit2$genes[,"symbol"])
  dev.off()
  
  # # #GSEA analysis
  setwd(anal_address)
  comp_list_table = LogFC
  comp_vec = colnames(comp_list_table)
  print(comp_vec)
  for(i in comp_vec){
    res_df = dplyr::select(stats, i, symbol)
    genelist = pull(res_df, i)
    names(genelist) = pull(res_df, symbol)
    genelist = sort(genelist, decreasing = TRUE)
    print(head(genelist))
    dir.create(paste0(i), showWarnings = TRUE)
    comp = as.character(i)
    gsea = GSEA(geneList = genelist,
                exponent = 1,
                minGSSize = 5,
                maxGSSize = 500,
                TERM2GENE = term2gene,
                pAdjustMethod = "none",
                pvalueCutoff = 0.05,
                by = "fgsea")
    gsea_df = as.data.frame(gsea)
    write.csv(gsea_df, paste0(i, "/", i, "_sig_gsea.csv"))
    gsea2 = GSEA(geneList = genelist,
                 exponent = 1,
                 minGSSize = 5,
                 maxGSSize = 500,
                 pvalueCutoff = 1,
                 pAdjustMethod = "none",
                 TERM2GENE = term2gene,
                 by = "fgsea")
    gsea_all = as.data.frame(gsea2)
    write.csv(gsea_all, paste0(i, "/", i, "_all_gsea.csv"))
    if(nrow(gsea_df) > 0){
      # Dotplot
      g1 = dotplot(gsea, color = "pvalue",  showCategory=nrow(gsea), split=".sign") + facet_grid(.~.sign)
      ggsave(paste0(i, "/", i, "_dotplot.png"), g1)
      # Ridge plot
      g3 = ridgeplot(gsea, fill = "pvalue") + labs(x = "enrichment distribution")
      ggsave(paste0(i, "/", i, "_ridgeplot.png"), g3)
    }
  }
}

####MegaP####
##LP I have created an MegaP value which combines the direction of the LogFC and the P value 
#Then get the mega p value
#To transform the p values associated with negative log fold change:
Transform_p <- function(P_value, LogFC, zero) {
  for(i in 1:nrow(P_value)) {
    for(j in 1:ncol(P_value)) {
      if(LogFC[i,j] < zero) {
        P_value[i,j] = -1/P_value[i,j]
      } else if(LogFC[i,j] > zero) {
        P_value[i,j] = 1/P_value[i,j]
      }
    }
  }
  return(P_value)
}

GetMegaP <- function(Pathway_Pvals, Pathway_logFC) {
  MegP_pathway <- Transform_p(Pathway_Pvals, Pathway_logFC, 0)
  #MegP_pathway <- Limit_p(MegP_pathway, 99999)
  return(MegP_pathway)
}

####create_comparison_data####
Create_Comparison_data = function(){
  P_values = read.csv(paste0(anal_address, StudyName, "_limma_Pvalues.csv"), row.names = 1)
  P_names = colnames(P_values)
  genes = rownames(P_values)
  col1 = P_names[1]
  P_values = as.data.frame(P_values[, 1:(length(P_values) -2)])
  P_cols = as.integer(ncol(P_values))
  if(P_cols == 1){
    colnames(P_values) = col1
    rownames(P_values) = genes
  }
  grch38 = grch38[!duplicated(grch38$symbol),]
  Unincluded_genes = grch38[which(!(grch38$symbol %in% rownames(P_values))),]
  Uninc_no = nrow(Unincluded_genes)
  NA_list = rep(NA, Uninc_no)
  NA_list = as.data.frame(NA_list)
  NA_table =  as.data.frame(NA_list[, rep(1, P_cols)])
  rownames(NA_table) = Unincluded_genes$symbol
  colnames(NA_table) = colnames(P_values)
  
  Total_Ps = rbind(P_values, NA_table)
  Total_Ps = cbind(Total_Ps, symbol = rownames(Total_Ps))
  Total_Ps = arrange(Total_Ps, symbol)
  rownames(Total_Ps) = Total_Ps$symbol
  if(P_cols > 1){
    Total_Ps = Total_Ps[, 1:(length(Total_Ps) -1)]
    TransPs = t(Total_Ps)
  } else if(P_cols == 1){
    Total_P_genes = rownames(Total_Ps)
    Total_Ps = as.data.frame(Total_Ps[, 1:(length(Total_Ps) -1)])
    rownames(Total_Ps) = Total_P_genes
    colnames(Total_Ps) = col1
    TransPs = t(Total_Ps)
  }
  rownames(TransPs)
  head(colnames(TransPs))
  
  Comparison_info = read.csv(paste0(anal_address, StudyName, "_comparison_info.csv"), row.names = 1)
  limma_line = rep("Limma", nrow(Comparison_info))
  Comparison_info = cbind(Comparison_info, Analysis = limma_line)
  Comparison_Pvalues = cbind(Comparison_info, TransPs)
  write.csv(Comparison_Pvalues, paste0(anal_address, StudyName, "_Comparison_limma_Pvalues.csv"))
  
  stats = read.csv(paste0(anal_address, StudyName, "_limma_stats.csv"), row.names = 1)
  rownames(stats) = stats$symbol
  if(P_cols > 1){
    stats = stats[, 10:(length(stats) -4)]
  } else if(P_cols == 1){
    gene_names = stats$symbol
    stats = as.data.frame(stats[, 10:(length(stats) -5)])
    colnames(stats) = col1
    rownames(stats) = gene_names
  }
  Total_log = rbind(stats, NA_table)
  Total_log = cbind(Total_log, symbol = rownames(Total_log))
  Total_log = arrange(Total_log, symbol)
  rownames(Total_log) = Total_log$symbol
  if(P_cols > 1){
    Total_log = Total_log[, 1:(length(Total_log) -1)]
    Translog = t(Total_log)
  } else if(P_cols == 1){
    Total_log_genes = rownames(Total_log)
    Total_log = as.data.frame(Total_log[, 1:(length(Total_log) -1)])
    rownames(Total_log) = Total_log_genes
    colnames(Total_log) = col1
    Translog = t(Total_log)
  }
  rownames(Translog)
  head(colnames(Translog))
  
  Comparison_log = cbind(Comparison_info, Translog)
  write.csv(Comparison_log, paste0(anal_address, StudyName, "_Comparison_limma_log.csv"))
  
  MegaP <- GetMegaP(P_values, stats)
  Total_MegP = rbind(MegaP, NA_table)
  Total_MegP = cbind(Total_MegP, symbol = rownames(Total_MegP))
  Total_MegP = arrange(Total_MegP, symbol)
  rownames(Total_MegP) = Total_MegP$symbol
  if(P_cols > 1){
    Total_MegP = Total_MegP[, 1:(length(Total_MegP) -1)]
    TransMegP = t(Total_MegP)
  } else if(P_cols == 1){
    Total_log_genes = rownames(Total_log)
    Total_MegP = as.data.frame(Total_MegP[, 1:(length(Total_MegP) -1)])
    rownames(Total_MegP) = Total_log_genes
    colnames(Total_MegP) = col1
    TransMegP = t(Total_MegP)
  }
  Comparison_MegP = cbind(Comparison_info, TransMegP)
  write.csv(Comparison_MegP, paste0(anal_address, StudyName, "_Comparison_limma_MegP.csv"))
}

####once all functions created#####
for(i in study_list){
  StudyName <- i
  print(StudyName)
  anal_address = paste0("/Volumes/Bex5TB/Sen23/AwiaitngResponseBeforeESET/",
                        StudyName, "/analysis/")
  eset = Create_and_save_eset()
  Design = Get_Design(eset)
  GetPandlogfromLimma(eset)
  Create_Comparison_data()
}

#For when script crashes so don't need to run EVERYTHING again
undesired = c("McHu")
study_list = study_list[!grepl(undesired, study_list)]

#for when running script for specific study because came back to it
desired ="Papa"
study_list = study_list[grep(desired, study_list)]
print(study_list)

desired <- c("Skea", "Hasegawa")
study_list <- study_list[grep(paste(desired, collapse = "|"), study_list, ignore.case = TRUE)]
print(study_list)

#identidying the name of the missing quant file
exists = file.exists(myfiles)
false = myfiles[!exists]
print(false)

####Guerrero eset ####
#because Guerrero has technical repeats, the quants are grouped by GSM number rather than SRR number
#therefore need to change 'run' to 'sample_acc' which is done below
study_list = list.dirs("/Volumes/Bex5TB/Sen23/Guerrero", full.names = FALSE, recursive = FALSE)
print(study_list)

Create_and_save_eset = function(){
  # wd_address = paste0("/media/njw262/MASSIVE/Transcriptomic_Analysis/RNAseq_Analysis4/", StudyName, "/Salmon/quant/")
  # setwd(wd_address)
  #mydir = paste0("/Volumes/Bex5TB/Sen23/test/", StudyName, "/salmon/")
  address = paste0("/Volumes/Bex5TB/Sen23/Guerrero/", StudyName, "/", StudyName, "_data_info.csv")
  #address = paste0("/Volumes/Bex5TB/Sen23/test/", StudyName, "/reference_folder/", StudyName, "_data_info.csv")
  p <- read.csv(address)
  pop <- rep("TSI", nrow(p))
  centre <- rep("UNIGE", nrow(p))
  p3 <- as.data.frame(cbind(p, pop, centre))
  #p3 = Add_replicate_row(p2)
  p3 <- mutate(p3, Fullname = paste(condition, replicate, sep = "_"))
  run <- paste0(p3$sample_acc, "_quant")
  condition <- p3$condition
  rownames(p3) <- p3$sample_acc
  mydir = paste0("/Volumes/Bex5TB/Sen23/Guerrero/", StudyName, "/quants")
  #Locate the quant.sf files and create an array of expression for expression eset:
  myfiles <- file.path(mydir, run, "quant.sf")
  names(myfiles) <- p3$sample_acc
  mygene = read.csv("/Volumes/Bex5TB/SenAnalysisFiles/tx2gene_ensembl_v87.csv")
  #mygene <- read_tsv("/media/njw262/DATA/Ageing_review/gencode_v32.gene.map.tsv", col_names = TRUE)
  # anal_address = paste0("/media/njw262/MASSIVE/Transcriptomic_Analysis/RNAseq_Analysis3/",
  #                       StudyName, "/analysis/")
  #NB lengthScaledTPM should only be needed if comparing transcripts across genes
  print(file.exists(myfiles))
  mytxi <- tximport(myfiles, type="salmon", tx2gene=mygene, ignoreTxVersion = TRUE, countsFromAbundance = "scaledTPM")
  #Use this one for DESeq (without TPM normalisation)
  mytxi2 <- tximport(myfiles, type="salmon", tx2gene=mygene, ignoreTxVersion = TRUE)
  file.exists(myfiles)
  
  #SAVE THE DESEQ DATA
  Save_data_for_DESeq(mytxi2, p3)
  
  #Create Eset
  dge <- DGEList(mytxi$counts)
  # filtering
  keep <- filterByExpr(dge)
  dge <- dge[keep, ]
  dge <- calcNormFactors(dge)
  #turn into alterable matrix and prepare col and rownames
  k <- as.data.frame(dge$counts)
  res = as.data.frame(rownames(k))
  colnames(res) = "ensgene"
  res_annot = left_join(res, grch38, by = "ensgene")
  res_annot2 <- res_annot[!duplicated(res_annot$ensgene),]
  res_annot2 <- res_annot2[!duplicated(res_annot2$symbol),]
  res_annot2 = res_annot2[which(res_annot2$ensgene != "NA" & res_annot2$symbol != "NA"),]
  rownames(res_annot2) = res_annot2$ensgene
  k <- k[which(rownames(k) %in% rownames(res_annot2) ==T),]
  raw_k_name = paste0(anal_address, StudyName, "_k_raw.csv")
  write.csv(k, raw_k_name)
  #Then normalise the data with cpm or voom depending on the variance
  totals = c(colSums(k))
  Red_line = max(totals)/ min(totals)
  print(Red_line)
  if(Red_line <= 3){
    dge = cpm(k, log=TRUE, prior.count = 3)
    cpm_name = paste0(anal_address, StudyName, "_cpm_density_plot.png")
    png(cpm_name)
    plotDensities(dge, legend = FALSE)
    dev.off()
  } else{
    vdes = as.matrix(p3$replicate)
    dge = voom(dge, vdes, plot = T)
    voom_name = paste0(anal_address, StudyName, "_voom_density_plot.png")
    png(voom_name)
    plotDensities(dge, legend = FALSE)
    dev.off()
  }
  #Remake k with normalised values
  k <- as.data.frame(dge)
  k <- k[which(rownames(k) %in% rownames(res_annot2) ==T),]
  p4 = p3
  rownames(p4) = p4$Fullname
  krun = colnames(k)
  knames = rownames(p4)
  p4run = as.character(p4$sample_acc)
  #To test you can randomise the vector
  #p5run = sample(p4run)
  true_vector = c()
  #See if all the run names match
  for(i in seq_along(krun)){
    if(krun[i] == p4run[i]){
      true_vector = append(true_vector, TRUE)
    } else {
      true_vector = append(true_vector, FALSE)
    }
  }
  #This tells you if all the values in true_vector are true
  m = all(true_vector)
  if(m == TRUE){
    colnames(k) = knames
  } else {
    print("You have got your run names mixed up!!!!!!")
  }
  
  grch2 <- grch38[!duplicated(grch38$ensgene),]
  grch2 <- grch2[!duplicated(grch2$symbol),]
  grch2 = grch2[which(grch2$ensgene != "NA" & grch2$symbol != "NA"),]
  
  k2 = cbind(k, symbol = res_annot2$symbol)
  rownames(k2) = k2$symbol
  k2 = k2[, 1:length(k)]
  other_genes = grch2[which(!(grch2$symbol %in% rownames(k2))),]
  NA_vec = c(rep(NA, nrow(other_genes)))
  NA_vec = as.data.frame(NA_vec)
  NA_matrix = NA_vec[rep(names(NA_vec), ncol(k2))]
  rownames(NA_matrix) = other_genes$symbol
  colnames(NA_matrix) = colnames(k2)
  total_genes = rbind(k2, NA_matrix)
  total_genes2 = arrange(total_genes, rownames(total_genes))
  total_genes3 = t(total_genes2)
  p5 = mutate(p4, Fullername = paste(study, condition, replicate, sep = "_"))
  combined_data = cbind(p5, total_genes3)
  combined_data2 = combined_data %>% pivot_longer(`5_8S_rRNA`:`ZZZ3`, names_to = "Gene", values_to = "norm_express")
  combined_data_name = paste0(anal_address, StudyName, "_normalised_expression.csv")
  write.csv(combined_data2, combined_data_name)
  Filtered_genes = combined_data2[which(combined_data2$Gene %in% rownames(k2)) ,]
  filtered_data_name = paste0(anal_address, StudyName, "_normalised_expression_NAfiltered.csv")
  write.csv(Filtered_genes, filtered_data_name)
  #Create an expression set
  eset <- ExpressionSet(assayData = as.matrix(k),
                        phenoData = AnnotatedDataFrame(p4),
                        featureData = AnnotatedDataFrame(res_annot2))
  SaveEset(k, p4, res_annot2, StudyName)
  
  MDSname = paste0(anal_address, StudyName, "_limma_MDS.png")
  png(MDSname)
  plotMDS(eset, labels = pData(eset)[, "label"], gene.selection = "common", cex = 0.6)
  dev.off()
  return(eset)
}
#also remember to re-run other functions to the right wd


#after running this function RS has moved Guerrero back to NewAnalysis folder to be with the rest of the data

####removing irregular curves####
#06dec23 Hasegawa and Skea both have one irregular curve each
#once removed, will need to rerun analysis
data = exprs(eset)
for (i in colnames(data)){
  x = data[,i]
  png(paste0(i, '.png'))
  plotDensities(x)
  dev.off()
}
#To more easily find which row nums are required
rownumbs = nrow(p4)
k2 = cbind(c(1:rownumbs), p4)
#georgilis_1 and marthadan_2 done by JW Skea and Hasagawa done by RS
#RS: to find my outlier sample i looked at PCA and MDS plot and then looked to see if obvious outlier in 'data' df
#For Skea
m = k[,c(1:7, 9:24)]
p3 = p4[c(1:7, 9:24),]
#For georgilis_1
m <- k[, c(1:124, 126, 128:256, 258:273, 275:276)]
p3 <- p4[c(1:124, 126, 128:256, 258:273, 275:276),]
#For marthadan_2
m <- k[, c(1:45, 47:60)]
p3 <- p4[c(1:45, 47:60),]
#Then:
eset <- ExpressionSet(assayData = as.matrix(m),
                      phenoData = AnnotatedDataFrame(p3),
                      featureData = AnnotatedDataFrame(res_annot2))
plotDensities(eset, legend = FALSE)
Dens_name = paste0(anal_address, StudyName, "_density_plot_redone.png")
png(Dens_name)
plotDensities(m, legend = FALSE)
dev.off()


grch2 <- grch38[!duplicated(grch38$ensgene),]
grch2 <- grch2[!duplicated(grch2$symbol),]
grch2 = grch2[which(grch2$ensgene != "NA" & grch2$symbol != "NA"),]

k2 = cbind(m, symbol = res_annot2$symbol)
rownames(k2) = k2$symbol
k2 = k2[, 1:length(m)]
other_genes = grch2[which(!(grch2$symbol %in% rownames(k2))),]
NA_vec = c(rep(NA, nrow(other_genes)))
NA_vec = as.data.frame(NA_vec)
NA_matrix = NA_vec[rep(names(NA_vec), ncol(k2))]
rownames(NA_matrix) = other_genes$symbol
colnames(NA_matrix) = colnames(k2)
total_genes = rbind(k2, NA_matrix)
total_genes2 = arrange(total_genes, rownames(total_genes))
total_genes3 = t(total_genes2)
p5 = mutate(p3, Fullername = paste(study, condition, replicate, sep = "_"))
combined_data = cbind(p5, total_genes3)
combined_data2 = combined_data %>% pivot_longer(`5_8S_rRNA`:`ZZZ3`, names_to = "Gene", values_to = "norm_express")
combined_data_name = paste0(anal_address, StudyName, "_normalised_expression.csv")
write.csv(combined_data2, combined_data_name)
Filtered_genes = combined_data2[which(combined_data2$Gene %in% rownames(k2)) ,]
filtered_data_name = paste0(anal_address, StudyName, "_normalised_expression_NAfiltered.csv")
write.csv(Filtered_genes, filtered_data_name)
write.csv(m, "k.csv")
write.csv(p3, "p4.csv")

