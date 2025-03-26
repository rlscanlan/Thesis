options(warn=-1)
library(clusterProfiler)
library(wordcloud)
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
BiocManager::install("clusterProfiler")
library(pathview)
# reading in input from deseq2

#setwd("/Users/rebekahscanlon/Desktop/Rebekah/from_ind_folders")
setwd("/Users/rebekahscanlon/Desktop/Copenhagen/APTX/Rebekah_PC/from_ind_folders/2024")
#df = read.csv("Adult_vs_Old.sig.genes.tsv", header=TRUE)
#df = read.table("genes_in_ballgown_and_salmon_for_APTX1_NS_vs_APTX0_NS.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
#df = read.table("genes_in_ballgown_and_salmon_for_APTX1_NS_vs_APTX0_NS.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
#df = read.csv(file="APTX0_NS_vs_APTX0_IS.sig.genes.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
#sigs <- list.files(pattern=".sig.genes.tsv")
#sigs <- list.files(pattern="^module_table_")
#sigs <- list.files(pattern="_node_table.csv$")
#names(sigs) <- gsub("_node_table.csv", "", sigs)
#sigs <- sigs[order(names(sigs))]
express <- list.files(path="/Users/rebekahscanlon/Desktop/Copenhagen/APTX/Rebekah_PC/from_ind_folders/", pattern="^ALL_significantly_differentially_ent_gene")
names(express) <- gsub(".sig.genes.tsv", "", express)
express <- express[order(names(express))]
#if(names(sigs) == names(expr)){
#s <- sigs[which(names(sigs) == names(expr))]
#print("TRUE")
#}
#all_hubs <- matrix(, nrow=length(sigs), ncol=21)
#rownames(all_hubs) <- names(sigs)
#all_hubs <- as.data.frame(all_hubs)
#all_subs <- matrix(, nrow=length(sigs), ncol=21)
#rownames(all_hubs) <- names(sigs)
#all_subs <- as.data.frame(all_subs)
#colnames(all_hubs) <- colnames(sigs[s])
#edf = read.csv(path="/media/louise/Scratch/Mansour_APTX/results", file=paste0(sigs[s]), sep="\t", stringsAsFactors=FALSE)
#for (s in 1:length(sigs)){
#e <- express[which(names(sigs[s]) == names(express))]
#df = read.csv(file=paste0(sigs[s]), header=TRUE, sep="\t", stringsAsFactors=FALSE)
#df = read.csv(file=paste0(sigs[s]), header=TRUE, sep=",", stringsAsFactors=FALSE)
#df <- df[order(df$Betweenness, decreasing=TRUE),]
#df <- df[order(df$Degree, decreasing=TRUE),]
#plot(df$Betweennes)
#plot(df$Degree)
#if(names(sigs[s]) == names(expr)){
for (e in 1:length(express)){
edf = read.csv(file=paste0("/Users/rebekahscanlon/Desktop/Copenhagen/APTX/Rebekah_PC/from_ind_folders/", express[e]), header=TRUE, sep="\t", stringsAsFactors=FALSE)
#mrgd <- merge(edf, df, by.x="symbol", by.y="Label")
df <- edf
df$log2FoldChange <- log2(df$fc)
#df$Comparison <- names(sigs[s])
nm <- gsub("ALL_significantly_differentially_ent_transcript_expression_in_bg_group", "", names(express[e]))
nm <- gsub(".txt", "", nm)
nm <- paste0("transcripts", nm)
#df <- df[which(df$entrez != "NA"),]
if(nrow(df) > 0){
#df$description <- mget(as.character(df$entrez), org.Hs.egGENENAME, ifnotfound=NA)
##df$ensg <- mget(as.character(df$entrez), org.Hs.egENSEMBL, ifnotfound=NA)
df <- na.omit(df)
df$entrez <- mget(as.character(df$gene_name), org.Hs.egENSEMBL, ifnotfound=NA)
#df$entrez <- lapply(df$entrez, '[[', 1)
df$description2 <- mget(as.character(df$entrez), org.Hs.egGENENAME, ifnotfound=NA)
#df$description <- mget(as.character(df$entrez), org.Hs.egGENENAME, ifnotfound=NA)
#df$ensg <- mget(as.character(df$entrez), org.Hs.egENSEMBL, ifnotfound=NA)
#df$ensgene <- mget(as.character(df$entrez), org.Hs.egENSEMBL, ifnotfound=NA)
#df$ensgene <- lapply(df$ensgene, '[[', 1)
if(nrow(df) > 0){
#df$description <- mget(as.character(df$entrez), org.Hs.egGENENAME, ifnotfound=NA)
df2 <- df
	df2$name <- substr(df$gene_name, start=1, stop=3)
	dffreq <- as.data.frame(table(unlist(strsplit(tolower(df2$name), " "))))
	#pdf(file=paste0("word_cloud_name", nm, ".pdf"))
	#print(wordcloud(words = dffreq$Var1, freq = dffreq$Freq, scale=(c(4, .5)), colors=brewer.pal(8, "Dark2"), max.words = 25))
	#dev.off()
	#pdf(file=paste0("barplot_name", nm, ".pdf"))
	#print(barplot(dffreq$Freq, colors=brewer.pal(8, "Dark2"))
	#dev.off()
	df2 <- df
	#df2$description <- substr(df$gene, start=1, stop=3)
	dffreq <- as.data.frame(table(unlist(strsplit(tolower(df2$description2), " "))))
	#pdf(file=paste0("word_cloud_description", nm, ".pdf"))
	#print(wordcloud(words = dffreq$Var1, freq = dffreq$Freq, scale=(c(4, .5)), colors=brewer.pal(8, "Dark2"), max.words = 25))
	#dev.off()


#M <- mean(df$Betweenness)
#D <- mean(df$Degree)
#hub <- df[which(max(df$Betweenness),]
#x <- max(df$Betweenness)
#x <- mean(df$Betweenness)*6
#hub<- as.data.frame(df[which(df$Betweenness == x),])
#hub<- as.data.frame(df[which(df$Betweenness > x),])
#colnames(all_hubs) <- colnames(df)
#all_hubs <- rbind(hub, all_hubs)
#suby <- df[which(df$symbol %in% ints),]
#su[, s] <- suby$log2FoldChange[which(rownames(su) == rownames(suby)),]
#rownames(suby) <- suby$symbol
#colnames(all_subs) <- colnames(suby)
#all_subs <- rbind(suby, all_subs)
#int_sub <- merge(su, suby, by = "row.names", all.x=TRUE, all.y=TRUE)
#int_sub <- merge(su, suby, by = "row.names", all.x=TRUE, all.y=TRUE)
#if(colnames(su) == suby$Comparison & rownames(suby) == rownames(su)){

#for(rw in 1:nrow(su)){
#su[rw, s] <- suby$log2FoldChange[which(rownames(su) %in% rownames(suby))]
#su[, s] <- suby$log2FoldChange[which(suby$symbol == rownames(su),]
#merge(suby$log2FoldChange, su, by="row.names", all=TRUE)

#if(suby$symbol %in% rownames(su)){
#su[, s] <- suby$log2FoldChange
#suby$log2FoldChange[which(suby$symbol %in% rownames(su))]
#merge(suby$log2FoldChange, su, by="row.names", all=FALSE)
#su <- merge(all_subs, su, by.x="symbol", by.y="rownames")
#su[, s]
#heat_dat <- merge(suby$log2FoldChange,
#original_gene_list <- df$Degree
#original_gene_list <- df$Betweenness
original_gene_list <- df$log2FoldChange
# name the vector
#names(original_gene_list) <- df$X
#names(original_gene_list) <- df$t_name
names(original_gene_list) <- df$id
#names(original_gene_list) <- df$sym

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# Exctract significant results (padj < 0.05)
#sig_genes_df = subset(df, Degree > 10)
#sig_genes_df = subset(df, Betweenness > M)
sig_genes_df = subset(df, log2FoldChange > 2 | log2FoldChange < -2)
#sig_genes_df = subset(df, qval < 0.05)
#sig_genes_df = subset(sig_genes_df, Betweenness > M)
#sig_genes_df = subset(sig_genes_df, Degree > D)
# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$log2FoldChange
#genes <- sig_genes_df$log2FoldChanges
# Name the vector
#names(genes) <- sig_genes_df$t_name
names(genes) <- sig_genes_df$id
# omit NA values
genes <- na.omit(genes)

# filter on min log2fold change (log2FoldChange > 2)
#genes <- names(genes)[abs(genes) > 1]
#genes <- genes[which(genes != "Inf" & genes != "-Inf")]


go_enrich <- enrichGO(gene = names(genes),
                      universe = names(gene_list),
                      OrgDb = organism, 
                      keyType = 'ENSEMBL',
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.05)

library(enrichplot)
if(!is.null(go_enrich)){
if(nrow(summary(go_enrich)) > 0){
	upsetplot(go_enrich)
	wcdf<-read.table(text=go_enrich$GeneRatio, sep = "/")[1]
	wcdf$term<-go_enrich[,2]
	#pdf(file=paste0("word_cloud_go_enrich", nm, ".pdf"))
	#print(wordcloud(words = wcdf$term, freq = wcdf$V1, scale=(c(2, .5)), colors=brewer.pal(8, "Dark2"), max.words = 25))
	#dev.off()
	pdf(file=paste0("barplot_go_enrich", nm, ".pdf"))
	print(barplot(go_enrich, 
        	drop = TRUE, 
        	showCategory = 10, 
        	title = "GO Biological Pathways",
        	font.size = 8))
	dev.off()
	pdf(file=paste0("dot_plot_go_enrich", nm, ".pdf"))
	print(dotplot(go_enrich))
	dev.off()
	#pdf(file=paste0("ema_go_enrich", nm, ".pdf"))
	#print(emapplot(go_enrich))
	#dev.off()

###Enriched GO induced graph:
	pdf(file=paste0("goplot_go_enrich", nm, ".pdf"))
	print(goplot(go_enrich, showCategory = 10))
	dev.off()

# categorySize can be either 'pvalue' or 'geneNum'
	pdf(file=paste0("cne_go_enrich", nm, ".pdf"))
#cnetplot(go_enrich, categorySize="pvalue", foldChange=gene_list)
	print(cnetplot(go_enrich, categorySize="geneNum", foldChange = gene_list))
	dev.off()
	#} else {
	#break;
}
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db") # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
dedup_ids = ids[!duplicated(ids[c("ENTREZID")]),]
# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
#df2 = df[df$ensgene %in% dedup_ids$ENSEMBL,]
#df2 = df2[df2$entrez %in% dedup_ids$ENTREZID,]
df2 <- df 
# Create a new column in df2 with the corresponding ENTREZ IDs
#df2 <- df2[which(df2$entrez != "NA"),]
#df2$Y = dedup_ids$ENTREZID
df2$Y = df2$entrez

# Create a vector of the gene unuiverse
#kegg_gene_list <- df2$Betweenness
kegg_gene_list <- df2$log2FoldChange
# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
#kegg_gene_list <- kegg_gene_list[which(kegg_gene_list != "Inf" & kegg_gene_list != "-Inf")]
# Exctract significant results from df2
#M <- mean(df$Betweenness)
#[1] 2285.04

#kegg_sig_genes_df = subset(df2, Betweenness > M)
#kegg_sig_genes_df = subset(kegg_sig_genes_df, Degree > D)
kegg_sig_genes_df = subset(df2, qval < 0.05)
#kegg_sig_genes_df = subset(df2, log2FoldChange > 1 | log2FoldChange < -1)
# From significant results, we want to filter on log2fold change
#kegg_genes <- kegg_sig_genes_df$Betweenness
kegg_genes <- kegg_sig_genes_df$log2FoldChange
# Name the vector with the CONVERTED ID!
names(kegg_genes) <- kegg_sig_genes_df$entrez

# omit NA values
kegg_genes <- na.omit(kegg_genes)

# filter on log2fold change (PARAMETER)
#kegg_genes <- names(kegg_genes)[abs(kegg_genes) > 1]
#kegg_genes <- kegg_genes[which(kegg_genes != "Inf" & kegg_genes != "-Inf")]
kegg_organism = "hsa"
kk <- enrichKEGG(gene=names(kegg_genes), universe=names(kegg_gene_list),organism=kegg_organism, pvalueCutoff = 0.05, keyType = "ncbi-geneid")
if(!is.null(kk)){
if(nrow(summary(kk)) > 0){

#sel = which(names(gs.annots[["kegg"]]@original) %in% problematic_pathways) 

pdf(file=paste0("barplot_pathways", nm, ".pdf"))
print(barplot(kk, 
        showCategory = 10, 
        title = "Enriched Pathways",
        font.size = 8))
dev.off()
pdf(file=paste0("dotplot_pathways", nm, ".pdf"))
print(dotplot(kk, 
        showCategory = 10, 
        title = "Enriched Pathways",
        font.size = 8))
dev.off()

library(pathview)
`%!in%` = Negate(`%in%`)
problematic_pathways <- c("Kaposi sarcoma-associated herpesvirus infection", "MicroRNAs in cancer") #c("hsa05167", "hsa05206")
kk <- kk[which(kk$Description %!in% problematic_pathways),]
# Produce the native KEGG plot (PNG)
options(bitmapType='cairo')
#png(file=paste0("pathway_cGAS_", nm, ".png"))
hsa <- pathview(gene.data=kegg_gene_list, pathway.id="hsa04623", species = kegg_organism, gene.idtype="entrez", out.suffix=paste0("cGAS_pathway_", nm))
knitr::include_graphics("hsa04623.pathview.png", error = FALSE)
#hsa <- pathview(gene.data=kegg_gene_list, pathway.id=kk$ID, species = kegg_organism, gene.idtype="entrez", out.suffix=paste0("pathway_", nm))
#knitr::include_graphics(paste0(kk$ID, ".pathview.png"), error=FALSE)
#knitr::include_graphics(paste0("pathway_interactions_", nm, "hsa04623.pathview.png"))
#dev.off()
# Produce a different plot (PDF) (not displayed here)
#options(bitmapType='cairo')
#pdf(file=paste0("pathway_interactions_cGAS_", nm, ".pdf"))
#hsa <- pathview(gene.data=kegg_gene_list, pathway.id=kk$ID, species = kegg_organism, gene.idtype="entrez", kegg.native = F, out.suffix=paste0("interactions_", nm))
#knitr::include_graphics(paste0(kk$ID, ".pathview.png"), error=FALSE)
#knitr::include_graphics(file=paste0("pathway_interactions_", nm, "hsa04623.pathview.png"), "hsa04623.pathview.png")
#dev.off()
}
}
}
}
}
}
}
print("F")

library(clusterProfiler)
#all_hubs_df <- as.data.frame(do.call(rbind, all_hubs))
all_hubs_df <- all_hubs[which(all_hubs$symbol != "<NA>"),]
all_hubs_df <- all_hubs_df[which(all_hubs_df$Degree > 10),]
#heatmap(as.matrix(all_hubs_df$log2FoldChange))
rownames(all_hubs_df) <- all_hubs_df$Comparison
all_subs_df <- all_subs[which(all_subs$symbol != "<NA>"),]
all_subs_df$code <- paste0(all_subs_df$symbol, "_", all_subs_df$Comparison)
st <- substr(all_subs_df$Comparison, 5, 8)
ed <- substr(all_subs_df$Comparison, 17, 20)


all_subs_df$shtnm <- paste0(st, "_v_", ed)

#rownames(all_subs_df) <- all_subs_df$Comparison
#all_hubs_df <- as.data.frame(do.call(rbind, all_hubs_df))
dfu <- apply(all_hubs_df,2,as.character)
#write.table(all_hubs_df, file="hubs_for_APTX_networks.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
write.table(dfu, file="hubs_for_APTX_networks.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
ints <- all_hubs_df$symbol
su <- matrix(, nrow=length(ints), ncol=length(sigs))
colnames(su) <- names(sigs)
rownames(su) <- ints

#merge(all_subs_df, by = "symbol", all.x=TRUE, all.y=TRUE)
#merge(su, suby, by = "row.names", all.x=TRUE, all.y=TRUE)
#tall_subs <- as.matrix(do.call((all_subs)
mycols <- c("darkblue","red", "deeppink", "green", "purple")


colourCount = length(unique(all_subs_df$symbol))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
colrs<- getPalette(colourCount)
names(colrs) <- unique(all_subs_df$symbol)
all_subs_df$colour <- colrs[all_subs_df$symbol]  


pdf(file="barplot_of_log2_FC_for_identified_hubs_in_APTX_cell_lines.pdf")
barplot(all_subs_df$log2FoldChange, col=all_subs_df$colour, names.arg=all_subs_df$shtnm, cex.names=0.7, las=2, legend = as.factor(all_subs_df$symbol))
dev.off()

#barplot(j$dbh, names.arg = j$plot, ylab = "dbh" , 
#    col = mycols[j$year], ylim=c(0,70))

#heatmap(dfu$log2FoldChange)
#barplot(all_subs_df$log2FoldChange, col=c("darkblue","red", "deeppink", "green", "black"), names.arg=all_subs_df$Comparison)
#barplot(all_subs_df$log2FoldChange, names.arg=all_subs_df$Comparison)
#barplot(all_subs_df$log2FoldChange, col=c("darkblue","red", "green"), names.arg=all_subs_df$Comparison)
#barplot(all_subs_df$log2FoldChange, col=all_subs_df$Comparison, names.arg=all_subs_df$symbol, cex.names=0.5)
#barplot(all_subs_df$log2FoldChange, col=c("darkblue","red", "deeppink", "green", "black"), names.arg=all_subs_df$Comparison, cex.names=0.5, legend = all_subs_df$symbol)
pdf(file="barplot_of_log2_FC_for_identified_hubs_in_APTX_cell_lines.pdf")
barplot(all_subs_df$log2FoldChange, col=mycols[as.factor(all_subs_df$symbol)], names.arg=all_subs_df$shtnm, cex.names=0.7, las=2, legend = as.factor(all_subs_df$symbol))
dev.off()

pdf(file="barplot_of_log2_FC_for_identified_hubs_in_APTX_cell_lines.pdf")
barplot(all_subs_df$log2FoldChange, col=mycols[as.factor(all_subs_df$shtnm)], names.arg=all_subs_df$shtnm, cex.names=0.7, las=2, legend = as.factor(all_subs_df$shtnm))
dev.off()

barplot(all_subs_df$log2FoldChange, col=c("darkblue","red", "deeppink", "green", "black"), names.arg=all_subs_df$code, cex.names=0.5, las=2)
