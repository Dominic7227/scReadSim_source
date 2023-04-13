library(Matrix)
library(Rsubread)
library(pscl)
library(parallel)
library(MASS)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(reshape2)  # melt
library(grid)
library(gridExtra)
library(plyr)
library(RColorBrewer)
library(cowplot) # ggsave

samplename <- "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.UMI.countmatrix"
filename_pre <- sprintf("%s_VerifyCount", samplename)
sample.format <- "txt"
directory <- "/home/gayan/Projects/scATAC_Simulator/results/20221011_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT"

barcode.file <- sprintf("%s/synthetic_cell_barcode.txt",directory) 

out_directory <- "/home/gayan/Projects/scATAC_Simulator/results/20221011_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/verification_UMI_count.noAlevin"
dir.create(out_directory)


## Generate synthetic count matrix
# Read in count matrix
cat(sprintf("Reading scReadSim count matrix %s.%s...\n", samplename, sample.format))
scDesign2_matrix <- read.table(sprintf("%s/e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.UMI.countmatrix.scDesign2Simulated.txt", directory), sep="\t",header = TRUE)
matrix_scDesign2 <-scDesign2_matrix
# matrix_scDesign2 <- scDesign2_matrix[,2:ncol(scDesign2_matrix)]
cat("Reading scDesign2 count matrix...\n")
chr1_genes <- read.table("/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE/gene_region/gencode.vM10.annotation.gene_region.chr1.gtf", sep="\t",header = FALSE)
scDesign2_genes_merged <- read.table("/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE/gene_region/gencode.vM10.annotation.gene_region.chr1.merged.bed", sep="\t",header = FALSE)
scDesign2_genecode <- unlist(lapply(sapply(chr1_genes[match(scDesign2_genes_merged[,2], chr1_genes[,4]),9], strsplit, split=";", fixed=TRUE), function(x) {
  return(strsplit(x[1], split=" ",fixed=TRUE)[[1]][2])
}), use.names=FALSE)
scDesign2_genecode <- unlist(lapply(sapply(scDesign2_genecode, strsplit, split=".", fixed = TRUE), function(x) return(x[1])), use.names = FALSE)

## Construct gene geneid list
# Gene region could overlap (different strands)
# scDesign2 count matrix row refers to merged genes, hence one row links to multiple genes potentially 
chr1_geneid <-  unlist(lapply(sapply(chr1_genes[,9], strsplit, split=";", fixed=TRUE), function(x) {
  return(strsplit(x[1], split=" ",fixed=TRUE)[[1]][2])
}), use.names=FALSE)
chr1_geneid <- unlist(lapply(sapply(chr1_geneid, strsplit, split=".", fixed = TRUE), function(x) return(x[1])), use.names = FALSE)
chr1_genename <- unlist(lapply(sapply(chr1_genes[,9], strsplit, split=";", fixed=TRUE), function(x) {
  return(strsplit(x[4], split=" ",fixed=TRUE)[[1]][3])
}), use.names=FALSE)
chr1_gene_geneid_list <- cbind(chr1_geneid, chr1_genename)
# Create a gene to region index list. Same length to chr1 gene list
scReadSim_gene2region_indice <- rep(1:length(scDesign2_genecode), c(diff(match(scDesign2_genes_merged[,2], chr1_genes[,4])), 1))
# write.table(chr1_gene_geneid_list[,1], file="/home/gayan/Projects/scATAC_Simulator/results/20221011_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/verification_alevin/ref_data_chr1/chr1_geneID.txt",col.names=FALSE, row.names=FALSE,quote = FALSE)

## UMItools
umitool.out.directory <- sprintf("%s/verification_umitools", directory)
matrix_umitools <- read.csv(sprintf("%s/verification_umitools/counts.tsv", directory), sep="\t", header = TRUE)
rownames(matrix_umitools) <- matrix_umitools[,1]
matrix_umitools <- matrix_umitools[,-1]
# Filter
umitools_filter_id <- which(rowSums(matrix_umitools)>0)
matrix_umitools <- matrix_umitools[umitools_filter_id,]
umitools_genecode <- unlist(lapply(sapply(rownames(matrix_umitools), strsplit, split=".", fixed = TRUE), function(x) return(x[1])), use.names = FALSE)
rownames(matrix_umitools) <- umitools_genecode

## Alevin
alevin.out.directory <- sprintf("%s/verification_alevin", directory)
matrix_Alevin <- as.matrix(Matrix::readMM(sprintf("%s/alevin/quants_mat.mtx", alevin.out.directory))) # dim 4877 by 31460
alevin_barcode <- unlist(read.table(sprintf("%s/alevin/quants_mat_rows.txt", alevin.out.directory), sep="\t",header = FALSE), use.names=FALSE)
alevin_genes <- unlist(read.table(sprintf("%s/alevin/quants_mat_cols.txt", alevin.out.directory), sep="\t",header = FALSE), use.names=FALSE)
# Filter
Alevin_filter_id <- which(colSums(matrix_Alevin)>0)
matrix_Alevin <- t(matrix_Alevin[,Alevin_filter_id])
alevin_genes <- alevin_genes[Alevin_filter_id]
alevin_genes <- unlist(lapply(sapply(alevin_genes, strsplit, split=".", fixed = TRUE), function(x) return(x[1])), use.names = FALSE)
colnames(matrix_Alevin) <- alevin_barcode
rownames(matrix_Alevin) <- alevin_genes
## Whole genome annotation
whole_genome_genes <- read.table("/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE/gencode.vM10.primary_assembly.annotation.gene_only.gtf", sep="\t",header = FALSE)
# linux
# awk -F'\t' '$3=="gene"' /home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE/gencode.vM10.primary_assembly.annotation.gtf > /home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE/gencode.vM10.primary_assembly.annotation.gene_only.gtf
whole_genome_geneid <-  unlist(lapply(sapply(whole_genome_genes[,9], strsplit, split=";", fixed=TRUE), function(x) {
  return(strsplit(x[1], split=" ",fixed=TRUE)[[1]][2])
}), use.names=FALSE)
whole_genome_geneid <- unlist(lapply(sapply(whole_genome_geneid, strsplit, split=".", fixed = TRUE), function(x) return(x[1])), use.names = FALSE)
alevin_id_inWhole <- sapply(alevin_genes, function(x){
  which(whole_genome_geneid == x)
})
alevin_gene_chr <- whole_genome_genes[alevin_id_inWhole,1]
table(alevin_gene_chr)

## Cellranger
matrix_cellranger <- as.matrix(Matrix::readMM(sprintf("%s/verification_cellranger/cellranger_results/outs/filtered_feature_bc_matrix/matrix.mtx", directory)))
cellranger_barcode_tmp <- unlist(read.table(sprintf("%s/verification_cellranger/cellranger_results/outs/filtered_feature_bc_matrix/barcodes.tsv", directory), sep="\t",header = FALSE), use.names=FALSE)
cellranger_barcode <- unlist(lapply(strsplit(cellranger_barcode_tmp, "-"), function(x){x[[1]]}))
cellranger_genes <- read.table(sprintf("%s/verification_cellranger/cellranger_results/outs/filtered_feature_bc_matrix/features.tsv", directory), sep="\t",header = FALSE)
# gene_genecode_list <- cellranger_genes
# Filter
cellranger_filter_id <- rowSums(matrix_cellranger) > 0
matrix_cellranger <- matrix_cellranger[cellranger_filter_id,]
cellranger_genecode <- cellranger_genes[cellranger_filter_id,1]
cellranger_genecode <- unlist(lapply(sapply(cellranger_genecode, strsplit, split=".", fixed = TRUE), function(x) return(x[1])), use.names = FALSE)
colnames(matrix_cellranger) <- cellranger_barcode
rownames(matrix_cellranger) <- cellranger_genecode

########################### Matrix matching ###########################
## Matching columns
# Extract real clustering to match synthetic UMI matrix
clustering_synthetic <- read.table(sprintf("%s/synthetic_cell_barcode.txt.withSynthCluster", directory), sep="\t",header = FALSE)
clustering_synthetic_cluster <- clustering_synthetic[,2]
names_withdots <- clustering_synthetic_cluster[grepl(".", clustering_synthetic_cluster, fixed=TRUE)]
clustering_synthetic_cluster[grepl(".", clustering_synthetic_cluster, fixed=TRUE)] <- unlist(lapply(names_withdots, function(x){strsplit(toString(x), ".",  fixed=TRUE)[[1]][1]}))
# Extract synthetic clustering to match UMItools matrix
# simu_umi_clusters <- clustering_synthetic_cluster[match(colnames(matrix_umitools), clustering_synthetic[,1])]
# # Extract synthetic clustering to match Alevin matrix
# simu_alevin_clusters <- clustering_synthetic_cluster[match(colnames(matrix_Alevin), clustering_synthetic[,1])]
# # Extract synthetic clustering to match celranger matrix
# simu_cellranger_clusters <- clustering_synthetic_cluster[match(colnames(matrix_cellranger), clustering_synthetic[,1])]

## Matching rows
# Compare with chr1 gene list
length(chr1_gene_geneid_list[,2]) # 3428
length(scDesign2_genecode) # 2165
# UMItools
length(intersect(umitools_genecode, chr1_gene_geneid_list[,1])) # 2194/2259
# Alevin
length(intersect(alevin_genes, chr1_gene_geneid_list[,1])) # 965/2210
# Cellranger
length(intersect(cellranger_genecode, chr1_gene_geneid_list[,1])) # 2216/2423
# Combine the processed gene by cell count matrix to region by cell, corresponding to scDesign2's count matrix.
merge_genebycell_mat <- function(input_genebycell_mat, chr1_gene_geneid_list, scReadSim_gene2region_indice){
  # input_genebycell_mat <- matrix_umitools
  input_genecode <- rownames(input_genebycell_mat)
  num_intersect_gene <- length(intersect(input_genecode, chr1_gene_geneid_list[,1]))
  print(num_intersect_gene)
  # Match columns
  matched_mat_col <- input_genebycell_mat[,match(clustering_synthetic[,1], colnames(input_genebycell_mat))]
  matched_mat_col[is.na(matched_mat_col)] <- 0
  colnames(matched_mat_col) <- clustering_synthetic[,1]
  # Match rows
  matched_mat <- matched_mat_col[match(chr1_gene_geneid_list[,1], input_genecode), ]
  matched_mat[is.na(matched_mat)] <- 0
  rownames(matched_mat) <- chr1_gene_geneid_list[,1]
  # Gene 2 regions
  merged_mat <- matrix(0, nrow=nrow(scDesign2_matrix), ncol=ncol(scDesign2_matrix))
  for (i in 1:nrow(scDesign2_matrix)){
    gene_ind <- which(scReadSim_gene2region_indice == i)
    if (length(gene_ind) > 1){
      merged_mat[i,] <- colSums(matched_mat[gene_ind,] )
    } else {
      merged_mat[i,] <- unlist(matched_mat[gene_ind,], use.names = FALSE)
    }
  }
  colnames(merged_mat) <- clustering_synthetic[,1]
  rownames(merged_mat) <- rownames(matrix_scDesign2)
  return(list(matched_mat=matched_mat, merged_mat=merged_mat))
}

maplist_cellranger <- merge_genebycell_mat(matrix_cellranger, chr1_gene_geneid_list, scReadSim_gene2region_indice)
maplist_Alevin <- merge_genebycell_mat(matrix_Alevin, chr1_gene_geneid_list, scReadSim_gene2region_indice)
maplist_umitools <- merge_genebycell_mat(matrix_umitools, chr1_gene_geneid_list, scReadSim_gene2region_indice)
load("/home/gayan/Projects/scATAC_Simulator/results/20221011_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/verification_UMI_count/MatchedProcessedUMIMat_20221013.Rdata")
save(maplist_cellranger, maplist_Alevin, maplist_umitools, file=sprintf("%s/MatchedProcessedUMIMat_20221013.Rdata", out_directory))

load(file=sprintf("%s/MatchedProcessedUMIMat_20221013.Rdata", out_directory))
## Keep overlapped nonzero genes
nz_scDesign2_geneReg <- which(rowSums(matrix_scDesign2) > 0)
nz_cellranger_geneReg <- which(rowSums(maplist_cellranger$merged_mat) > 0)
# nz_Alevin_geneReg <- which(rowSums(maplist_Alevin$merged_mat) > 0)
nz_umitools_geneReg <- which(rowSums(maplist_umitools$merged_mat) > 0)
# Include Alevin
# intersect_gene_ID <- Reduce(intersect, list(nz_scDesign2_geneReg, nz_cellranger_geneReg, nz_umitools_geneReg, nz_Alevin_geneReg))
# Not Include Alevin
intersect_gene_ID <- Reduce(intersect, list(nz_scDesign2_geneReg, nz_cellranger_geneReg, nz_umitools_geneReg))
## Keep overlapped nonzero cells
nz_scDesign2_cell <- which(colSums(matrix_scDesign2) > 0)
nz_cellranger_cell <- which(colSums(maplist_cellranger$merged_mat) > 0)
# nz_Alevin_cell <- which(colSums(maplist_Alevin$merged_mat) > 0)
nz_umitools_cell <- which(colSums(maplist_umitools$merged_mat) > 0)
# Include Alevin
# intersect_cell_ID <- Reduce(intersect, list(nz_scDesign2_cell, nz_cellranger_cell, nz_Alevin_cell, nz_umitools_cell))
# Not Include Alevin
intersect_cell_ID <- Reduce(intersect, list(nz_scDesign2_cell, nz_cellranger_cell, nz_umitools_cell))

## processed Matrices
processed_scDesign2 <- matrix_scDesign2[intersect_gene_ID, intersect_cell_ID]
processed_cellranger <- maplist_cellranger$merged_mat[intersect_gene_ID, intersect_cell_ID]
# processed_Alevin <- maplist_Alevin$merged_mat[intersect_gene_ID, intersect_cell_ID]
processed_umitools <- maplist_umitools$merged_mat[intersect_gene_ID, intersect_cell_ID]
save(processed_scDesign2, processed_cellranger, processed_umitools, file=sprintf("%s/Filtered.noAlevin.MatchedProcessedUMIMat_20221015.Rdata", out_directory))
save(intersect_gene_ID, intersect_cell_ID, file=sprintf("%s/Filtered.noAlevin.MatchedProcessedUMIMat.GeneIDandCellID_20221015.Rdata", out_directory))

## TODO
## 1. Add two metrics: # identified chr1 gene and running time
load(sprintf("%s/Filtered.noAlevin.MatchedProcessedUMIMat_20221015.Rdata", out_directory))

########################### Quantity. Correlation per gene ###########################
filter_id_gene <- 1:nrow(processed_scDesign2)
## Pearson Cor.
cor_cellranger <- unlist(lapply(filter_id_gene, function(id){
  return(cor(unlist(processed_scDesign2[id,]), unlist(processed_cellranger[id,])))
}), use.names=FALSE) 
# cor_Alevin <- unlist(lapply(filter_id_gene, function(id){
#   return(cor(unlist(processed_scDesign2[id,]), processed_Alevin[id,]))
# }), use.names=FALSE) 
cor_umitools <- unlist(lapply(filter_id_gene, function(id){
  return(cor(unlist(processed_scDesign2[id,]), processed_umitools[id,]))
}), use.names=FALSE) 
## Kendall's tau
Kcor_cellranger <- unlist(lapply(filter_id_gene, function(id){
  return(cor(unlist(processed_scDesign2[id,]), unlist(processed_cellranger[id,]), method="kendall"))
}), use.names=FALSE) 
# Kcor_Alevin <- unlist(lapply(filter_id_gene, function(id){
#   return(cor(unlist(processed_scDesign2[id,]), processed_Alevin[id,], method="kendall"))
# }), use.names=FALSE)
Kcor_umitools <- unlist(lapply(filter_id_gene, function(id){
  return(cor(unlist(processed_scDesign2[id,]), processed_umitools[id,], method="kendall"))
}), use.names=FALSE) 
print(mean(cor_cellranger))
# print(mean(cor_Alevin))
print(mean(cor_umitools))
print(mean(Kcor_cellranger))
# print(mean(Kcor_Alevin))
print(mean(Kcor_umitools))
save(cor_cellranger, cor_umitools, Kcor_cellranger, Kcor_umitools, file=sprintf("%s/VerifyPipeline.noAlevin.Quantity.CorrelationPerGene_20221015.Rdata", out_directory))

########################### Quantity. Correlation per cell ###########################
filter_id_cell <- 1:ncol(processed_scDesign2)
## Pearson Cor.
cor_cellranger.percell <- unlist(lapply(filter_id_cell, function(id){
  return(cor(unlist(processed_scDesign2[,id]), processed_cellranger[,id]))
}), use.names=FALSE) 
# cor_Alevin.percell <- unlist(lapply(filter_id_cell, function(id){
#   return(cor(unlist(processed_scDesign2[,id]), processed_Alevin[,id]))
# }), use.names=FALSE) # 517 nonzeros
cor_umitools.percell <- unlist(lapply(filter_id_cell, function(id){
  return(cor(unlist(processed_scDesign2[,id]), processed_umitools[,id]))
}), use.names=FALSE) 
## Kendall's tau
Kcor_cellranger.percell <- unlist(lapply(filter_id_cell, function(id){
  return(cor(unlist(processed_scDesign2[,id]), processed_cellranger[,id], method="kendall"))
}), use.names=FALSE) 
# Kcor_Alevin.percell <- unlist(lapply(filter_id_cell, function(id){
#   return(cor(unlist(processed_scDesign2[,id]), processed_Alevin[,id], method="kendall"))
# }), use.names=FALSE)
Kcor_umitools.percell <- unlist(lapply(filter_id_cell, function(id){
  return(cor(unlist(processed_scDesign2[,id]), processed_umitools[,id], method="kendall"))
}), use.names=FALSE) 
print(mean(cor_cellranger.percell))
# print(mean(cor_Alevin.percell))
print(mean(cor_umitools.percell))
print(mean(Kcor_cellranger.percell))
# print(mean(Kcor_Alevin.percell))
print(mean(Kcor_umitools.percell))
save(cor_cellranger.percell, cor_umitools.percell, Kcor_cellranger.percell, Kcor_umitools.percell, file=sprintf("%s/VerifyPipeline.noAlevin.Quantity.CorrelationPerCell_20221015.Rdata", out_directory))

########################### Quantity. Correlation Global ###########################
cor_cellranger.global <- cor(c(as.matrix(processed_scDesign2)), c(as.matrix(processed_cellranger)))
# cor_Alevin.global <- cor(c(as.matrix(processed_scDesign2)), c(as.matrix(processed_Alevin)))
cor_umitools.global <- cor(c(as.matrix(processed_scDesign2)), c(as.matrix(processed_umitools)))
# Kcor_cellranger.global <- cor(c(as.matrix(processed_scDesign2)), c(as.matrix(processed_cellranger)), method="kendall")
# Kcor_Alevin.global <- cor(c(as.matrix(processed_scDesign2)), c(as.matrix(processed_Alevin)), method="kendall")
# Kcor_umitools.global <- cor(c(as.matrix(processed_scDesign2)), c(as.matrix(processed_umitools)),method="kendall")
# save(cor_cellranger.global, cor_Alevin.global, cor_umitools.global, Kcor_cellranger.percell, Kcor_Alevin.percell, Kcor_umitools.percell, file=sprintf("%s/VerifyPipeline.Quantity.Global.Correlation_20221013.Rdata", out_directory))
save(cor_cellranger.global, cor_umitools.global, file=sprintf("%s/VerifyPipeline.noAlevin.Quantity.Global.Correlation_20221015.Rdata", out_directory))

########################### Fig A. Correlation  ###########################
cor_df <- data.frame(value=c(cor_cellranger, cor_umitools, cor_cellranger.percell,cor_umitools.percell,
                             Kcor_cellranger, Kcor_umitools, Kcor_cellranger.percell,Kcor_umitools.percell),
                     measure=c(rep("gene-wise", 2*length(cor_cellranger)), rep("cell-wise", 2*length(cor_cellranger.percell)),
                               rep("gene-wise", 2*length(cor_cellranger)), rep("cell-wise", 2*length(cor_cellranger.percell))),
                     method=c(rep("cellranger", length(cor_cellranger)),  rep("UMI-tools", length(cor_cellranger)), rep("cellranger", length(cor_cellranger.percell)), rep("UMI-tools", length(cor_cellranger.percell)),
                              rep("cellranger", length(cor_cellranger)),  rep("UMI-tools", length(cor_cellranger)), rep("cellranger", length(cor_cellranger.percell)), rep("UMI-tools", length(cor_cellranger.percell))),
                     cor=c(rep("Pearson cor.", 2 * (length(cor_cellranger) + length(cor_cellranger.percell))),
                           rep("Kendall's tau", 2 * (length(cor_cellranger) + length(cor_cellranger.percell)))))
cor_df$method <- factor(cor_df$method, levels=c("cellranger", "UMI-tools"))
cor_df$measure <- factor(cor_df$measure, levels=c("gene-wise", "cell-wise"))
cor_df$cor <- factor(cor_df$cor, levels=c("Pearson cor.", "Kendall's tau"))

save(cor_df, file=sprintf("%s/BenchmarkUMI.noAlevin.Fig1_Correlation_20221015.Rdata", out_directory))

load(file=sprintf("%s/BenchmarkUMI.noAlevin.Fig1_Correlation_20221015.Rdata", out_directory))
group.colors <- c("#D55E00", "#E69F00", "#56B4E9", "#009E73") # groudtruth, cellranger, Alevin, UMI-tools
cor_plot <- cor_df %>%
  ggplot( aes(x=method, y=value, fill=method)) + 
  stat_boxplot(geom = "errorbar",width=0.3,size=1)+
  geom_boxplot(size=1) +
  facet_grid(measure~cor) +
  theme_bw() +
  scale_fill_manual(values=group.colors[-c(1,3)]) +
  theme(
    # panel.border = element_blank(),
    # axis.line = element_line(colour = "black",size = 1,lineend = "square"),
    # panel.border = element_rect(colour = "black", size=1),
    strip.text = element_text(size=25),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    axis.text.y = element_text(size = 25),
    # axis.text.x=element_blank(),
    # axis.text.x = element_text(size = 25),
    axis.text.x = element_text(size = 25, angle = 45,hjust = 1),
    axis.title=element_text(size=25),
    legend.text=element_text(size=25),
    legend.position="none") +
  ylab("") +
  xlab("")
ggsave(file=sprintf("%s/BenchmarkUMI.noAlevin.Fig1_Correlation_20230319.pdf",out_directory), cor_plot, height=8.4, width = 8.4)


########################### Fig B. Summary statistics ###########################
# two functions for calculating the correlation matrix(-ces) of selected genes ----------
get_stats <- function(mat, group, log_trans = TRUE){
  mean <- rowMeans(mat)
  var <- apply(mat,1,var)
  cv <- sqrt(var)/mean
  zero_gene <- rowSums(mat < 1e-5)/ncol(mat)
  zero_cell <- colSums(mat < 1e-5)/nrow(mat)
  libsize <- colSums(mat)
  
  if(log_trans){
    mean <- log2(mean + 1)
    var <- log2(var + 1)
    libsize <- log2(libsize + 1)
  }
  
  summs <- list(mean = mean, var = var, cv = cv, drop_gene = zero_gene,
                drop_cell = zero_cell, libsize = libsize)
  summs = lapply(1:length(summs), function(i){
    data.frame(value = summs[[i]], measure = names(summs)[i], group = group,
               stringsAsFactors = FALSE)
  })
  summs = Reduce(rbind, summs)
  return(summs)
}

stats_scDesign2 <- get_stats(processed_scDesign2, 'Ground truth')
stats_umitools <- get_stats(processed_umitools , 'UMI-tools')
# stats_alevin <- get_stats(processed_Alevin, 'Alevin')
stats_cellranger <- get_stats(processed_cellranger, 'cellranger')

stats_dat <- rbind(stats_scDesign2, stats_umitools, stats_cellranger)
stats_dat$group <- factor(stats_dat$group, levels = c('Ground truth','cellranger', 'UMI-tools'))
measures1 <-  c("mean", "var", "cv", "drop_gene",
                "drop_cell", "libsize")
measures2 <-  c("gene mean", "gene variance", "gene cv",
                "gene zero prop.", "cell zero prop.", "cell library size")
stats_dat$measure <- factor(stats_dat$measure, levels = measures1)
stats_dat$measure <- mapvalues(stats_dat$measure, from = measures1, to = measures2)
save(stats_dat, file=sprintf("%s/BenchmarkUMI.noAlevin.Fig2_SummaryStats_20221015.Rdata", out_directory))

# create violin-plots to compare the marginal stat values -------------------------------
load(sprintf("%s/BenchmarkUMI.noAlevin.Fig2_SummaryStats_20221015.Rdata", out_directory))
group.colors <- c("#D55E00", "#E69F00", "#56B4E9", "#009E73") # groudtruth, cellranger, Alevin, UMI-tools
stats_plot <- ggplot(stats_dat, aes(x = group, y = value, fill=group)) +
  geom_violin(scale = 'width', trim = TRUE) +
  scale_fill_manual(values=group.colors[-3]) +
  facet_wrap(~measure, scales = "free_y", nrow = 2) +
  theme_bw() +
  theme(
    # panel.border = element_rect(colour = "black", size=1),
    strip.text = element_text(size=25),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    axis.text = element_text(size = 25),
    # axis.text.x=element_blank(),
    axis.text.x = element_text(size = 25, angle = 45,hjust = 1),
    legend.text=element_text(size=25),
    legend.position="none") +
  xlab("") + ylab("")
ggsave(plot=stats_plot, file=sprintf("%s/BenchmarkUMI.noAlevin.Fig2_SummaryStats_20230319.pdf",out_directory), height=8, width = 11)


########################### Fig C. Dimension reduction ###########################
# Important: perform PCA on real matrix, using the same projection onto other synthetic matrices
# Then perform UMAP on real PCA components, then project other synthetic matrices into the real UMAP space
filter_nonzero_id <- function(data_mat, filter_thr=0){
  nonzero_prop <- apply(data_mat, 1, function(x){
    sum(x>0)
  }) / ncol(data_mat)
  nonzero_id <- which(nonzero_prop > filter_thr)
  return(nonzero_id)
}

logtrans <- function(Y){
  return(log(Y+1))
}


# Specify parameters
pca_dim <- 30
filter_thr <- 0
# Filter zero genes
scDesign2.filtered.id <- filter_nonzero_id(processed_scDesign2, filter_thr)
umitools.filtered.id <- filter_nonzero_id(processed_umitools, filter_thr)
cellranger.filtered.id <- filter_nonzero_id(processed_cellranger, filter_thr)
# Keep consensus peaks
gene.filterd.id <- Reduce(intersect, list(scDesign2.filtered.id, umitools.filtered.id, cellranger.filtered.id)) 
scDesign2.processed.matrix <- logtrans(processed_scDesign2[gene.filterd.id,])
umitools.processed.matrix <- logtrans(processed_umitools[gene.filterd.id,])
cellranger.processed.matrix <- logtrans(processed_cellranger[gene.filterd.id,])
save(gene.filterd.id, scDesign2.processed.matrix, umitools.processed.matrix, cellranger.processed.matrix, file=sprintf("%s/BenchmarkUMI.noAlevin.Fig3_UMAP.ManualDR.RealPCASpace.ProcessedData.20230320.Rdata", out_directory))

# PCA
scDesign2_pca <- irlba::irlba(t(scDesign2.processed.matrix), nv=pca_dim)
scDesign2.pcatfidf <- t(scDesign2_pca$u %*% diag(scDesign2_pca$d)) # K by n
umitools.pcatfidf <- t(t(umitools.processed.matrix) %*% scDesign2_pca$v) # K by n
cellranger.pcatfidf <- t(t(cellranger.processed.matrix) %*% scDesign2_pca$v) # K by n
save(scDesign2_pca, scDesign2.pcatfidf, umitools.pcatfidf, cellranger.pcatfidf, file=sprintf("%s/BenchmarkUMI.noAlevin.Fig3_UMAP.ManualDR.RealPCASpace.PCAtfidfprojected.20230320.Rdata", out_directory))

# UMAP
custom.config = umap::umap.defaults
custom.config$random_state = 123
scDesign2_umap_fit <- umap::umap(t(scDesign2.pcatfidf), config=custom.config)
scDesign2.UMAPtfidf <- scDesign2_umap_fit$layout
umitools.UMAPtfidf.projected <- predict(scDesign2_umap_fit, t(umitools.pcatfidf))
cellranger.UMAPtfidf.projected <- predict(scDesign2_umap_fit, t(cellranger.pcatfidf))
save(scDesign2_umap_fit, scDesign2.UMAPtfidf, umitools.UMAPtfidf.projected, cellranger.UMAPtfidf.projected, file=sprintf("%s/BenchmarkUMI.noAlevin.Fig3_UMAP.ManualDR.RealPCASpace.UMAPtfidfprojected.20230320.Rdata", out_directory))

# Dataframe for projected UMAP
umap_dat.projected.tfidf <- data.frame(
  rbind(scDesign2.UMAPtfidf, umitools.UMAPtfidf.projected, cellranger.UMAPtfidf.projected))
colnames(umap_dat.projected.tfidf) <- c('x', 'y')
umap_dat.projected.tfidf$labels <- rep(clustering_synthetic_cluster[intersect_cell_ID], 3)
umap_dat.projected.tfidf$labels <- factor(umap_dat.projected.tfidf$labels, levels=1:length(unique(umap_dat.projected.tfidf$labels)), ordered = TRUE)
umap_dat.projected.tfidf$panels <- c(rep("ground truth", length(clustering_synthetic_cluster[intersect_cell_ID])),
                                     rep("UMI-tools", length(clustering_synthetic_cluster[intersect_cell_ID])),
                                     rep("cellranger", length(clustering_synthetic_cluster[intersect_cell_ID])))
umap_dat.projected.tfidf$panels <- factor(umap_dat.projected.tfidf$panels,
                                          levels = c("ground truth", "cellranger", "UMI-tools"))
save(umap_dat.projected.tfidf, file=sprintf("%s/BenchmarkUMI.noAlevin.Fig3_UMAP.ManualDR.RealPCASpace.UMAPtfidfprojected.20230320.Rdata", out_directory))

# Calculate lisi by catenating low-dimensional cells
library(lisi)
get_lisi_labelset <- function(x, subset_labelset){
  emb_dat <- subset(x, subset = panels %in% subset_labelset,
                    select = c("x", "y"))
  emb_label <- subset(x, subset = panels %in% subset_labelset,
                      select = "panels", drop = FALSE)
  compute_lisi(emb_dat, emb_label, "panels")
}
subset_label_list <- list(c("ground truth", "cellranger"), c("ground truth", "UMI-tools"))
mlisi_umap_mixed <- unlist(lapply(subset_label_list, function(subset_labelset){
  median(unlist(get_lisi_labelset(umap_dat.projected.tfidf, subset_labelset)))
}))
print(mlisi_umap_mixed)
save(mlisi_umap_mixed, file=sprintf("%s/BenchmarkUMI.noAlevin.Fig3_UMAP.ManualDR.RealPCASpace.UMAPtfidfprojected.miLISIvalues.20230320.Rdata", out_directory))

# Plot
load(sprintf("%s/BenchmarkUMI.noAlevin.Fig3_UMAP.ManualDR.RealPCASpace.UMAPtfidfprojected.20230320.Rdata", out_directory))
point_size <- 0.25; point_transp <- 0.8
get_colScale <- function(labels){
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  n_cell_type <- length(levels(labels)) - 2
  # print(str(labels))
  if(n_cell_type > 8)
  {
    colourCount <- n_cell_type
    getPalette <- colorRampPalette(brewer.pal(8, "Dark2"))
    cbPalette <- getPalette(colourCount)
  }else{
    cbPalette <- cbPalette[1:n_cell_type]
  }
  cbPalette <- c(cbPalette, '#FF0000', '#000000')
  names(cbPalette) <- levels(labels)
  scale_colour_manual(name = "labels",values = cbPalette)
}

umap_mixed_plot1 <- umap_dat.projected.tfidf %>% 
  # filter( x < 10 & x > -10 ) %>%
  filter(y > -10 & y < 10 &  x < 10 & x > -10 ) %>%
  ggplot(aes(x = x, y = y, color = labels)) +
  ggrastr::rasterize(geom_point(cex = point_size, alpha = point_transp),dpi = 300) +
  facet_wrap(~panels, nrow = 2) +
  theme_bw() +
  theme(
    # panel.border = element_rect(colour = "black", size=1),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    # legend.title = element_blank(),
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 25),
    legend.text=element_text(size=25),
    strip.text = element_text(size = 25), 
    legend.title=element_text(size=25)) +
  guides(color=guide_legend(override.aes = list(size = 5),
                            ncol=2, byrow=FALSE,
                            title="cell cluster")) +
  get_colScale(umap_dat.projected.tfidf$labels) +
  xlab("UMAP 1") + ylab("UMAP 2")
ggsave(file=sprintf("%s/BenchmarkUMI.noAlevin.Fig3_UMAP.ManualDR.RealPCASpace.UMAPtfidfprojected.noLISI.20230320.pdf",out_directory), umap_mixed_plot1, height=9.4, width = 9.4)

# PCA + UMAP individually
# # library(dplyr)
# # library(Seurat)
# # library(patchwork)
# 
# Manual_dimreduc <- function(matrix, filter_thr){
#   filter_nonzero <- function(data_mat, filter_thr=0){
#     nonzero_prop <- apply(data_mat, 1, function(x){
#       sum(x>0)
#     }) / ncol(data_mat)
#     data_mat_filtered <- data_mat[nonzero_prop > filter_thr,]
#     return(data_mat_filtered)
#   }
#   logtrans <- function(Y){
#     return(log(Y+1))
#   }
#   pca_fun <- function(Y, pca_dim=30){
#     res_pca <- irlba::irlba(t(Y), nv=pca_dim)
#     score_pca <- res_pca$u %*% diag(res_pca$d)
#     return(t(score_pca))
#   }
#   matrix_filtered <- filter_nonzero(matrix, filter_thr=filter_thr)
#   pca_cluster_label_combined <-
#     pca_fun(logtrans(matrix_filtered))
#   return(pca=pca_cluster_label_combined)
# }
# 
# # Seurat: Extract submatrix and cluster results for umap
# # real_seurat_list <- Seurat_cluster(t(matrix_num))
# # umitools_seurat_list <- Seurat_cluster(t(matrix_umitools))
# # alevin_seurat_list <- Seurat_cluster(t(matrix_Alevin))
# # real_clusters <- real_seurat_list$clusters
# # simu_umi_clusters <- umitools_seurat_list$clusters
# # simu_alevin_clusters <- alevin_seurat_list$clusters
# # real_pca <- real_seurat_list$pca
# # simu_umi_pca <- umitools_seurat_list$pca
# # simu_alevin_pca <- alevin_seurat_list$pca
# # save(real_seurat_list, umitools_seurat_list, alevin_seurat_list, file=sprintf("%s/VerifyCount_Fig4_UMAP.Seurat.clustering.Rdata", umitool.out.directory))
# 
# # Manual: Extract submatrix and cluster results for umap
# scDesign2_dr_list <- Manual_dimreduc(processed_scDesign2, 0)
# umitools_dr_list <- Manual_dimreduc(processed_umitools, 0)
# # alevin_dr_list <- Manual_dimreduc(processed_Alevin, 0)
# cellranger_dr_list <- Manual_dimreduc(processed_cellranger, 0)
# 
# scDesign2_pca <- scDesign2_dr_list
# simu_umi_pca <- umitools_dr_list
# # simu_alevin_pca <- alevin_dr_list
# simu_cellranger_pca <- cellranger_dr_list
# 
# # UMAP
# custom.config = umap::umap.defaults
# custom.config$random_state = 123
# scDesign2_umap_fit <- umap::umap(t(scDesign2_pca), config=custom.config)
# umap_result_scDesign2_data <- scDesign2_umap_fit$layout
# # umap_result_simu_umi <- predict(scDesign2_umap_fit, t(simu_umi_pca))
# # umap_result_simu_alevin <- predict(scDesign2_umap_fit, t(simu_alevin_pca))
# umap_result_simu_umi <- umap::umap(t(simu_umi_pca), config=custom.config)$layout
# # umap_result_simu_alevin <- umap::umap(t(simu_alevin_pca), config=custom.config)$layout
# umap_result_simu_cellranger <- umap::umap(t(simu_cellranger_pca), config=custom.config)$layout
# 
# # Rotate some figures for a better comparison
# umap_result_simu_cellranger_rotate <- cbind(umap_result_simu_cellranger[,1], -umap_result_simu_cellranger[,2])
# umap_result_simu_umi_rotate <- cbind(-umap_result_simu_umi[,2], umap_result_simu_umi[,1])
# # umap_result_simu_alevin_rotate <- cbind(-umap_result_simu_alevin[,2], -umap_result_simu_alevin[,1])
# 
# # Prepare dataframe for plot
# umap_dat <- data.frame(
#   rbind(umap_result_scDesign2_data, umap_result_simu_umi, umap_result_simu_cellranger))
# colnames(umap_dat) <- c('x', 'y')
# umap_dat$labels <- rep(clustering_synthetic_cluster[intersect_cell_ID], 3)
# umap_dat$labels <- factor(umap_dat$labels, levels=1:length(unique(umap_dat$labels)), ordered = TRUE)
# umap_dat$panels <- c(rep("Ground truth", length(clustering_synthetic_cluster[intersect_cell_ID])),
#                      rep("UMI-tools", length(clustering_synthetic_cluster[intersect_cell_ID])),
#                      rep("cellranger", length(clustering_synthetic_cluster[intersect_cell_ID])))
# umap_dat$panels <- factor(umap_dat$panels,
#                           levels = c("Ground truth", "cellranger", "UMI-tools")) 
# # Prepare dataframe for rotated plot
# umap_dat_rotate <- data.frame(
#   rbind(umap_result_scDesign2_data, umap_result_simu_umi_rotate, umap_result_simu_cellranger_rotate))
# colnames(umap_dat_rotate) <- c('x', 'y')
# umap_dat_rotate$labels <- rep(clustering_synthetic_cluster[intersect_cell_ID], 3)
# umap_dat_rotate$labels <- factor(umap_dat_rotate$labels, levels=1:length(unique(umap_dat_rotate$labels)), ordered = TRUE)
# umap_dat_rotate$panels <- c(rep("Ground truth", length(clustering_synthetic_cluster[intersect_cell_ID])),
#                             rep("UMI-tools", length(clustering_synthetic_cluster[intersect_cell_ID])),
#                             rep("cellranger", length(clustering_synthetic_cluster[intersect_cell_ID])))
# umap_dat_rotate$panels <- factor(umap_dat_rotate$panels,
#                                  levels = c("Ground truth", "cellranger", "UMI-tools"))
# 
# # Add global Pearson correlation to UMAP
# lisi_panel_vec <- c(levels(umap_dat_rotate$panels)[1], 
#                     paste0(levels(umap_dat_rotate$panels)[-1], "\n",
#                            "[Pearson Cor. = ",
#                            format(round(c(cor_cellranger.global, cor_umitools.global), 3), nsmall = 3), "]"))
# umap_dat_rotate$panels <- mapvalues(umap_dat_rotate$panels, from = levels(umap_dat_rotate$panels),
#                                     to = lisi_panel_vec)
# 
# # save(umap_dat, file=sprintf("%s/VerifyCount_Fig4_UMAP.SeuratDR.noguide.Rdata", umitool.out.directory))
# # save(umap_dat, file=sprintf("%s/VerifyCount_Fig4_UMAP.ManualDR.filter03.noguide.Rdata", umitool.out.directory))
# save(umap_dat, file=sprintf("%s/BenchmarkUMI.noAlevin.Fig3_UMAP.ManualDR.fixseed_20221015.Rdata", out_directory))
# save(umap_dat_rotate, file=sprintf("%s/BenchmarkUMI.noAlevin.Fig3_UMAP.ManualDR.fixseed.rotated_20221015.Rdata", out_directory))
# 
# 
# ### plotting
# load(file=sprintf("%s/BenchmarkUMI.noAlevin.Fig3_UMAP.ManualDR.fixseed.rotated_20221015.Rdata", out_directory))
# 
# point_size <- 0.25; point_transp <- 0.8
# get_colScale <- function(labels){
#   cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
#                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#   n_cell_type <- length(levels(labels)) - 2
#   # print(str(labels))
#   if(n_cell_type > 8)
#   {
#     colourCount <- n_cell_type
#     getPalette <- colorRampPalette(brewer.pal(8, "Dark2"))
#     cbPalette <- getPalette(colourCount)
#   }else{
#     cbPalette <- cbPalette[1:n_cell_type]
#   }
#   cbPalette <- c(cbPalette, '#FF0000', '#000000')
#   names(cbPalette) <- levels(labels)
#   scale_colour_manual(name = "labels",values = cbPalette)
# }
# 
# umap_mixed_plot1 <- umap_dat_rotate %>% 
#   # filter( x < 10 & x > -10 ) %>%
#   filter(y > -10 & y < 10 &  x < 10 & x > -10 ) %>%
#   ggplot(aes(x = x, y = y, color = labels)) +
#   ggrastr::rasterize(geom_point(cex = point_size, alpha = point_transp),dpi = 300) +
#   facet_wrap(~panels, nrow = 1) +
#   theme_bw() +
#   theme(
#     # panel.border = element_rect(colour = "black", size=1),
#     plot.title = element_text(hjust = 0.5),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     legend.position = "bottom",
#     # legend.title = element_blank(),
#     axis.text = element_text(size = 25),
#     axis.title = element_text(size = 25),
#     legend.text=element_text(size=25),
#     strip.text = element_text(size = 25), 
#     legend.title=element_text(size=25)) +
#   guides(color=guide_legend(override.aes = list(size = 5),
#                             nrow=2, byrow=FALSE,
#                             title="cell cluster")) +
#   get_colScale(umap_dat_rotate$labels) +
#   xlab("UMAP 1") + ylab("UMAP 2")
# ggsave(file=sprintf("%s/BenchmarkUMI.noAlevin.Fig3_UMAP.ManualDR.fixseed.rotated_20230315.pdf",out_directory), umap_mixed_plot1,height=5.5, width = 12)
# # ggsave(file=sprintf("%s/VerifyCount_Fig4_UMAP.SeuratDR.noguide.pdf",umitool.out.directory), umap_mixed_plot1, height = 4.8, width = 9.6, dpi = 600)


########################### Fig D. Time Complexity ###########################
varyCell <- c(1219, 2438, 4877, 9754, 19508)
varySeq <- c("5M", "10M", "20M", "40M", "82M")

cellranger.varyCell <- c(3825, 3936, 3930, 3858, 3831)
cellranger.varySeq <- c(1671, 2235, 5118, 6722, 7007)
umitools.varyCell <- c(10104, 10350, 10428, 10321, 10544)
umitools.varySeq <- c(2603, 5191, 10198, 20519, 39972)

time_dat <- data.frame(time=c(cellranger.varyCell, cellranger.varySeq,
                              umitools.varyCell, umitools.varySeq),
                       case=rep(c(varyCell, varySeq), 2),
                       method=c(rep("cellranger", 10),
                                rep("UMI-tools", 10)),
                       measure=rep(c(rep("Cell number", 5), rep("Seq. depth", 5)), 2))
time_dat$method <- factor(time_dat$method, levels=c("cellranger", "UMI-tools"))
time_dat$measure <- factor(time_dat$measure, levels=c("Cell number", "Seq. depth"))
time_dat$case <- factor(time_dat$case, levels=c(varyCell, varySeq))
save(time_dat, file=sprintf("%s/BenchmarkUMI_Fig4_TIME.noAlevin.20230315.Rdata", out_directory))

load(file=sprintf("%s/BenchmarkUMI_Fig4_TIME.noAlevin.20230315.Rdata", out_directory))
group.colors <- c("#D55E00", "#E69F00", "#56B4E9", "#009E73") # groudtruth, cellranger, Alevin, UMI-tools
time_plot <- time_dat %>%
  ggplot( aes(x=case, y=time, color=method, shape=method, group=method)) +
  scale_y_log10() +
  facet_wrap(~measure, ncol=2, scales = "free_x") +
  geom_point(size=6) +
  geom_line(aes(color=method), size=1) +
  theme_bw() +
  theme(
    # axis.line = element_line(colour = "black",size = 1,lineend = "square"),
    # panel.border = element_rect(colour = "black", size=1),
    strip.text = element_text(size=25),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    axis.text.y = element_text(size = 25),
    # axis.text.x=element_blank(),
    # axis.text.x = element_text(size = 12),
    axis.text.x = element_text(size = 25, angle = 45,hjust = 1),
    axis.title=element_text(size=25),
    legend.text=element_text(size=25),
    legend.position="bottom") +
  scale_shape_manual(labels = c("cellranger", "UMI-tools"),
                     values = c(16, 17)) +
  scale_color_manual(labels=c("cellranger", "UMI-tools"), values=group.colors[-c(1,3)]) +
  scale_size_manual(values=c(5,5)) + 
  ylab("Elapsed time (sec)") +
  xlab("")
ggsave(file=sprintf("%s/BenchmarkUMI.noAlevin.Fig4_Time.20230315.pdf",out_directory), time_plot, height=6.5, width = 10)





### Combine
# gs <- list(arrangeGrob(cor_plot, left = textGrob("A", x = unit(1, "npc"), 
#                                                    y = unit(.95, "npc"))),
#            arrangeGrob(stats_plot, left = textGrob("B", x = unit(1, "npc"), 
#                                                    y = unit(.95, "npc"))),
#            arrangeGrob(umap_mixed_plot1, left = textGrob("C", x = unit(1, "npc"), 
#                                                          y = unit(.95, "npc"))))
# lay <- rbind(c(1,1,2,2),
#              c(1,1,2,2),
#              c(3,3,3,3),
#              c(3,3,3,3),
#              c(4,4,4,4),
#              c(4,4,4,4),
#              c(4,4,4,4))

# With Alevin: fill in time complexity later
gs <- list(arrangeGrob(time_plot, left = textGrob("A", x = unit(1, "npc"),
                                                  y = unit(.95, "npc"))),
           arrangeGrob(cor_plot, left = textGrob("B", x = unit(1, "npc"),
                                                 y = unit(.95, "npc"))),
           arrangeGrob(stats_plot, left = textGrob("C", x = unit(1, "npc"),
                                                   y = unit(.95, "npc"))),
           arrangeGrob(umap_mixed_plot1, left = textGrob("D", x = unit(1, "npc"),
                                                         y = unit(.95, "npc"))))
lay <- rbind(c(1,1,2,2),
             c(1,1,2,2),
             c(1,1,2,2),
             c(3,3,3,3),
             c(3,3,3,3),
             c(3,3,3,3),
             c(4,4,4,4),
             c(4,4,4,4),
             c(4,4,4,4))
p_final <-arrangeGrob(grobs = gs, layout_matrix = lay)
# ggsave(file=sprintf("%s/BenchmarkUMI.noAlevin.Combined_20221108.pdf",out_directory), p_final, height = 9.6, width = 9.6, dpi = 600)
ggsave(file=sprintf("%s/BenchmarkUMI.noAlevin.Combined_20230224.pdf",out_directory), p_final, height = 9.6, width = 9.6, dpi = 600)

