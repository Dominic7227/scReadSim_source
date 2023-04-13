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

out_directory <- "/home/gayan/Projects/scATAC_Simulator/results/20230304_Minnow"
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
# scDesign2_genecode <- unlist(lapply(sapply(scDesign2_genecode, strsplit, split=".", fixed = TRUE), function(x) return(x[1])), use.names = FALSE)

# ## Construct gene geneid list
# # Gene region could overlap (different strands)
# # scDesign2 count matrix row refers to merged genes, hence one row links to multiple genes potentially 
# chr1_geneid <-  unlist(lapply(sapply(chr1_genes[,9], strsplit, split=";", fixed=TRUE), function(x) {
#   return(strsplit(x[1], split=" ",fixed=TRUE)[[1]][2])
# }), use.names=FALSE)
# chr1_geneid <- unlist(lapply(sapply(chr1_geneid, strsplit, split=".", fixed = TRUE), function(x) return(x[1])), use.names = FALSE)
# chr1_genename <- unlist(lapply(sapply(chr1_genes[,9], strsplit, split=";", fixed=TRUE), function(x) {
#   return(strsplit(x[4], split=" ",fixed=TRUE)[[1]][3])
# }), use.names=FALSE)
# chr1_gene_geneid_list <- cbind(chr1_geneid, chr1_genename)
# # Create a gene to region index list. Same length to chr1 gene list
# scReadSim_gene2region_indice <- rep(1:length(scDesign2_genecode), c(diff(match(scDesign2_genes_merged[,2], chr1_genes[,4])), 1))
# # write.table(chr1_gene_geneid_list[,1], file="/home/gayan/Projects/scATAC_Simulator/results/20221011_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/verification_alevin/ref_data_chr1/chr1_geneID.txt",col.names=FALSE, row.names=FALSE,quote = FALSE)


# # Match scDesign2 meta-gene gene_id to a random transcript_id 
# # Select 200 genes for analysis
# ngene_scDesign2 <- nrow(scDesign2_matrix)
# # Obatin the transcript and gene id list
# library(rtracklayer)
# my_obj <- import("/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE/gencode.vM10.annotation.gtf")
# chr1_obj <- my_obj %>% as.data.frame() %>% filter(seqnames == "chr1") %>% dplyr::select(c(gene_id, transcript_id))
# rm(my_obj)
# # Find a random transcript id for each gene id
# set.seed(2023)
# assigned_transcript_id_list <- vector("integer", length=ngene_scDesign2)
# for (id in 1:ngene_scDesign2){
#   tmp <- unlist(drop_na(unique(chr1_obj %>% filter(gene_id == scDesign2_genecode[id]) %>% dplyr::select(transcript_id))))
#   if (length(tmp) > 0) {
#     assigned_transcript_id_list[id] <- sample(tmp)[1]
#   } else {
#     print(sprintf("No Transcript for Gene: %s", scDesign2_genecode[id]))
#   }
# }
# fileConn<-file(sprintf("%s/ground_truth_scReadSim_countmat/quants_mat_rows.txt", out_directory))
# writeLines(assigned_transcript_id_list, fileConn)
# close(fileConn)

fileConn<-file(sprintf("%s/ground_truth_scReadSim_countmat/quants_mat_rows.txt", out_directory))
writeLines(scDesign2_genecode, fileConn)
close(fileConn)
write.table(matrix_scDesign2, sprintf("%s/ground_truth_scReadSim_countmat/quants_mat.csv", out_directory), quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")


