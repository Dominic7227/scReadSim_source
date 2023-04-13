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

directory <- "/home/gayan/Projects/scATAC_Simulator/results/20230304_SCANATACSim"
out_directory <- "/home/gayan/Projects/scATAC_Simulator/results/20230304_SCANATACSim/verification_output"
dir.create(out_directory)

####################################################################################
############################ Read in Real Count Matrix file #######################
####################################################################################
real.barcode.file <- sprintf("%s/CBwithCellType.noCBZ.HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.Real.txt",directory) 
real_celltype <- unlist(read.table(real.barcode.file, sep="\t",header = FALSE)[,2])
real.count.mat <- read.table(sprintf("%s/HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.Real.CountMatrix.txt", directory), sep="\t",header = FALSE)
peak_vec <- unlist(real.count.mat[,1])
real.count.mat <- real.count.mat[,-1]
####################################################################################
############################ Read in scReadSim Count Matrix file ###################
####################################################################################
scReadSim.barcode.file <- sprintf("%s/CBwithCellType.noCBZ.HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.scReadSim.txt",directory) 
scReadSim_celltype <- unlist(read.table(scReadSim.barcode.file, sep="\t",header = FALSE)[,2])
scReadSim.count.mat <- read.table(sprintf("%s/HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.scReadSim.CountMatrix.txt", directory), sep="\t",header = FALSE)
scReadSim.count.mat <- scReadSim.count.mat[,-1]
####################################################################################
######################## Read in SCAN-ATAC-Sim Count Matrix files #################
####################################################################################
SCANATACSim_celltype <- c(rep("hematopoietic Progenitors", 1621), rep("erythroblasts", 1153), rep("monocytes", 531), rep("immature B cells", 296))
SCANATACSim_HematopoieticProgenitors.count.mat <- read.table(sprintf("%s/SCANATACSim_output/Hematopoieticprogenitors.SCAN-ATAC-Sim.CountMatrix.txt", directory), sep="\t",header = FALSE)
SCANATACSim_Erythroblasts.count.mat <- read.table(sprintf("%s/SCANATACSim_output/Erythroblasts.SCAN-ATAC-Sim.CountMatrix.txt", directory), sep="\t",header = FALSE)
SCANATACSim_Monocytes.count.mat <- read.table(sprintf("%s/SCANATACSim_output/Monocytes.SCAN-ATAC-Sim.CountMatrix.txt", directory), sep="\t",header = FALSE)
SCANATACSim_ImmatureBcells.count.mat <- read.table(sprintf("%s/SCANATACSim_output/ImmatureBcells.SCAN-ATAC-Sim.CountMatrix.txt", directory), sep="\t",header = FALSE)
SCANATACSim.count.mat <- cbind(SCANATACSim_HematopoieticProgenitors.count.mat[,-1], SCANATACSim_Erythroblasts.count.mat[,-1], SCANATACSim_Monocytes.count.mat[,-1], SCANATACSim_ImmatureBcells.count.mat[,-1])
save(real.count.mat, scReadSim.count.mat, SCANATACSim.count.mat, real_celltype, scReadSim_celltype, SCANATACSim_celltype, file=sprintf("%s/Comparison_CountMatrix.Rdata", out_directory))

load(sprintf("%s/Comparison_CountMatrix.Rdata", out_directory))
########################### Fig Palettes Setting ###########################
method_pallete <- c("#F8766D", "#00BFC4", "#C77CFF") # wesandersan palette: Royal2/ Method: real, scReadSim, SCANATACSim
cell_type_pallete <- c("#899DA4", "#C93312", "#FAEFD1", "#DC863B") # wesandersan palette: Royal1
########################### Fig A. Summary statistics ###########################
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

stats.Real <- get_stats(real.count.mat, 'real')
stats.scReadSim <- get_stats(scReadSim.count.mat , 'scReadSim')
stats.SCANATACSim <- get_stats(SCANATACSim.count.mat, 'SCAN-ATAC-Sim')

stats_dat <- rbind(stats.Real, stats.scReadSim, stats.SCANATACSim)
stats_dat$group <- factor(stats_dat$group, levels = c('real','scReadSim', 'SCAN-ATAC-Sim'))
measures1 <-  c("mean", "var", "cv", "drop_gene",
                "drop_cell", "libsize")
measures2 <-  c("peak mean", "peak variance", "peak cv",
                "peak zero prop.", "cell zero prop.", "cell library size")
stats_dat$measure <- factor(stats_dat$measure, levels = measures1)
stats_dat$measure <- mapvalues(stats_dat$measure, from = measures1, to = measures2)
save(stats_dat, file=sprintf("%s/Fig1_SummaryStats_20230310.Rdata", out_directory))

# # Method KS test: Both distinct from null
# scReadSim
KS_test_scReadSim <- unlist(lapply(1:6, function(summary_statistic_id){
  real_samples <- unlist(stats.Real %>% filter(measure==measures1[summary_statistic_id]) %>% select(value))
  synthetic_samples <- unlist(stats.scReadSim %>% filter(measure==measures1[summary_statistic_id]) %>% select(value))
  return(ks.test(synthetic_samples, real_samples)$statistic)
}))
names(KS_test_scReadSim) <- measures1
# SCAN-ATAC-Sim
KS_test_SCANATACSim <- unlist(lapply(1:6, function(summary_statistic_id){
  real_samples <- unlist(stats.Real %>% filter(measure==measures1[summary_statistic_id]) %>% select(value))
  synthetic_samples <- unlist(stats.SCANATACSim %>% filter(measure==measures1[summary_statistic_id]) %>% select(value))
  return(ks.test(synthetic_samples, real_samples)$statistic)
}))
names(KS_test_SCANATACSim) <- measures1
save(KS_test_scReadSim, KS_test_SCANATACSim, file=sprintf("%s/Fig1_SummaryStats_KS_Statistics_20230313.Rdata", out_directory))
format(round(KS_test_scReadSim, 3), nsmall = 3)
format(round(KS_test_SCANATACSim, 3), nsmall = 3)

# create violin-plots to compare the marginal stat values -------------------------------
load(sprintf("%s/Fig1_SummaryStats_20230310.Rdata", out_directory))
group.colors <- method_pallete # groudtruth, scReadSim, SCAN-ATAC-Sim
stats_plot <- ggplot(stats_dat, aes(x = group, y = value, fill=group)) +
  geom_violin(scale = 'width', trim = TRUE) +
  scale_fill_manual(values=group.colors) +
  facet_wrap(~measure, scales = "free", ncol = 3) +
  theme_bw() +
  theme(
    strip.text = element_text(size=15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title=element_blank(),
    axis.text = element_text(size = 12),
    axis.text.x=element_blank(),
    legend.text=element_text(size=15),
    legend.position="right") +
  xlab("") + ylab("")
ggsave(plot=stats_plot, file=sprintf("%s/Fig1_SummaryStats_20230312.pdf",out_directory), height = 4.8, width = 9.6, dpi = 600)


########################### Fig B. Dimension reduction ###########################
# Important: perform PCA on real matrix, using the same projection onto other synthetic matrices
# Then perform UMAP on real PCA components, then project other synthetic matrices into the real UMAP space
filter_nonzero_id <- function(data_mat, filter_thr=0){
  nonzero_prop <- apply(data_mat, 1, function(x){
    sum(x>0)
  }) / ncol(data_mat)
  nonzero_id <- which(nonzero_prop > filter_thr)
  return(nonzero_id)
}
tf_idf <- function(Y){
  frequences <- colSums(Y)
  nfreqs <- t(apply(Y, 1, function(x){x/frequences}))
  nfreqs[is.na(nfreqs)] <- 0
  idf <- log(1 + ncol(Y)) - log(rowSums(Y > 0) + 1) + 1
  Y_idf <- apply(nfreqs, 2, function(x){x * idf})
  return(Y_idf)
}
logtrans <- function(Y){
  return(log(Y+1))
}


# Specify parameters
pca_dim <- 30
filter_thr <- 0.01
# Filter zero peaks
real.filtered.id <- filter_nonzero_id(real.count.mat, filter_thr)
scReadSim.filtered.id <- filter_nonzero_id(scReadSim.count.mat, filter_thr)
SCANATACSim.filtered.id <- filter_nonzero_id(SCANATACSim.count.mat, filter_thr)
# Keep consensus peaks
peak.filterd.id <- Reduce(intersect, list(real.filtered.id, scReadSim.filtered.id, SCANATACSim.filtered.id))
real.processed.matrix <- logtrans(tf_idf(real.count.mat[peak.filterd.id,]))
scReadSim.processed.matrix <- logtrans(tf_idf(scReadSim.count.mat[peak.filterd.id,]))
SCANATACSim.processed.matrix <- logtrans(tf_idf(SCANATACSim.count.mat[peak.filterd.id,]))
save(peak.filterd.id, real.processed.matrix, scReadSim.processed.matrix, SCANATACSim.processed.matrix, file=sprintf("%s/Fig2_DimensionReduction.RealPCASpace.ProcessedData.20230313.Rdata", out_directory))

# PCA
res_pca <- irlba::irlba(t(real.processed.matrix), nv=pca_dim)
real.pcatfidf <- t(res_pca$u %*% diag(res_pca$d)) # K by n
scReadSim.pcatfidf <- t(t(scReadSim.processed.matrix) %*% res_pca$v) # K by n
SCANATACSim.pcatfidf <- t(t(SCANATACSim.processed.matrix) %*% res_pca$v) # K by n
save(res_pca, real.pcatfidf, scReadSim.pcatfidf, SCANATACSim.pcatfidf, file=sprintf("%s/Fig2_DimensionReduction.RealPCASpace.PCAtfidfprojected.20230313.Rdata", out_directory))

# UMAP
custom.config = umap::umap.defaults
custom.config$random_state = 123
real_umap_fit <- umap::umap(t(real.pcatfidf), config=custom.config)
real.UMAPtfidf <- real_umap_fit$layout
scReadSim.UMAPtfidf.projected <- predict(real_umap_fit, t(scReadSim.pcatfidf))
SCANATACSim.UMAPtfidf.projected <- predict(real_umap_fit, t(SCANATACSim.pcatfidf))
save(real_umap_fit, real.UMAPtfidf, scReadSim.UMAPtfidf.projected, SCANATACSim.UMAPtfidf.projected, file=sprintf("%s/Fig2_DimensionReduction.RealPCASpace.UMAPtfidfprojected.20230313.Rdata", out_directory))

# Dataframe for projected UMAP
umap_dat.projected.tfidf <- data.frame(
  rbind(real.UMAPtfidf, scReadSim.UMAPtfidf.projected, SCANATACSim.UMAPtfidf.projected))
colnames(umap_dat.projected.tfidf) <- c('x', 'y')
umap_dat.projected.tfidf$labels <- c(real_celltype, scReadSim_celltype, SCANATACSim_celltype)
umap_dat.projected.tfidf$labels <- factor(umap_dat.projected.tfidf$labels, levels=c("hematopoietic progenitors", "erythroblasts", "monocytes", "immature B cells"), ordered = TRUE)
umap_dat.projected.tfidf$panels <- c(rep("real", length(real_celltype)),
                                     rep("scReadSim", length(scReadSim_celltype)),
                                     rep("SCAN-ATAC-Sim", length(SCANATACSim_celltype)))
umap_dat.projected.tfidf$panels <- factor(umap_dat.projected.tfidf$panels,
                                          levels = c("real", "scReadSim", "SCAN-ATAC-Sim")) 
save(umap_dat.projected.tfidf, file=sprintf("%s/Fig2_DimensionReduction.RealPCASpace.UMAPtfidfprojected.dataframe.20230313.Rdata", out_directory))


# Plot
load(file=sprintf("%s/Fig2_DimensionReduction.RealPCASpace.UMAPtfidfprojected.dataframe.20230313.Rdata", out_directory))

point_size <- 0.25; point_transp <- 0.8
cbPalette <- cell_type_pallete  # wesanderson Royal1
names(cbPalette) <- levels(umap_dat.projected.tfidf$labels)

umap_plot.tfidf <- umap_dat.projected.tfidf %>% 
  # filter( x < 10 & x > -10 ) %>%
  # filter(y > -10 & y < 10 &  x < 10 & x > -10 ) %>%
  ggplot(aes(x = x, y = y, color = labels)) +
  ggrastr::rasterize(geom_point(cex = point_size, alpha = point_transp),dpi = 300) +
  facet_wrap(~panels, nrow = 1) +
  theme_bw() +
  theme(
    # panel.border = element_rect(colour = "black", size=1),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    # legend.title = element_blank(),
    axis.text = element_text(size = 25),
    axis.title = element_text(size = 25),
    legend.text=element_text(size=25),
    strip.text = element_text(size = 25), 
    legend.title=element_text(size=25)) +
  guides(color=guide_legend(override.aes = list(size = 5),
                            nrow=2, byrow=FALSE,
                            title="cell type")) +
  # get_colScale(umap_dat.tfidf$labels) +
  scale_colour_manual(values = cbPalette) + # wesanderson Royal1
  xlab("UMAP 1") + ylab("UMAP 2")
ggsave(file=sprintf("%s/Fig2_DimensionReduction.RealPCASpace.UMAPtfidfprojected.dataframe.20230315.pdf",out_directory), umap_plot.tfidf, height=5.5, width = 13)


# # # Manual: PCA
# # Real.pca <- Manual_dimreduc(real.count.mat, 0.01)
# # scReadSim.pca <- Manual_dimreduc(scReadSim.count.mat, 0.01)
# # SCANATACSim.pca <- Manual_dimreduc(SCANATACSim.count.mat, 0.01)
# # save(Real.pca, scReadSim.pca, SCANATACSim.pca, file=sprintf("%s/Fig2_DimensionReduction.PCA.20230310.Rdata", out_directory))
# # Manual: TFIDF + PCA
# real.pcatfidf <- Manual_dimreduc(real.count.mat, 0.01)
# scReadSim.pcatfidf <- Manual_dimreduc(scReadSim.count.mat, 0.01)
# SCANATACSim.pcatfidf <- Manual_dimreduc(SCANATACSim.count.mat, 0.01)
# save(real.pcatfidf, scReadSim.pcatfidf, SCANATACSim.pcatfidf, file=sprintf("%s/Fig2_DimensionReduction.PCAtfidf.20230310.Rdata", out_directory))
# 
# # # UMAP
# # custom.config = umap::umap.defaults
# # custom.config$random_state = 123
# # real_umap_fit <- umap::umap(t(Real.pca), config=custom.config)
# # Real.UMAP <- real_umap_fit$layout
# # scReadSim.UMAP <- umap::umap(t(scReadSim.pca), config=custom.config)$layout
# # SCANATACSim.UMAP <- umap::umap(t(SCANATACSim.pca), config=custom.config)$layout
# # save(Real.UMAP, scReadSim.UMAP, SCANATACSim.UMAP, file=sprintf("%s/Fig2_DimensionReduction.UMAP.20230310.Rdata", out_directory))
# # scReadSim.UMAP.projected <- predict(real_umap_fit, t(scReadSim.pca))
# # SCANATACSim.UMAP.projected <- predict(real_umap_fit, t(SCANATACSim.pca))
# # save(Real.UMAP, scReadSim.UMAP.projected, SCANATACSim.UMAP.projected, file=sprintf("%s/Fig2_DimensionReduction.UMAPprojected.20230310.Rdata", out_directory))
# # UMAP (TFIDF)
# custom.config = umap::umap.defaults
# custom.config$random_state = 123
# real_umap_fit <- umap::umap(t(real.pcatfidf), config=custom.config)
# real.UMAPtfidf <- real_umap_fit$layout
# scReadSim.UMAPtfidf <- umap::umap(t(scReadSim.pcatfidf), config=custom.config)$layout
# SCANATACSim.UMAPtfidf <- umap::umap(t(SCANATACSim.pcatfidf), config=custom.config)$layout
# save(real.UMAPtfidf, scReadSim.UMAPtfidf, SCANATACSim.UMAPtfidf, file=sprintf("%s/Fig2_DimensionReduction.UMAPtfidf.20230310.Rdata", out_directory))
# scReadSim.UMAPtfidf.projected <- predict(real_umap_fit, t(scReadSim.pcatfidf))
# SCANATACSim.UMAPtfidf.projected <- predict(real_umap_fit, t(SCANATACSim.pcatfidf))
# save(real.UMAPtfidf, scReadSim.UMAPtfidf.projected, SCANATACSim.UMAPtfidf.projected, file=sprintf("%s/Fig2_DimensionReduction.UMAPtfidfprojected.20230310.Rdata", out_directory))
# 
# # Rotate some figures for a better comparison
# scReadSim.UMAPtfidf.rotated <- cbind(-scReadSim.UMAPtfidf[,1], -scReadSim.UMAPtfidf[,2])
# 
# # # Prepare dataframe for plot
# # umap_dat <- data.frame(
# #   rbind(Real.UMAP, scReadSim.UMAP, SCANATACSim.UMAP))
# # colnames(umap_dat) <- c('x', 'y')
# # umap_dat$labels <- c(real_celltype, scReadSim_celltype, SCANATACSim_celltype)
# # umap_dat$labels <- factor(umap_dat$labels, levels=c("Hematopoietic progenitors", "Erythroblasts", "Monocytes", "Immature B cells"), ordered = TRUE)
# # umap_dat$panels <- c(rep("real", length(real_celltype)),
# #                      rep("scReadSim", length(scReadSim_celltype)),
# #                      rep("SCAN-ATAC-Sim", length(SCANATACSim_celltype)))
# # umap_dat$panels <- factor(umap_dat$panels,
# #                           levels = c("real", "scReadSim", "SCAN-ATAC-Sim")) 
# # save(umap_dat, file=sprintf("%s/Fig2_DimensionReduction.UMAP.dataframe.20230310.Rdata", out_directory))
# # # Dataframe for projected UMAP
# # umap_dat.projected <- data.frame(
# #   rbind(Real.UMAP, scReadSim.UMAP.projected, SCANATACSim.UMAP.projected))
# # colnames(umap_dat.projected) <- c('x', 'y')
# # umap_dat.projected$labels <- c(real_celltype, scReadSim_celltype, SCANATACSim_celltype)
# # umap_dat.projected$labels <- factor(umap_dat.projected$labels, levels=c("Hematopoietic progenitors", "Erythroblasts", "Monocytes", "Immature B cells"), ordered = TRUE)
# # umap_dat.projected$panels <- c(rep("real", length(real_celltype)),
# #                      rep("scReadSim", length(scReadSim_celltype)),
# #                      rep("SCAN-ATAC-Sim", length(SCANATACSim_celltype)))
# # umap_dat.projected$panels <- factor(umap_dat$panels,
# #                           levels = c("real", "scReadSim", "SCAN-ATAC-Sim")) 
# # save(umap_dat.projected, file=sprintf("%s/Fig2_DimensionReduction.UMAPprojected.dataframe.20230310.Rdata", out_directory))
# 
# # Prepare dataframe for plot (TFIDF)
# umap_dat.tfidf <- data.frame(
#   rbind(real.UMAPtfidf, scReadSim.UMAPtfidf.rotated, SCANATACSim.UMAPtfidf))
# colnames(umap_dat.tfidf) <- c('x', 'y')
# umap_dat.tfidf$labels <- c(real_celltype, scReadSim_celltype, SCANATACSim_celltype)
# umap_dat.tfidf$labels <- factor(umap_dat.tfidf$labels, levels=c("Hematopoietic progenitors", "Erythroblasts", "Monocytes", "Immature B cells"), ordered = TRUE)
# umap_dat.tfidf$panels <- c(rep("real", length(real_celltype)),
#                      rep("scReadSim", length(scReadSim_celltype)),
#                      rep("SCAN-ATAC-Sim", length(SCANATACSim_celltype)))
# umap_dat.tfidf$panels <- factor(umap_dat.tfidf$panels,
#                           levels = c("real", "scReadSim", "SCAN-ATAC-Sim")) 
# save(umap_dat.tfidf, file=sprintf("%s/Fig2_DimensionReduction.UMAPtfidf.dataframe.20230310.Rdata", out_directory))
# # Dataframe for projected UMAP
# umap_dat.projected.tfidf <- data.frame(
#   rbind(Real.UMAP, scReadSim.UMAPtfidf.projected, SCANATACSim.UMAPtfidf.projected))
# colnames(umap_dat.projected.tfidf) <- c('x', 'y')
# umap_dat.projected.tfidf$labels <- c(real_celltype, scReadSim_celltype, SCANATACSim_celltype)
# umap_dat.projected.tfidf$labels <- factor(umap_dat.projected.tfidf$labels, levels=c("Hematopoietic progenitors", "Erythroblasts", "Monocytes", "Immature B cells"), ordered = TRUE)
# umap_dat.projected.tfidf$panels <- c(rep("real", length(real_celltype)),
#                                rep("scReadSim", length(scReadSim_celltype)),
#                                rep("SCAN-ATAC-Sim", length(SCANATACSim_celltype)))
# umap_dat.projected.tfidf$panels <- factor(umap_dat$panels,
#                                     levels = c("real", "scReadSim", "SCAN-ATAC-Sim")) 
# save(umap_dat.projected.tfidf, file=sprintf("%s/Fig2_DimensionReduction.UMAPtfidfprojected.dataframe.20230310.Rdata", out_directory))
# 
# ### plotting
# load(file=sprintf("%s/Fig2_DimensionReduction.UMAPtfidf.dataframe.20230310.Rdata", out_directory))
# 
# point_size <- 0.25; point_transp <- 0.8
# cbPalette <- cell_type_pallete  # wesanderson Royal1
# names(cbPalette) <- levels(umap_dat.tfidf$labels)
# 
# umap_plot.tfidf <- umap_dat.tfidf %>% 
#   # filter( x < 10 & x > -10 ) %>%
#   # filter(y > -10 & y < 10 &  x < 10 & x > -10 ) %>%
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
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 15),
#     legend.text=element_text(size=15),
#     strip.text = element_text(size = 15), 
#     legend.title=element_text(size=15)) +
#   guides(color=guide_legend(override.aes = list(size = 3),
#                             nrow=1, byrow=FALSE,
#                             title="cell type")) +
#   # get_colScale(umap_dat.tfidf$labels) +
#   scale_colour_manual(values = cbPalette) + # wesanderson Royal1
#   xlab("UMAP 1") + ylab("UMAP 2")
# ggsave(file=sprintf("%s/Fig2_DimensionReduction.20230312.pdf",out_directory), umap_plot.tfidf, height = 4.8, width = 9.6, dpi = 600)
# # ggsave(file=sprintf("%s/VerifyCount_Fig4_UMAP.SeuratDR.noguide.pdf",umitool.out.directory), umap_mixed_plot1, height = 4.8, width = 9.6, dpi = 600)



########################### Fig C. Mixed Matrix Dimension reduction ###########################
# PCA + TFIDF
# Preprocess using filtered ID obtained in Figure 2
mixed.scReadSim.processed.matrix <- logtrans(tf_idf(cbind(real.count.mat,scReadSim.count.mat)[peak.filterd.id,]))
mixed.SCANATACSim.processed.matrix <- logtrans(tf_idf(cbind(real.count.mat,SCANATACSim.count.mat)[peak.filterd.id,]))

# PCA: Projected using real projection matrix
mixed.scReadSim.pcatfidf.realprojected <- t(t(mixed.scReadSim.processed.matrix) %*% res_pca$v) # K by n
mixed.SCANATACSim.pcatfidf.realprojected <- t(t(mixed.SCANATACSim.processed.matrix) %*% res_pca$v) # K by n
save(mixed.scReadSim.pcatfidf.realprojected, mixed.SCANATACSim.pcatfidf.realprojected, file=sprintf("%s/Fig3_MixedMatDimensionReduction.RealPCASpace.PCAtfidf.20230313.Rdata", out_directory))

# UMAP
mixed.scReadSim.UMAPtfidf.realprojected <- predict(real_umap_fit, t(mixed.scReadSim.pcatfidf.realprojected))
mixed.SCANATACSim.UMAPtfidf.realprojected <- predict(real_umap_fit, t(mixed.SCANATACSim.pcatfidf.realprojected))
save(mixed.scReadSim.UMAPtfidf.realprojected, mixed.SCANATACSim.UMAPtfidf.realprojected, file=sprintf("%s/Fig3_MixedMatDimensionReduction.RealPCASpace.UMAPtfidf.20230313.Rdata", out_directory))

# Prepare dataframe for plot (TFIDF)
mixed.umap_dat.tfidf.realprojected <- data.frame(
  rbind(mixed.scReadSim.UMAPtfidf.realprojected, mixed.SCANATACSim.UMAPtfidf.realprojected))
colnames(mixed.umap_dat.tfidf.realprojected) <- c('x', 'y')
mixed.umap_dat.tfidf.realprojected$labels <- c(rep("real", ncol(real.count.mat)), rep("scReadSim", ncol(real.count.mat)),
                                               rep("real", ncol(real.count.mat)), rep("SCAN-ATAC-Sim", ncol(real.count.mat)))
mixed.umap_dat.tfidf.realprojected$labels <- factor(mixed.umap_dat.tfidf.realprojected$labels, levels=c("real", "scReadSim", "SCAN-ATAC-Sim"), ordered = TRUE)
mixed.umap_dat.tfidf.realprojected$panels <- c(rep("real + scReadSim", 2*length(real_celltype)),
                                               rep("real + SCAN-ATAC-Sim", 2*length(real_celltype)))
mixed.umap_dat.tfidf.realprojected$panels <- factor(mixed.umap_dat.tfidf.realprojected$panels,
                                                    levels = c("real + scReadSim", "real + SCAN-ATAC-Sim")) 
save(mixed.umap_dat.tfidf.realprojected, file=sprintf("%s/Fig3_MixedMatDimensionReduction.RealPCASpace.UMAPtfidf.dataframe.20230313.Rdata", out_directory))

# Calculate LISI
library(lisi)
get_lisi <- function(x, subset_label){
  emb_dat <- subset(x, subset = panels == subset_label,
                    select = c("x", "y"))
  # str(emb_dat)
  emb_label <- subset(x, subset = panels == subset_label,
                      select = "labels", drop = FALSE)
  # str(emb_label)
  # str(compute_lisi(emb_dat, emb_label, "labels"))
  compute_lisi(emb_dat, emb_label, "labels")
}
subset_label_list <- c("real + scReadSim", "real + SCAN-ATAC-Sim")

mlisi_umap_mixed <- sapply(subset_label_list, function(subset_label){
  median(unlist(get_lisi(mixed.umap_dat.tfidf.realprojected, subset_label)))
})
print(mlisi_umap_mixed)
lisi_panel_vec <- paste0(levels(mixed.umap_dat.tfidf.realprojected$panels), "\n",
                         "[miLISI = ",
                         format(round(mlisi_umap_mixed, 3), nsmall = 3), "]")
mixed.umap_dat.tfidf.realprojected$panels <- mapvalues(mixed.umap_dat.tfidf.realprojected$panels, from = levels(mixed.umap_dat.tfidf.realprojected$panels),
                                                       to = lisi_panel_vec)
save(mixed.umap_dat.tfidf.realprojected, file=sprintf("%s/Fig3_MixedMatDimensionReduction.RealPCASpace.UMAPtfidf.withLISI.dataframe.20230313.Rdata", out_directory))



# Plot
load(file=sprintf("%s/Fig3_MixedMatDimensionReduction.RealPCASpace.UMAPtfidf.dataframe.20230313.Rdata", out_directory))
mixed.cbPalette <- method_pallete # wesanderson Royal2
names(mixed.cbPalette) <- levels(mixed.umap_dat.tfidf.realprojected$labels)

mixed.umap_plot.tfidf <- mixed.umap_dat.tfidf.realprojected %>% 
  # filter( x < 10 & x > -10 ) %>%
  # filter(y > -10 & y < 10 &  x < 10 & x > -10 ) %>%
  ggplot(aes(x = x, y = y, color = labels)) +
  ggrastr::rasterize(geom_point(cex = point_size, alpha = point_transp),dpi = 300) +
  facet_wrap(~panels, nrow = 1) +
  theme_bw() +
  theme(
    # panel.border = element_rect(colour = "black", size=1),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    # legend.title = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    legend.text=element_text(size=15),
    strip.text = element_text(size = 15), 
    legend.title=element_text(size=15)) +
  guides(color=guide_legend(override.aes = list(size = 3),
                            nrow=1, byrow=FALSE,
                            title="method")) +
  # get_colScale(umap_dat.tfidf$labels) +
  scale_colour_manual(values = mixed.cbPalette) + # wesanderson Royal1
  xlab("UMAP 1") + ylab("UMAP 2")
ggsave(file=sprintf("%s/Fig3_MixedMatDimensionReduction.RealPCASpace.noLISI.20230313.pdf",out_directory), mixed.umap_plot.tfidf, height = 4.8, width = 9.6, dpi = 600)

