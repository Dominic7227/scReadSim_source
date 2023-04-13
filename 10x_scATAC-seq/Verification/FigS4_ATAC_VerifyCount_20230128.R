library(Matrix)
library(Rsubread)
library(pscl)
library(parallel)
library(MASS)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(hrbrthemes)
library(SnapATAC)
library(viridis)
library(reshape2)  # melt
library(grid)
library(gridExtra)
library(plyr)

samplename <- "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.countmatrix"
filename_pre <- sprintf("%s_VerifyCount", samplename)

sample.format <- "txt"
# directory <- "/home/gayan/Projects/scATAC_Simulator/results/20211126_e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1_NONINPUT"
directory <- "/home/gayan/Projects/scATAC_Simulator/results/20220314_e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1_NONINPUT"
# out.directory <- "/home/gayan/Projects/scATAC_Simulator/results/20220407_e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1_NONINPUT"

out.directory <- directory
barcode.file <- sprintf("%s/synthetic_cell_barcode.txt",out.directory) 
peak_file <- sprintf("%s/e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.MACS3_peaks.bed",out.directory) 

# Generate synthetic count matrix
## Read in count matrix
cat(sprintf("Reading real count matrix %s.%s...\n", samplename, sample.format))
count_matrix <- read.table(sprintf("%s/%s.%s", directory, samplename, sample.format), sep="\t",header = FALSE)
matrix_num <- count_matrix[,2:ncol(count_matrix)]
cat("Reading synthetic count matrix...\n")
simu_matrix <- read.table(sprintf("%s/%s.scDesign2Simulated.%s", directory, samplename, sample.format), sep="\t",header = TRUE)

clustering_result <- unlist(read.table(sprintf("%s/%s.LouvainClusterResults.txt", directory, samplename), sep="\t",header = FALSE),use.names=FALSE)
peak_list <- count_matrix[,1]
rownames(matrix_num) <- peak_list
rownames(simu_matrix) <- peak_list
n_cell_type <- length(unique(clustering_result))
n_cell_each <- table(clustering_result)
tmp <- unlist(lapply(strsplit(colnames(simu_matrix),".", fixed=TRUE), function(x){
  if (length(x) == 2){
    return(x[1])
  } else {
    return(x)
  }
}))
colnames(simu_matrix) <- unlist(lapply(tmp, function(x) {return(substring(x, 2))} ), use.names = FALSE)
# colnames(simu_matrix) <- unlist(lapply(1:n_cell_type, function(x){rep(x, n_cell_each[x])}))
colnames(matrix_num) <- clustering_result



# ####################################################################################
# ################################## Summary statistics ###############################
# ####################################################################################

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


# sim.feature.depth.df <- data.frame(group = "scReadSim", location = 1:nrow(simu_matrix),value = sim.feature.depth)
# real.feature.depth.df <- data.frame(group = "Real",location = 1:nrow(matrix_num), value = real.feature.depth)

# df <- rbind(sim.feature.depth.df, real.feature.depth.df)
# p_featuredepth <-ggplot(df, aes(x=location, y=value, group=group, color=group)) +
# geom_line() +
# scale_color_manual(values=c("#F8766D", "#00BFC4")) +
# theme(legend.position="none") +
# # ggtitle("count per feature") +
# theme(
#   legend.position="none",
#   panel.spacing = unit(0.1, "lines"),
#   strip.text.x = element_text(size = 8),
#   plot.title = element_text(size=14)
# ) +
# facet_wrap(~group,ncol=1)
########################### Fig A. Summary statistics ###########################
stats_real <- get_stats(matrix_num, 'Real')
stats_simu <- get_stats(simu_matrix, 'scReadSim')
stats_dat <- rbind(stats_real, stats_simu)
stats_dat$group <- factor(stats_dat$group, levels = c('Real', 'scReadSim'))
measures1 <-  c("mean", "var", "cv", "drop_gene",
                "drop_cell", "libsize")
measures2 <-  c("feature mean", "feature variance", "feature ßcv",
                "feature zero prop.", "cell zero prop.", "cell library size")
stats_dat$measure <- factor(stats_dat$measure, levels = measures1)
stats_dat$measure <- mapvalues(stats_dat$measure, from = measures1, to = measures2)
save(stats_dat, file=sprintf("%s/VerifyCount_Fig1_SummaryStats.Rdata", out.directory))

load(file=sprintf("%s/VerifyCount_Fig1_SummaryStats.Rdata", out.directory))
stats_dat$group <- mapvalues(stats_dat$group, from = c("Real", "scReadSim"), to = c("real", "scReadSim"))

# create violin-plots to compare the marginal stat values -------------------------------
stats_plot <- ggplot(stats_dat, aes(x = group, y = value, fill=group)) +
  geom_violin(scale = 'width', trim = TRUE) +
  scale_fill_manual(values=c("#F8766D", "#00BFC4")) +
  facet_wrap(~measure, scales = "free", ncol = 3) +
  theme_bw() +
  theme(strip.text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x=element_blank(),
        legend.text=element_text(size=12)) +
  xlab("") + ylab("")

########################### Fig B. RPKM ###########################
setwd(out.directory) # outdirectory 20221130_e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1_NONINPUT

ref_feature_set_dir <- "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.MACS3_peaks.bed"
bg_ref_feature_set_dir <- "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.MACS3_thr10_peaks.COMPLE.bed"
ref_feature_set <- read.table(ref_feature_set_dir)
bg_ref_feature_set <- read.table(bg_ref_feature_set_dir)

# Read in real read coverage over real peaks and non-peaks
count_pergene_vec <- read.table(sprintf("%s/e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.countmatrix.MarginalFeatureCount.txt", out.directory))
count_pergene_vec <- as.numeric(count_pergene_vec[,-1])
bg_count_pergene_vec <- read.table(sprintf("%s/e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.COMPLE.countmatrix.MarginalFeatureCount.txt", out.directory))
bg_count_pergene_vec <- as.numeric(bg_count_pergene_vec[,-1])

# Read in MACS3 synthetic read coverage over real peaks and non-peaks
simu_count_MACS3peak <- read.table(sprintf("%s/e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.countmatrix.scDesign2Simulated.MarginalFeatureCount.txt", out.directory))
simu_count_MACS3peak <- as.numeric(simu_count_MACS3peak[,-1])
bg_simu_count_MACS3peak <- read.table(sprintf("%s/e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.COMPLE.countmatrix.scDesign2Simulated.MarginalFeatureCount.txt", out.directory))
bg_simu_count_MACS3peak <- as.numeric(bg_simu_count_MACS3peak[,-1])

# Calculate RPKM
trun_95 <- function(vec){
  return(vec[vec < quantile(vec, 0.95)])
}

# Calculate read coverage density of reference feature for real data
read_density_real <- 1000 * (1000000 / sum(count_pergene_vec)) * count_pergene_vec / (ref_feature_set[,3] - ref_feature_set[,2])
bg_read_density_real <- 1000 * (1000000 / sum(bg_count_pergene_vec)) * bg_count_pergene_vec / (bg_ref_feature_set[,3] - bg_ref_feature_set[,2])
read_density_real_trun <- trun_95(read_density_real)
bg_read_density_real_trun <- trun_95(bg_read_density_real)

# Calculate read coverage density of MACS3 synthetic feature for synthetic data
read_density_MACS3 <- 1000 * (1000000 / sum(simu_count_MACS3peak)) *  simu_count_MACS3peak / (ref_feature_set[,3] - ref_feature_set[,2])
bg_read_density_MACS3 <- 1000 * (1000000 / sum(bg_simu_count_MACS3peak)) *  bg_simu_count_MACS3peak / (bg_ref_feature_set[,3] - bg_ref_feature_set[,2])
read_density_MACS3_trun <- trun_95(read_density_MACS3)
bg_read_density_MACS3_trun <- trun_95(bg_read_density_MACS3)

create_dataftame <- function(inform_RPKM, bg_RPKM, group_name){
  df_tmp <- data.frame(read_density=c(inform_RPKM, bg_RPKM), 
                       group=factor(rep(group_name, length(inform_RPKM) + length(bg_RPKM))),
                       category=factor(c(rep("Foreground Features", length(inform_RPKM)), rep("Background Features", length(bg_RPKM))), levels = c("Foreground Features", "Background Features"), ordered=TRUE))
  return(df_tmp)
}

df_INPUT_Real <- create_dataftame(read_density_real_trun, bg_read_density_real_trun, "Real")
df_INPUT_MACS3 <- create_dataftame(read_density_MACS3_trun, bg_read_density_MACS3_trun, "scReadSim")
list_trun <- list(df_INPUT_Real, df_INPUT_MACS3)
df_trun <- Reduce(rbind, list_trun)
save(df_trun, file=sprintf("%s/RPKM_plot.Rdata", fig.directory))

load(("/home/gayan/Projects/scATAC_Simulator/results/20221130_e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1_NONINPUT/VerifyReadsPlots/RPKM_plot.Rdata"))
measures1 <-  c("Foreground Features", "Background Features")
measures2 <-  c("peak", "non-peak")
df_trun$category <- mapvalues(df_trun$category, from = measures1, to = measures2)
df_trun$group <- mapvalues(df_trun$group, from = c("Real", "scReadSim"), to = c("real", "scReadSim"))

read_des_p <- df_trun %>%
  ggplot( aes(x=category, y=read_density, fill=category, color=category)) + 
  stat_boxplot(geom = "errorbar",width=0.3,size=1)+
  geom_boxplot(size=1, outlier.size=0.3) +
  scale_fill_manual(values=c("#D9D9D9", "#D9D9D9")) +
  scale_color_manual(values=c("black", "black")) +
  theme_bw() +
  theme(legend.position="none",
        legend.title=element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        strip.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("") + ylab("RPKM") +
  facet_wrap(~group, scales = "free_x", ncol=3)



########################### Fig C. Feature correlation ###########################
gene_mean <- apply(matrix_num, 1, mean)
cutoff <- 100
gene_sel <- order(gene_mean, decreasing = TRUE)[1:cutoff]
get_cor_mat <- function(x, cor_fun){
  sub_mat <- x[gene_sel, ]
  cor_fun(t(sub_mat))
}
get_heatmap_dat <- function(mat_list, cor_fun){
  cor_mat_list <- lapply(mat_list, get_cor_mat, cor_fun)
  # reorder cor_mat entries according to hierarchical clustering result
  cor_mat_list <- lapply(cor_mat_list, function(x){
    x[hclust_result$order, hclust_result$order]})
  # organize the cor values as input for ggplot2
  cor_melted <- lapply(cor_mat_list, melt)
  cor_dat <- Reduce(rbind, cor_melted)
  cor_dat$group <- unlist(lapply(1:length(group_list), function(x){
    rep(group_list[[x]], nrow(cor_melted[[x]]))
  }))
  return(cor_dat)
}

# calculate the correlations and organize as input for ggplot2 --------------------------
rownames(simu_matrix) <- rownames(matrix_num)
mat_list <- list(real = matrix_num,  scReadSim = simu_matrix)
hclust_result <- hclust(as.dist(1-get_cor_mat(mat_list$real, cor)))

cor_dat <- get_heatmap_dat(mat_list, cor)
cor_tau_dat <- cor_dat
cor_tau_dat$group <- factor(cor_tau_dat$group, levels = group_list)


save(cor_tau_dat, file=sprintf("%s/VerifyCount_Fig3_FeatureCorrelation.Rdata", out.directory))

load(file=sprintf("%s/VerifyCount_Fig3_FeatureCorrelation.Rdata", out.directory))
group_list <- c('Real', 'scReadSim')
cor_tau_dat$group <- mapvalues(cor_tau_dat$group, from = group_list, to = c('real', 'scReadSim'))

# create heatmaps to display the correlation values -------------------------------------
cor_tau_plot <- ggplot(cor_tau_dat, aes(Var2, Var1, fill = value))+
  facet_wrap(~group, ncol = 2) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1,1), space = "Lab",
                       name="Pearson") +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_text(size=15),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) +
  xlab("") + ylab("") + coord_fixed()


########################### Fig D. Dimension reduction ###########################
# # use t-SNE to visualize the simulated data based on cluster label ------------------------
# library(Rtsne)
filter_nonzero <- function(data_mat, filter_thr=0){
  nonzero_prop <- apply(data_mat, 1, function(x){
    sum(x>0)
  }) / ncol(data_mat)
  data_mat_filtered <- data_mat[nonzero_prop > filter_thr,]
  return(data_mat_filtered)
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

pca_fun <- function(Y, pca_dim=30){
  res_pca <- irlba::irlba(t(Y), nv=pca_dim)
  score_pca <- res_pca$u %*% diag(res_pca$d)
  return(t(score_pca))
}
# tsne_fun <- function(mat, rand_seed = 2022){
#   set.seed(rand_seed)
#   pca_tsne <- Rtsne::Rtsne(t(mat),check_duplicates = FALSE)$Y
#   pca_tsne <- as.data.frame(pca_tsne)
#   colnames(pca_tsne) <- c('tsne1','tsne2')
#   return(pca_tsne)
# }

## Try Seurat DE peaks identification
# library(Signac)
# library(Seurat)

# rownames_new <- unlist(lapply(rownames(simu_matrix), function(x){
#   paste0(unlist(strsplit(x, split = "_")), collapse="-")
#   }))
# rownames(simu_matrix) <- rownames_new
# rownames(matrix_num) <- rownames_new
# brain_assay_combined <- CreateChromatinAssay(
#   counts = matrix_num
# )
# brain <- CreateSeuratObject(
#   counts = brain_assay_combined,
#   assay = 'peaks',
#   project = 'ATAC'
# )
# da_peaks <- FindMarkers(
#   object = brain,
#   ident.1 = c("1"), 
#   ident.2 = c("2", "3", "4"),
#   min.pct = 0.05,
#   test.use = 'LR',
#   latent.vars = 'peak_region_fragments'
# )

## Filter features
# matrix_num_sub <- filter_nonzero(matrix_num,filter_thr=0.1)
# simu_matrix_sub <- filter_nonzero(simu_matrix,filter_thr=0.1)
matrix_combined_sub <- filter_nonzero(cbind(matrix_num, simu_matrix),filter_thr=0.01)
# count_features <- rowSums(matrix_num)
# filter_level <- 0.25
# filter_thr <- quantile(count_features, filter_level)
# filtered_ind <- which(count_features >= filter_thr)
# matrix_num_sub <- matrix_num[filtered_ind,]
# simu_matrix_sub <- simu_matrix[filtered_ind,]
# matrix_num_sub <- matrix_num_sub[,which(colSums(matrix_num[filtered_ind,]) > 0)] # remove zero cells
# simu_matrix_sub <- simu_matrix_sub[,which(colSums(matrix_num[filtered_ind,]) > 0)]
# clustering_result_sub <- clustering_result[which(colSums(matrix_num[filtered_ind,]) > 0)]


# PCA
# pca_real_data <- pca_fun(tf_idf(matrix_num_sub))
# pca_cluster_label <- pca_fun(tf_idf(simu_matrix_sub))
pca_cluster_label_combined <-
  pca_fun(logtrans(tf_idf(matrix_combined_sub)))

# # tsne on pca of subnmatrix
# tsne_result_real_data <- tsne_fun(pca_real_data)
# tsne_result_cluster_label <- tsne_fun(pca_cluster_label)
# tsne_result_cluster_label_combined <-
#    tsne_fun(pca_cluster_label_combined)
# save(pca_real_data, pca_cluster_label, pca_cluster_label_combined, tsne_result_real_data, tsne_result_cluster_label, tsne_result_cluster_label_combined, file=sprintf("%s/%s.selectedCellType.filtered.tf_idf_pca_tsne.Rdata",out.directory, filename_pre))

## UMAP
# real_umap_fit <- umap::umap(t(pca_real_data))
# umap_result_real_data <- real_umap_fit$layout
# umap_result_cluster_label <- predict(real_umap_fit, t(pca_cluster_label))
# umap_result_cluster_label_combined <-
#    predict(real_umap_fit, t(pca_cluster_label_combined))
# save(umap_result_real_data, umap_result_cluster_label, umap_result_cluster_label_combined, file=sprintf("%s/%s.selectedCellType.filtered.tfidf_pca_umap.Rdata",out.directory, filename_pre))

## UMAP mixed only
umap_result_cluster_label_combined_notprojected <- umap::umap(t(pca_cluster_label_combined))$layout
umap_dat_mixed <- as.data.frame(rbind(umap_result_cluster_label_combined_notprojected,umap_result_cluster_label_combined_notprojected))
colnames(umap_dat_mixed) <- c('x', 'y')
umap_dat_mixed$labels <- c(rep("Real", ncol(matrix_num)), rep("scReadSim", ncol(simu_matrix)),
                           colnames(matrix_num),
                           colnames(simu_matrix))
umap_dat_mixed$labels <- factor(umap_dat_mixed$labels,
                                levels = c(sort(unique(colnames(matrix_num))),
                                           "scReadSim", "Real"))
umap_dat_mixed$panels <- c(rep("Real + scReadSim",
                               ncol(matrix_num) + ncol(simu_matrix)),
                           rep("Cell Types",
                               ncol(matrix_num) + ncol(simu_matrix)))
umap_dat_mixed$panels <- factor(umap_dat_mixed$panels,
                                levels = c("Real + scReadSim","Cell Types"))
# save(umap_dat_mixed, file=sprintf("%s/VerifyCount_Fig4_UMAP.20220422.1percent.Rdata", out.directory))
save(umap_dat_mixed, file=sprintf("%s/VerifyCount_Fig4_UMAP.20220422.1percent.TFIDFlog.Rdata", out.directory))

load(sprintf("%s/VerifyCount_Fig4_UMAP.20220422.1percent.TFIDFlog.Rdata", out.directory))

# which(umap_dat_mixed %>% filter(panels == "Real + scReadSim") %>% select(y) > 10)
# colSums(matrix_combined_sub[, which(umap_dat_mixed %>% filter(panels == "Real + scReadSim") %>% select(y) > 10)])
# colSums(matrix_combined_sub)

# # plotting --------------------------------------------------------------------------------
library(ggplot2); theme_set(theme_bw());
library(RColorBrewer)
library(gridExtra)
library(plyr) # mapvalues
library(cowplot) # ggsave
library(ggpubr) # as_ggplot


### dat
# pca_dat <- data.frame(
#   rbind(t(pca_real_data[1:2,]), t(pca_cluster_label[1:2,]),
#         t(pca_cluster_label_combined[1:2,])))
# colnames(pca_dat) <- c('x', 'y')
# pca_dat$labels <- c(colnames(matrix_num_sub),
#                        colnames(simu_matrix_sub),
#                        rep("scReadSim", ncol(simu_matrix_sub)),
#                        rep("Real", ncol(matrix_num_sub)))
# pca_dat$labels <- factor(pca_dat$labels,
#                              levels = c(sort(unique(colnames(matrix_num_sub))),
#                                         "scReadSim", "Real"))
# pca_dat$panels <- c(rep("Real", ncol(matrix_num_sub)),
#                         rep("scReadSim", ncol(simu_matrix_sub)),
#                         rep("Real + scReadSim",
#                             ncol(matrix_num_sub) + ncol(simu_matrix_sub)))
# pca_dat$panels <- factor(pca_dat$panels,
#                              levels = c("Real", "scReadSim",
#                                         "Real + scReadSim"))

# dim_red_dat <- data.frame(
#   rbind(tsne_result_real_data, tsne_result_cluster_label,
#         tsne_result_cluster_label_combined))
# colnames(dim_red_dat) <- c('x', 'y')
# dim_red_dat$labels <- c(colnames(matrix_num_sub),
#                        colnames(simu_matrix_sub),
#                        rep("scReadSim", ncol(simu_matrix_sub)),
#                        rep("Real", ncol(matrix_num_sub)))
# dim_red_dat$labels <- factor(dim_red_dat$labels,
#                              levels = c(sort(unique(colnames(matrix_num_sub))),
#                                         "scReadSim", "Real"))
# dim_red_dat$panels <- c(rep("Real", ncol(matrix_num_sub)),
#                         rep("scReadSim", ncol(simu_matrix_sub)),
#                         rep("Real + scReadSim",
#                             ncol(matrix_num_sub) + ncol(simu_matrix_sub)))
# dim_red_dat$panels <- factor(dim_red_dat$panels,
#                              levels = c("Real", "scReadSim",
#                                         "Real + scReadSim"))


# umap_dat <- data.frame(
#   rbind(umap_result_real_data, umap_result_cluster_label,
#         umap_result_cluster_label_combined))
# colnames(umap_dat) <- c('x', 'y')
# umap_dat$labels <- c(colnames(matrix_num_sub),
#                        colnames(simu_matrix_sub),
#                        rep("scReadSim", ncol(simu_matrix_sub)),
#                        rep("Real", ncol(matrix_num_sub)))
# umap_dat$labels <- factor(umap_dat$labels,
#                              levels = c(sort(unique(colnames(matrix_num_sub))),
#                                         "scReadSim", "Real"))
# umap_dat$panels <- c(rep("Real", ncol(matrix_num_sub)),
#                         rep("scReadSim", ncol(simu_matrix_sub)),
#                         rep("Real + scReadSim",
#                             ncol(matrix_num_sub) + ncol(simu_matrix_sub)))
# umap_dat$panels <- factor(umap_dat$panels,
#                              levels = c("Real", "scReadSim",
#                                         "Real + scReadSim"))

# ### lisi
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


# subset_label_list <- c("Real", "scReadSim", "Real + scReadSim")
# PCA
# mlisi_pca <- sapply(subset_label_list, function(subset_label){
#   median(unlist(get_lisi(pca_dat, subset_label)))
# })
# print(mlisi_pca)
# pca_dat$panels <- mapvalues(pca_dat$panels, from = levels(pca_dat$panels),
#                                 to = c(paste0(levels(pca_dat$panels), "\n",
#                                               c("[mcLISI = ", "[mcLISI = ", "[miLISI = "),
#                                               format(round(mlisi_pca, 3), nsmall = 3), "]")))

# mlisi_tsne <- sapply(subset_label_list, function(subset_label){
#   median(unlist(get_lisi(dim_red_dat, subset_label)))
# })
# print(mlisi_tsne)

# tSNE
# dim_red_dat$panels <- mapvalues(dim_red_dat$panels, from = levels(dim_red_dat$panels),
#                                 to = c(paste0(levels(dim_red_dat$panels), "\n",
#                                               c("[mcLISI = ", "[mcLISI = ", "[miLISI = "),
#                                               format(round(mlisi_tsne, 3), nsmall = 3), "]")))

# UMAP
# mlisi_umap <- sapply(subset_label_list, function(subset_label){
#   median(unlist(get_lisi(umap_dat, subset_label)))
# })
# print(mlisi_umap)
# umap_dat$panels <- mapvalues(umap_dat$panels, from = levels(umap_dat$panels),
#                                 to = c(paste0(levels(umap_dat$panels), "\n",
#                                               c("[mcLISI = ", "[mcLISI = ", "[miLISI = "),
#                                               format(round(mlisi_umap, 3), nsmall = 3), "]")))

# UMAP mixed
subset_label_list <- c("Real + scReadSim", "Cell Types")
mlisi_umap_mixed <- sapply(subset_label_list, function(subset_label){
  median(unlist(get_lisi(umap_dat_mixed, subset_label)))
})
print(mlisi_umap_mixed)
# umap_dat_mixed$panels <- mapvalues(umap_dat_mixed$panels, from = levels(umap_dat_mixed$panels),
#                                 to = c(paste0(levels(umap_dat_mixed$panels), "\n",
#                                               c("[miLISI = ", "[mcLISI = "),
#                                               format(round(mlisi_umap_mixed, 3), nsmall = 3), "]")))
lisi_panel_vec <- c(paste0(levels(umap_dat_mixed$panels)[1], "\n",
                           "[miLISI = ",
                           format(round(mlisi_umap_mixed[1], 3), nsmall = 3), "]"), levels(umap_dat_mixed$panels)[2])
umap_dat_mixed$panels <- mapvalues(umap_dat_mixed$panels, from = levels(umap_dat_mixed$panels),
                                   to = lisi_panel_vec)

### plotting


# pca_red_plot1 <- ggplot(pca_dat, aes(x = x, y = y, color = labels)) +
#   geom_point(cex = point_size, alpha = point_transp) +
#   facet_wrap(~panels, nrow = 1) +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "bottom",
#         legend.title = element_blank()) +
#   guides(color=guide_legend(override.aes = list(size = 1),
#                             nrow=2, byrow=FALSE)) +
#   get_colScale(pca_dat$labels) +
#   xlab("PC 1") + ylab("PC 2")

# dim_red_plot1 <- ggplot(dim_red_dat, aes(x = x, y = y, color = labels)) +
#   geom_point(cex = point_size, alpha = point_transp) +
#   facet_wrap(~panels, nrow = 1) +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "bottom",
#         legend.title = element_blank()) +
#   guides(color=guide_legend(override.aes = list(size = 1),
#                             nrow=2, byrow=FALSE)) +
#   get_colScale(dim_red_dat$labels) +
#   xlab("TSNE 1") + ylab("TSNE 2")

# umap_red_plot1 <- ggplot(umap_dat, aes(x = x, y = y, color = labels)) +
# geom_point(cex = point_size, alpha = point_transp) +
# facet_wrap(~panels, nrow = 1) +
# theme(plot.title = element_text(hjust = 0.5),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.position = "bottom",
#       legend.title = element_blank()) +
# guides(color=guide_legend(override.aes = list(size = 1),
#                           nrow=2, byrow=FALSE)) +
# get_colScale(umap_dat$labels) +
# xlab("UMAP 1") + ylab("UMAP 2")

save(umap_dat_mixed, file=sprintf("%s/VerifyCount_Fig4_UMAP_withLISI.20220515.1percent.TFIDFlog.Rdata", out.directory))
load(file=sprintf("%s/VerifyCount_Fig4_UMAP_withLISI.20220515.1percent.TFIDFlog.Rdata", out.directory))

umap_dat_mixedonly <- umap_dat_mixed %>% filter(panels != "Cell Types")
umap_dat_mixedonly$labels <- factor(umap_dat_mixedonly$labels, levels=c("scReadSim", "Real"))
umap_dat_mixedonly$panels <- factor(umap_dat_mixedonly$panels, levels=unique(umap_dat_mixedonly$panels))
save(umap_dat_mixedonly, file=sprintf("%s/VerifyCount_Fig4.MixtureOnly.UMAP_withLISI.20221015.1percent.TFIDFlog.Rdata.Rdata", out.directory))

load(file=sprintf("%s/VerifyCount_Fig4.MixtureOnly.UMAP_withLISI.20221015.1percent.TFIDFlog.Rdata.Rdata", out.directory))
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
umap_mixed_plot1 <- umap_dat_mixedonly %>% filter(y>-10 & y <10) %>%
  ggplot(aes(x = x, y = y, color = labels)) +
  ggrastr::rasterize(geom_point(cex = point_size, alpha = point_transp), dpi = 300) +
  facet_wrap(~panels, nrow = 1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        strip.text = element_text(size = 15)) +
  guides(color=guide_legend(override.aes = list(size = 3, alpha = 1),
                            ncol=1, byrow=FALSE)) +
  scale_colour_manual(values=c('#FF0000', '#000000')) + 
  xlab("UMAP 1") + ylab("UMAP 2")

# # With truncation
#  umap_mixed_plot1 <- umap_dat_mixed[-which(umap_dat_mixed$y > 10),] %>% 
#  ggplot(aes(x = x, y = y, color = labels)) +
#    ggrastr::rasterize(geom_point(cex = point_size, alpha = point_transp),dpi = 300) +
#    facet_wrap(~panels, nrow = 1) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "right",
#         legend.title = element_blank(), 
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 15),
#         legend.text=element_text(size=12),
#         strip.text = element_text(size = 15)) +
#   guides(color=guide_legend(override.aes = list(size = 1),
#                             ncol=1, byrow=FALSE)) +
#   get_colScale(umap_dat_mixed$labels) +
#   xlab("UMAP 1") + ylab("UMAP 2")
# 
# # No truncation
# # Save the dataframe with lisi values to plot with RNA-seq
# umap_mixed_plot2 <- ggplot(umap_dat_mixed, aes(x = x, y = y, color = labels)) +
#   ggrastr::rasterize(geom_point(cex = point_size, alpha = point_transp),dpi = 300) +
#   facet_wrap(~panels, nrow = 1) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = "right",
#         legend.title = element_blank(), 
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 15),
#         legend.text=element_text(size=12),
#         strip.text = element_text(size = 15)) +
#   guides(color=guide_legend(override.aes = list(size = 1),
#                             ncol=1, byrow=FALSE)) +
#   get_colScale(umap_dat_mixed$labels) +
#   xlab("UMAP 1") + ylab("UMAP 2")


# pdf(sprintf("%s/%s.VERIFYCOUNT.pdf",out.directory, samplename),width=5000, height=5000, res=500)
gs <- list(arrangeGrob(stats_plot, left = textGrob("A", x = unit(1, "npc"), 
                                                   y = unit(.95, "npc"))),
           arrangeGrob(read_des_p, left = textGrob("B", x = unit(1, "npc"),
                                                   y = unit(.95, "npc"))),
           arrangeGrob(cor_tau_plot, left = textGrob("B", x = unit(1, "npc"), 
                                                     y = unit(.95, "npc"))),
           arrangeGrob(umap_mixed_plot1, left = textGrob("C", x = unit(1, "npc"), 
                                                         y = unit(.95, "npc"))))
lay <- rbind(c(1,1,1,1,1),
             c(1,1,1,1,1),
             c(1,1,1,1,1),
             c(2,2,3,3,3),
             c(2,2,3,3,3),
             c(4,4,4,4,4),
             c(4,4,4,4,4),
             c(4,4,4,4,4))
p_final <-arrangeGrob(grobs = gs, layout_matrix = lay)
ggsave(file=sprintf("%s/%s.VERIFYCOUNT.20230128.zoomin.pdf",out.directory, filename_pre), p_final, height = 9.6, width = 9.6, dpi = 600)

# 
# png(paste(out.directory, "/", filename_pre, "_feature_depth.png", sep = ""))
# print(p_featuredepth+xlab("feature id")+ylab("log count"))
# dev.off()
# png(paste(out.directory, "/", filename_pre, "_stats_violinplot.png", sep = ""), width=5000, height=5000, res=500)
# print(stats_plot)
# dev.off()
# png(paste(out.directory, "/", filename_pre, "_featurecorrelation_heatmap.png", sep = ""), width=5000, height=5000, res=500)
# print(cor_tau_plot)
# dev.off()
# png(sprintf("%s/%s.real_synthetic_pca_comparison.png",out.directory, filename_pre),width=5000, height=5000, res=500)
# print(pca_red_plot1)
# dev.off()
# # png(sprintf("%s/%s.real_synthetic_pcatsne_comparison.png",out.directory, filename_pre))
# # print(dim_red_plot1)
# # dev.off()
# png(sprintf("%s/%s.real_synthetic_pcaumap_comparison.png",out.directory, filename_pre))
# print(umap_red_plot1)
# dev.off()
# png(sprintf("%s/%s.real_synthetic_pcaumap_mixed.png",out.directory, filename_pre))
# print(umap_mixed_plot1)
# dev.off()

