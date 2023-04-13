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

samplename <- "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.countmatrix"

filename_pre <- sprintf("%s_VerifyCount", samplename)

sample.format <- "txt"
# directory <- "/home/gayan/Projects/scATAC_Simulator/results/20211126_e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1_NONINPUT"
directory <- "/home/gayan/Projects/scATAC_Simulator/results/20220314_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT"

out.directory <- directory
barcode.file <- sprintf("%s/synthetic_cell_barcode.txt",out.directory) 

# Generate synthetic count matrix
## Read in count matrix
cat(sprintf("Reading real count matrix %s.%s...\n", samplename, sample.format))
count_matrix <- read.table(sprintf("%s/%s.%s", directory, samplename, sample.format), sep="\t",header = FALSE)
matrix_num <- count_matrix[,2:ncol(count_matrix)]
cat("Reading synthetic count matrix...\n")
# simu_matrix <- read.table(sprintf("%s/%s.scDesign2Simulated.%s", directory, samplename, sample.format), sep="\t",header = FALSE)
simu_matrix <- read.csv(sprintf("%s/%s.scDesign2Simulated.%s", directory, samplename, sample.format), sep="\t", header = TRUE)


clustering_result <- unlist(read.table(sprintf("%s/%s.LouvainClusterResults.txt", out.directory, samplename), sep="\t",header = FALSE),use.names=FALSE)
peak_list <- count_matrix[,1]
rownames(matrix_num) <- peak_list
rownames(simu_matrix) <- peak_list
colnames(matrix_num) <- clustering_result
# Convert simu count matrix colnames "X1.11" to "1"
names_withdots <- colnames(simu_matrix)[grepl(".", colnames(simu_matrix), fixed=TRUE)]
colnames(simu_matrix)[grepl(".", colnames(simu_matrix), fixed=TRUE)] <- unlist(lapply(names_withdots, function(x){strsplit(x, ".",  fixed=TRUE)[[1]][1]}))
colnames(simu_matrix) <- unlist(lapply(colnames(simu_matrix), function(string){substring(string, 2)}))

n_cell_type <- length(unique(clustering_result))
n_cell_each <- table(clustering_result)


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


########################### Fig A. Summary statistics ###########################
stats_real <- get_stats(matrix_num, 'Real')
stats_simu <- get_stats(simu_matrix, 'scReadSim')
stats_dat <- rbind(stats_real, stats_simu)
stats_dat$group <- factor(stats_dat$group, levels = c('Real', 'scReadSim'))
measures1 <-  c("mean", "var", "cv", "drop_gene",
                "drop_cell", "libsize")
measures2 <-  c("feature mean", "feature variance", "feature cv",
                "feature zero prop.", "cell zero prop.", "cell library size")
stats_dat$measure <- factor(stats_dat$measure, levels = measures1)
stats_dat$measure <- mapvalues(stats_dat$measure, from = measures1, to = measures2)
save(stats_dat, file=sprintf("%s/VerifyCount_Fig1_SummaryStats.Rdata", out.directory))

# create violin-plots to compare the marginal stat values -------------------------------
load(sprintf("%s/VerifyCount_Fig1_SummaryStats.Rdata", out.directory))
stats_dat$group <- mapvalues(stats_dat$group, from = c("Real", "scReadSim"), to = c("real", "scReadSim"))

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
group_list <- c('Real', 'scReadSim')

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
lognormal <- function(Y){
  Y_norm <- Y / colSums(Y) * 1e4
  return(log(Y_norm+1))
}
logtrans <- function(Y){
  return(log(Y+1))
}
pca_fun <- function(Y, pca_dim=30){
  res_pca <- irlba::irlba(t(Y), nv=pca_dim)
  score_pca <- res_pca$u %*% diag(res_pca$d)
  return(t(score_pca))
}

## Filter features
matrix_combined_sub <- filter_nonzero(cbind(matrix_num, simu_matrix),filter_thr=0)


# PCA on subnmatrix
# pca_real_data <- pca_fun(lognormal(matrix_num_sub))
# pca_cluster_label <- pca_fun(lognormal(simu_matrix_sub))
pca_cluster_label_combined <-
  pca_fun(logtrans(matrix_combined_sub))


## UMAP mixed only
umap_result_cluster_label_combined_notprojected <- umap::umap(t(pca_cluster_label_combined))$layout
umap_dat_mixed <- as.data.frame(rbind(umap_result_cluster_label_combined_notprojected,umap_result_cluster_label_combined_notprojected))
colnames(umap_dat_mixed) <- c('x', 'y')
umap_dat_mixed$labels <- c(rep("Real", ncol(matrix_num)), rep("scReadSim", ncol(simu_matrix)),
                           colnames(matrix_num),
                           colnames(simu_matrix))
umap_dat_mixed$labels <- factor(umap_dat_mixed$labels,
                                levels = c(as.character(sort(as.integer(unique(colnames(matrix_num))))),
                                           "scReadSim", "Real"))
umap_dat_mixed$panels <- c(rep("Real + scReadSim",
                               ncol(matrix_num) + ncol(simu_matrix)),
                           rep("Cell Types",
                               ncol(matrix_num) + ncol(simu_matrix)))
umap_dat_mixed$panels <- factor(umap_dat_mixed$panels,
                                levels = c("Real + scReadSim","Cell Types"))
# save(umap_dat_mixed, file=sprintf("%s/VerifyCount_Fig4_UMAP.Rdata", out.directory))
save(umap_dat_mixed, file=sprintf("%s/VerifyCount_Fig4_UMAP.lognonorm.Rdata", out.directory))

load(sprintf("%s/VerifyCount_Fig4_UMAP.lognonorm.Rdata", out.directory))
# which(umap_dat_mixed %>% filter(panels == "Real + scReadSim") %>% select(x) > 10)
# colSums(matrix_combined_sub[, which(umap_dat_mixed %>% filter(panels == "Real + scReadSim") %>% select(x) > 10)])
# colSums(matrix_combined_sub)
# # plotting --------------------------------------------------------------------------------
library(ggplot2); theme_set(theme_bw());
library(RColorBrewer)
library(gridExtra)
library(plyr) # mapvalues
library(cowplot) # ggsave
library(ggpubr) # as_ggplot

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


# UMAP mixed
subset_label_list <- c("Real + scReadSim", "Cell Types")
mlisi_umap_mixed <- sapply(subset_label_list, function(subset_label){
  median(unlist(get_lisi(umap_dat_mixed, subset_label)))
})
print(mlisi_umap_mixed)
lisi_panel_vec <- c(paste0(levels(umap_dat_mixed$panels)[1], "\n",
                           "[miLISI = ",
                           format(round(mlisi_umap_mixed[1], 3), nsmall = 3), "]"), levels(umap_dat_mixed$panels)[2])
umap_dat_mixed$panels <- mapvalues(umap_dat_mixed$panels, from = levels(umap_dat_mixed$panels),
                                   to = lisi_panel_vec)

### plotting
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

save(umap_dat_mixed, file=sprintf("%s/VerifyCount_Fig4_UMAP_withLISI.lognonorm.Rdata", out.directory))
load(file=sprintf("%s/VerifyCount_Fig4_UMAP_withLISI.lognonorm.Rdata", out.directory))
umap_dat_mixedonly <- umap_dat_mixed %>% filter(panels != "Cell Types")
umap_dat_mixedonly$labels <- factor(umap_dat_mixedonly$labels, levels=c("scReadSim", "Real"))
umap_dat_mixedonly$panels <- factor(umap_dat_mixedonly$panels, levels=unique(umap_dat_mixedonly$panels))
save(umap_dat_mixedonly, file=sprintf("%s/VerifyCount_Fig4.MixtureOnly.UMAP_withLISI.20230128.lognonorm.Rdata", out.directory))

load(file=sprintf("%s/VerifyCount_Fig4.MixtureOnly.UMAP_withLISI.20230128.lognonorm.Rdata", out.directory))
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
umap_mixed_plot1 <- umap_dat_mixedonly[-which(umap_dat_mixed$x > 10),] %>%
  ggplot(aes(x = x, y = y, color = labels)) +
  ggrastr::rasterize(geom_point(cex = point_size, alpha = point_transp),dpi = 300) +
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
  guides(color=guide_legend(override.aes = list(size = 1),
                            ncol=1, byrow=FALSE)) +
  scale_colour_manual(values=c('#FF0000', '#000000')) + 
  xlab("UMAP 1") + ylab("UMAP 2")


# No truncation
# Save the dataframe with lisi values to plot with ATAC-seq
# umap_mixed_plot2 <- umap_dat_mixed %>% 
#   ggplot(aes(x = x, y = y, color = labels)) +
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
#                             ncol=2, byrow=FALSE)) +
#   get_colScale(umap_dat_mixed$labels) +
#   xlab("UMAP 1") + ylab("UMAP 2")

# pdf(sprintf("%s/%s.VERIFYCOUNT.pdf",out.directory, samplename),width=5000, height=5000, res=500)
gs <- list(arrangeGrob(stats_plot, left = textGrob("A", x = unit(1, "npc"), 
                                                   y = unit(.95, "npc"))),
           # arrangeGrob(p_featurecorrelation, left = textGrob("B", x = unit(1, "npc"), 
                                                             # y = unit(.95, "npc"))),
           arrangeGrob(cor_tau_plot, left = textGrob("B", x = unit(1, "npc"), 
                                                     y = unit(.95, "npc"))),
           arrangeGrob(umap_mixed_plot1, left = textGrob("C", x = unit(1, "npc"), 
                                                         y = unit(.95, "npc"))))
lay <- rbind(c(1,1,1,1),
             c(1,1,1,1),
             c(1,1,1,1),
             c(2,2,2,2),
             c(2,2,2,2),
             c(4,4,4,4),
             c(4,4,4,4),
             c(4,4,4,4))
p_final <-arrangeGrob(grobs = gs, layout_matrix = lay)
# ggsave(file=sprintf("%s/%s.VERIFYCOUNT.20220422.noTruncation.pdf",out.directory, filename_pre), p_final, height = 9.6, width = 9.6, dpi = 600)
ggsave(file=sprintf("%s/%s.VERIFYCOUNT.20230128.zoomin.pdf",out.directory, filename_pre), p_final, height = 9.6, width = 9.6, dpi = 600)


