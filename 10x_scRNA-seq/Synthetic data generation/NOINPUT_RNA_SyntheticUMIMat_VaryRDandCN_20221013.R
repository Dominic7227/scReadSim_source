library(Matrix)
library(Rsubread)
library(pscl)
library(parallel)
library(MASS)
# library(ggpubr)
# library(ggplot2)
# library(hrbrthemes)
library(tidyverse)
library(ROGUE)
library(Seurat)

# Calculate marginal count for each region
# Arguments
## samplename: Count matrix file name (without format)
## sample.format: Count matrix file format, like csv
## directory: Count matrix directory
## out.directory: output directory
# args <- commandArgs(trailingOnly = TRUE) # Extract all arguments from Linux input

# samplename <- args[1]
# sample.format <- args[2]
# directory <- args[3]
# out.directory <- args[4]

# samplename <- "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.countmatrix.scDesign2functions"
# sample.format <- "txt"
# directory <- "/home/gayan/Projects/scATAC_Simulator/results/20211126_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT"
# out.directory <- "/home/gayan/Projects/scATAC_Simulator/results/20220116_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT_withCluster"
# dir.create(out.directory)

# samplename <- "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.countmatrix"
# sample.format <- "txt"
# directory <- "/home/gayan/Projects/scATAC_Simulator/results/20211126_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT"
# out.directory <- "/home/gayan/Projects/scATAC_Simulator/results/20211126_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT"


## Clustering
# function of getting a submatrix of selected cell types ----------------------------------
get_submat <- function(data_mat, cell_type_sel){
  if(is.null(colnames(data_mat))){
    data_mat
  }else{
    data_mat[, colnames(data_mat) %in% cell_type_sel]
  }
}
# function of using ROGUE scores to check cluster qualities -------------------------------
check_cluster_quality <- function(data_mat, cell_type_sel,
                                  platform = c("full-length", "UMI")){
  platform <- match.arg(platform)
  
  data_mat_sel <- get_submat(data_mat, cell_type_sel)
  
  expr <- data_mat_sel
  rogue.res <- rogue(expr, labels = colnames(expr),
                     samples = rep(1, ncol(expr)), platform = platform)
  unlist(rogue.res)
}
# function of using Seurat to cluster -----------------------------------------------------
get_cluster_seurat <- function(data_mat, platform = c("full-length", "UMI"),
                               dims = 1:10, res = 0.5){
  platform <- match.arg(platform)
  submat <- data_mat
  
  if(is.null(rownames(submat)))
    rownames(submat) <- 1:nrow(submat)
  count_seurat <- CreateSeuratObject(counts = submat)
  ### normalization
  count_seurat <- NormalizeData(count_seurat)
  ### select highly variable genes
  count_seurat <- FindVariableFeatures(count_seurat, selection.method = "vst", nfeatures = 4000)
  ### scale the data
  count_seurat <- ScaleData(count_seurat)
  ### PCA
  count_seurat <- RunPCA(count_seurat,
                         features = VariableFeatures(object = count_seurat),
                         verbose = F)
  ### clustering
  count_seurat <- FindNeighbors(count_seurat, dims = dims)
  count_seurat <- FindClusters(count_seurat, resolution = res)
  ### results
  cluster_predicted <- as.integer(Idents(count_seurat))
  
  colnames(submat) <- as.character(cluster_predicted)
  rogue_scores <- check_cluster_quality(submat, unique(cluster_predicted), platform)
  list(clustering_result = cluster_predicted, rogue_scores = rogue_scores)
}


# use scDesign2 to fit model and simulate data with cluster label -------------------------
fit_marginals_new <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                          pval_cutoff = 0.05, epsilon = 1e-5,
                          jitter = TRUE, DT = TRUE){
  p <- nrow(x)
  n <- ncol(x)

  marginal <- match.arg(marginal)
  if(marginal == 'auto_choose'){
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v){
        mle_Poisson <- glm(gene ~ 1, family = poisson)
        tryCatch({
          mle_ZIP <- zeroinfl(gene ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }else{
        mle_NB <- glm.nb(gene ~ 1)
        if(min(gene) > 0)
          c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
        else
          tryCatch({
            mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
            chisq_val <- 2 * (logLik(mle_ZINB) - logLik(mle_NB))
            pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
            if(pvalue < pval_cutoff)
              c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
            else
              c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          },
          error = function(cond){
            c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          })
      }
    }))
  }else if(marginal == 'zinb'){
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v)
      {
        mle_Poisson <- glm(gene ~ 1, family = poisson)
        tryCatch({
          mle_ZIP <- zeroinfl(gene ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }
      else
      {
        if(min(gene) > 0)
        {
          mle_NB <- glm.nb(gene ~ 1)
          c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
        }
        else
          tryCatch({
            mle_ZINB <- zeroinfl(gene ~ 1|1, dist = 'negbin')
            c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
          },
          error = function(cond){
            mle_NB <- glm.nb(gene ~ 1)
            c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          })
      }
    }))
  }else if(marginal == 'nb'){
    params <- t(apply(x, 1, function(gene){
      m <- mean(gene)
      v <- var(gene)
      if(m >= v){
        c(0.0, Inf, m)
      }else{
        mle_NB <- glm.nb(gene ~ 1)
        c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
      }
    }))
  }else if(marginal == 'poisson'){
    params <- t(apply(x, 1, function(gene){
      c(0.0, Inf, mean(gene))
    }))
  }

  if(DT){
    u <- t(sapply(1:p, function(iter){
      param <- params[iter, ]
      gene <- unlist(x[iter, ])
      prob0 <- param[1]
      u1 <- prob0 + (1 - prob0) * pnbinom(gene, size = param[2], mu = param[3])
      u2 <- (prob0 + (1 - prob0) * pnbinom(gene - 1, size = param[2], mu = param[3])) *
        as.integer(gene > 0)
      if(jitter)
        v <- runif(n)
      else
        v <- rep(0.5, n)
      r <- u1 * v + u2 * (1 - v)
      idx_adjust <- which(1-r < epsilon)
      r[idx_adjust] <- r[idx_adjust] - epsilon
      idx_adjust <- which(r < epsilon)
      r[idx_adjust] <- r[idx_adjust] + epsilon

      r
    }))
  }else{
    u <- NULL
  }

  return(list(params = params, u = u))
}

fit_Gaussian_copula_new <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                                jitter = TRUE, zp_cutoff = 0.8,
                                min_nonzero_num = 2){
  marginal <- match.arg(marginal)
  n <- ncol(x)
  p <- nrow(x)

  gene_zero_prop <- apply(x, 1, function(y){
    sum(y < 1e-5) / n
  })

  gene_sel1 <- which(gene_zero_prop < zp_cutoff)
  gene_sel2 <- which(gene_zero_prop < 1.0 - min_nonzero_num/n &
                       gene_zero_prop >= zp_cutoff)
  gene_sel3 <- (1:p)[-c(gene_sel1, gene_sel2)]

  if(length(gene_sel1) > 0){
    marginal_result1 <- fit_marginals_new(x[gene_sel1, , drop = FALSE], marginal, jitter = jitter, DT = TRUE)
    quantile_normal <- qnorm(marginal_result1$u)
    cov_mat <- cor(t(quantile_normal))
  }else{
    cov_mat = NULL
    marginal_result1 = NULL
  }

  if(length(gene_sel2) > 0){
    marginal_result2 <- fit_marginals_new(x[gene_sel2, , drop = FALSE], marginal, DT = FALSE)
  }else{
    marginal_result2 = NULL
  }
  return(list(cov_mat = cov_mat, marginal_param1 = marginal_result1$params,
              marginal_param2 = marginal_result2$params,
              gene_sel1 = gene_sel1, gene_sel2 = gene_sel2, gene_sel3 = gene_sel3,
              zp_cutoff = zp_cutoff, min_nonzero_num = min_nonzero_num,
              sim_method = 'copula', n_cell = n, n_read = sum(x)))
}

fit_wo_copula_new <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                          jitter = TRUE, min_nonzero_num = 2){
  marginal <- match.arg(marginal)
  n <- ncol(x)
  p <- nrow(x)

  gene_zero_prop <- apply(x, 1, function(y){
    sum(y < 1e-5) / n
  })

  gene_sel1 <- which(gene_zero_prop < 1.0 - min_nonzero_num/n)
  gene_sel2 <- (1:p)[-gene_sel1]

  if(length(gene_sel1) > 0){
    marginal_result1 <- fit_marginals_new(x[gene_sel1, ], marginal, jitter = jitter, DT = FALSE)
  }else{
    marginal_result1 = NULL
  }

  return(list(marginal_param1 = marginal_result1$params,
              gene_sel1 = gene_sel1, gene_sel2 = gene_sel2,
              min_nonzero_num = min_nonzero_num, sim_method = 'ind',
              n_cell = n, n_read = sum(x)))
}


fit_model_scDesign2_new <- function(data_mat, cell_type_sel, sim_method = c('copula', 'ind'),
                                marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                                jitter = TRUE, zp_cutoff = 0.8,
                                min_nonzero_num = 2, ncores = 1){
  sim_method <- match.arg(sim_method)
  marginal <- match.arg(marginal)

  if(sum(abs(data_mat - round(data_mat))) > 1e-5){
    warning('The entries in the input matrix are not integers. Rounding is performed.')
    data_mat <- round(data_mat)
  }

  if(sim_method == 'copula'){
    param <- mclapply(1:length(cell_type_sel), function(iter){
      fit_Gaussian_copula_new(data_mat[, colnames(data_mat) == cell_type_sel[iter]], marginal,
                          jitter = jitter, zp_cutoff = zp_cutoff,
                          min_nonzero_num = min_nonzero_num)
    }, mc.cores = ncores)
  }else if(sim_method == 'ind'){
    param <- mclapply(1:length(cell_type_sel), function(iter){
      fit_wo_copula_new(data_mat[, colnames(data_mat) == cell_type_sel[iter]], marginal,
                    jitter = jitter,
                    min_nonzero_num = min_nonzero_num)
    }, mc.cores = ncores)
  }

  names(param) <- cell_type_sel
  return(param)
}

get_scDesign2_result <- function(count_mat){
  n_cell_new <- ncol(count_mat)
  cell_type_sel <- unique(colnames(count_mat))
  cell_type_prop <- table(colnames(count_mat))[cell_type_sel]

  copula_result <- fit_model_scDesign2_new(count_mat, cell_type_sel, sim_method = 'copula',
                                       ncores = length(cell_type_sel))
  sim_count_copula <- scDesign2::simulate_count_scDesign2(copula_result, n_cell_new, sim_method = 'copula',
                                               cell_type_prop = cell_type_prop)
  list(copula_result = copula_result, sim_count_copula = sim_count_copula)
}


###############################################################################
###############################################################################
# samplename <- "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.UMI.countmatrix"
# cluster_prestep <- TRUE # FALSE for comple
samplename <- "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.UMI.COMPLE.countmatrix"
cluster_prestep <- FALSE # FALSE for comple

# samplename <- "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.countmatrix.scDesign2functions"
sample.format <- "txt"
directory <- "/home/gayan/Projects/scATAC_Simulator/results/20221011_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT"
out.directory <- "/home/gayan/Projects/scATAC_Simulator/results/20221013_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_VaryNewCell"
dir.create(out.directory)

## Model Fitting
cluster_prestep <- as.integer(cluster_prestep)
# Generate synthetic count matrix
## Read in count matrix
cat(sprintf("Reading count matrix %s.%s...\n", samplename, sample.format))
count_matrix <- read.table(sprintf("%s/%s.txt", directory, samplename), sep="\t",header = FALSE)
matrix_num <- data.matrix(count_matrix[,2:ncol(count_matrix)])
# count_pergene_vec <- rowSums(matrix_num)
# results.out.directory <- sprintf("%s/simulation_output", out.directory)
# dir.create(results.out.directory)
# write.table(count_pergene_vec, sprintf("%s/%s.real.nReadRegionmargional.txt",out.directory, samplename), row.names = FALSE,col.names = FALSE)

## Clustering
if (cluster_prestep==TRUE){
  cat("Louvain Clustering Before Simulation...\n")
  set.seed(2022)
  clustering_result <- get_cluster_seurat(matrix_num, 'UMI')
  # clustering_result <- unlist(read.table('/home/gayan/Projects/scATAC_Simulator/results/20220116_e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1_NONINPUT_withCluster/e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.countmatrix.LouvainClusterResults.txt'))
  print(clustering_result$rogue_scores)
  print(mean(clustering_result$rogue_scores))
  colnames(matrix_num) <- clustering_result$clustering_result
  write.table(clustering_result$clustering_result, file=sprintf("%s/%s.LouvainClusterResults.txt", out.directory, samplename), sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
} else {
  colnames(matrix_num) <- rep("Cluster1", ncol(matrix_num))
}

## Model Fitting
cat("Fitting Model of real data...\n")
set.seed(2022)
n_cell_new <- ncol(matrix_num)
cell_type_sel <- unique(colnames(matrix_num))
cell_type_prop <- table(colnames(matrix_num))[cell_type_sel]
copula_result <- fit_model_scDesign2_new(matrix_num, cell_type_sel, sim_method = 'copula',
                                       ncores = length(cell_type_sel))
save(copula_result, file=sprintf("%s/%s.scDesign2.modelfitting.Rdata", out.directory, samplename))

## Generate synthetic count matrix
load(sprintf("%s/%s.scDesign2.modelfitting.Rdata", out.directory, samplename))
## Read in count matrix
cat(sprintf("Reading count matrix %s.%s...\n", samplename, sample.format))
# data_mat <- read.table(sprintf("%s/%s.%s", directory, samplename, sample.format), sep="\t",header = TRUE)
# data_mat <- data_mat[,-1]
data_mat <- matrix_num

clustering_result <- unlist(read.table(sprintf("%s/%s.LouvainClusterResults.txt", out.directory, samplename), sep="\t",header = FALSE),use.names=FALSE)
cell_type_prop <- table(clustering_result)
total_count_old <- sum(data_mat)
n_cell_old <- ncol(data_mat)

## Vary read depth
adj_factor <- c(1/8, 1/4, 1/2, 1, 2, 4, 8)
set.seed(1)
print('simulating count (vary seq depth) ......\n')
sim_count <- mclapply(1:length(adj_factor), function(iter){
  scDesign2::simulate_count_scDesign2(copula_result, 
                           total_count_old = total_count_old,
                           n_cell_old = n_cell_old,
                           total_count_new = round(adj_factor[iter] * total_count_old),
                           n_cell_new = n_cell_old,
                           cell_type_prop = cell_type_prop,
                           reseq_method = 'mean_scale', cell_sample = TRUE)
}, mc.cores = length(adj_factor))
save(sim_count, file=sprintf('%s/sim_count_vary_seq_depth_%s.Rdata', out.directory, samplename))
sim_count_summary <- data.frame(seq_depth = sapply(sim_count, sum),
                                n_cell = sapply(sim_count, ncol))
sim_count_label <- lapply(sim_count, colnames)
save(sim_count_summary, sim_count_label, file=sprintf('%s/sim_count_info_vary_seq_depth_%s.Rdata', out.directory, samplename))


## Vary cell number 
set.seed(1)
print('simulating count (vary cell number) ......\n')
sim_count <- mclapply(1:length(adj_factor), function(iter){
  scDesign2::simulate_count_scDesign2(copula_result, 
                           total_count_old = total_count_old,
                           n_cell_old = n_cell_old,
                           total_count_new = total_count_old,
                           n_cell_new = round(adj_factor[iter] * n_cell_old),
                           cell_type_prop = cell_type_prop,
                           reseq_method = 'mean_scale', cell_sample = TRUE)
}, mc.cores = length(adj_factor))

save(sim_count, file=sprintf('%s/sim_count_vary_cell_number_%s.Rdata', out.directory, samplename))

sim_count_summary <- data.frame(seq_depth = sapply(sim_count, sum),
                                n_cell = sapply(sim_count, ncol))
sim_count_label <- lapply(sim_count, colnames)
save(sim_count_summary, sim_count_label, file=sprintf('%s/sim_count_info_vary_cell_number_%s.Rdata', out.directory, samplename))


## Write out tables
load(sprintf('%s/sim_count_vary_seq_depth_%s.Rdata', out.directory, samplename))
lapply(1:length(sim_count), function(iter){
    simu_matrix <- sim_count[[iter]]
    write.table(simu_matrix, sprintf("%s/%s.VarySeqDepth%s.scDesign2Simulated.%s", out.directory, samplename, adj_factor[iter], sample.format), sep="\t", row.names = TRUE,col.names = TRUE)
    # simu_count_pergene_vec <- rowSums(simu_matrix)
    # write.table(simu_count_pergene_vec, sprintf("%s/%s.VarySeqDepth%s.scDesign2Simulated.nReadRegionmargional.%s", out.directory, samplename, adj_factor[iter], sample.format), row.names = FALSE,col.names = FALSE)
  })

load(sprintf('%s/sim_count_vary_cell_number_%s.Rdata', out.directory, samplename))
lapply(1:length(sim_count), function(iter){
    simu_matrix <- sim_count[[iter]]
    write.table(simu_matrix, sprintf("%s/%s.VaryCellNumber%s.scDesign2Simulated.%s", out.directory, samplename, adj_factor[iter], sample.format), sep="\t", row.names = TRUE,col.names = TRUE)
    # simu_count_pergene_vec <- rowSums(simu_matrix)
    # write.table(simu_count_pergene_vec, sprintf("%s/%s.VaryCellNumber%s.scDesign2Simulated.nReadRegionmargional.%s", out.directory, samplename, adj_factor[iter], sample.format), row.names = FALSE,col.names = FALSE)
  })


### Test
out.directory <- "/home/gayan/Projects/scATAC_Simulator/results/20221013_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_VaryNewCell"
samplename <- "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.UMI.countmatrix"
load(sprintf('%s/sim_count_vary_cell_number_%s.Rdata', out.directory, samplename))
lapply(1:length(sim_count), function(iter){
    simu_matrix <- sim_count[[iter]]
    dim(simu_matrix)
     })
