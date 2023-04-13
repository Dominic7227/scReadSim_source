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

### Functions


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

# Calculate marginal count for each region
# Arguments
## samplename: Count matrix file name (without format)
## sample.format: Count matrix file format, like csv
## directory: Count matrix directory
## out.directory: output directory
## celltype_input: selected cell types for analysis
args <- commandArgs(trailingOnly = TRUE) # Extract all arguments from Linux input

samplename <- args[1]
sample.format <- args[2]
directory <- args[3]
out.directory <- args[4]
celltype_input <- args[5]

dir.create(out.directory)

## Read in count matrix
cat(sprintf("Reading count matrix %s.%s...\n", samplename, sample.format))
count_matrix <- read.table(sprintf("%s/%s.%s", directory, samplename, sample.format), sep="\t",header = FALSE)
matrix_num <- data.matrix(count_matrix[,2:ncol(count_matrix)])
count_pergene_vec <- rowSums(ceiling(matrix_num/2))
# results.out.directory <- sprintf("%s/simulation_output", out.directory)
# dir.create(results.out.directory)
write.table(count_pergene_vec, sprintf("%s/%s.real.nPairsRegionmargional.txt",out.directory, samplename), row.names = FALSE,col.names = FALSE)

## Read in cell types
clustering_result <- read.table(celltype_input, sep="\t")[-1,11]
print(which(table(clustering_result) == 1))
clustering_result[which(clustering_result == "Podocytes")] <- "Unknown"
clustering_result[which(clustering_result == "Endothelial I (glomerular)")] <- "Unknown"
celltype_sel_show <- c("Erythroblasts", "Hematopoietic progenitors", "Monocytes", "Immature B cells")

## Use scDesign2 for training countmatrix
colnames(matrix_num) <- clustering_result
# simu_matrix_selected <- simu_matrix[,colnames(simu_matrix)%in%celltype_sel_show]
matrix_num_selected <- matrix_num[,colnames(matrix_num)%in%celltype_sel_show]
colnames(matrix_num_selected) <- unlist(lapply(strsplit(gsub('[[:digit:]]+', '', colnames(matrix_num_selected)),".", fixed=TRUE), function(x){return(paste(x, collapse = ' '))}))

# simu_matrix_scDesign2 <- get_scDesign2_result(matrix_num)
set.seed(2022)
n_cell_new <- ncol(matrix_num_selected)
cell_type_sel <- unique(colnames(matrix_num_selected))
cell_type_prop <- table(colnames(matrix_num_selected))[cell_type_sel]
copula_result <- fit_model_scDesign2_new(matrix_num_selected, cell_type_sel, sim_method = 'copula',
                                       ncores = length(cell_type_sel))
simu_matrix <- scDesign2::simulate_count_scDesign2(copula_result, n_cell_new, sim_method = 'copula',
                                               cell_type_prop = cell_type_prop)
# simu_matrix <- simu_matrix_scDesign2$sim_count_copula
rownames(simu_matrix) <- count_matrix[,1]
write.table(simu_matrix, sprintf("%s/%s.scDesign2Simulated.%s", out.directory, samplename, sample.format), sep="\t", row.names = TRUE,col.names = TRUE)
simu_count_pergene_vec <- rowSums(ceiling(simu_matrix/2))
write.table(simu_count_pergene_vec, sprintf("%s/%s.scDesign2Simulated.nPairsRegionmargional.txt", out.directory, samplename), row.names = FALSE,col.names = FALSE)
save(count_matrix, simu_matrix, file=sprintf("%s/%s.scDesign2Simulated.Rdata", out.directory, samplename))
















