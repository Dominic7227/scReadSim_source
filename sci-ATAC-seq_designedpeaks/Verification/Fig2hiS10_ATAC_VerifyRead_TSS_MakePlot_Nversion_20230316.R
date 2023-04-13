library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(plyr)
library(GenomicRanges)

## Set directory
out.directory <- "/home/gayan/Projects/scATAC_Simulator/results/20221130_BoneMarrow_62016_chr1_INPUT_withCluster_mm9TSS"
# out.directory <- "/Users/gayan/Dropbox/Research/ATACseq_peak_calling/Fragment Simulator/Results/20220407_e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1_NONINPUT"
fig.directory <- paste0(out.directory, "/", "VerifyReadsPlots")
filename <- "BoneMarrow_62016_chr1"
dir.create(fig.directory)

############# Plot 1 kmer frequencies ############# 
for (kmer_length in c(11,21,31,41)){
  kmer_count_real_read1 <- read.table(sprintf("%s/VerifyCoverage/kmer_real/mer%s_counts_read1_hist.txt", out.directory, kmer_length), header = FALSE)
  kmer_count_real_read2 <- read.table(sprintf("%s/VerifyCoverage/kmer_real/mer%s_counts_read2_hist.txt", out.directory, kmer_length), header = FALSE)
  kmer_count_scReadSim_read1 <- read.table(sprintf("%s/VerifyCoverage/kmer_scReadSim/mer%s_counts_read1_hist.txt", out.directory, kmer_length), header = FALSE)
  kmer_count_scReadSim_read2 <- read.table(sprintf("%s/VerifyCoverage/kmer_scReadSim/mer%s_counts_read2_hist.txt", out.directory, kmer_length), header = FALSE)
  colnames(kmer_count_real_read1) <- c("frequency", "multiplicity")
  colnames(kmer_count_real_read2) <- c("frequency", "multiplicity")
  colnames(kmer_count_scReadSim_read1) <- c("frequency", "multiplicity")
  colnames(kmer_count_scReadSim_read2) <- c("frequency", "multiplicity")
  
  kmer_count_real_df <- rbind(kmer_count_real_read1, kmer_count_real_read2)
  kmer_count_real_df$group <- c(rep("read 1", nrow(kmer_count_real_read1)), rep("read 2", nrow(kmer_count_real_read2)))
  kmer_count_scReadSim_df <- rbind(kmer_count_scReadSim_read1, kmer_count_scReadSim_read2)
  kmer_count_scReadSim_df$group <- c(rep("read 1", nrow(kmer_count_scReadSim_read1)), rep("read 2", nrow(kmer_count_scReadSim_read2)))
  
  kmer_count_df <- rbind(kmer_count_real_df, kmer_count_scReadSim_df)
  kmer_count_df$method <- c(rep("real", nrow(kmer_count_real_df)), c(rep("scReadSim", nrow(kmer_count_scReadSim_df))))
  assign(sprintf("mer%s_count_df", kmer_length), kmer_count_df)
  save(kmer_count_df, file=sprintf("%s/Plot_%smer.Rdata", fig.directory, kmer_length))
  
  p_kmer <- ggplot() +
    geom_line(data = kmer_count_df, aes(x = frequency, y = multiplicity, color = method, group = method), size=1.2, alpha= 0.5) +
    xlab(sprintf("%s-mer frequency x", kmer_length)) +
    ylab(sprintf("# unique %s-mer with frequency x", kmer_length)) +
    facet_wrap(~group) +
    scale_y_log10() +
    theme_bw() + 
    theme(legend.position="right",
          legend.title = element_blank(), 
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          legend.text=element_text(size=10),
          strip.text = element_text(size = 15),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # Output 31mer for main plot in manuscript
  if (kmer_length == 31) {
    p_kmermain <- p_kmer
  } else {
    p_kmermain <- NULL
  }
  ggsave(plot=p_kmer, sprintf("%s/%smer_libeplot.pdf", fig.directory, kmer_length), width=10, height = 5)
  rm(list = setdiff(ls(), c("out.directory", "filename", "p_kmermain","fig.directory", "mer11_count_df", "mer21_count_df","mer31_count_df","mer41_count_df")))
}

combined_mer_df <- rbind(mer11_count_df, mer21_count_df, mer31_count_df, mer41_count_df)
combined_mer_df$case <- c(rep("11-mer", nrow(mer11_count_df)), rep("21-mer", nrow(mer21_count_df)), rep("31-mer", nrow(mer31_count_df)), rep("41-mer", nrow(mer41_count_df)))
save(combined_mer_df, file=sprintf("%s/Plot_kmer_combined.Rdata", fig.directory))

load(file=sprintf("%s/Plot_kmer_combined.Rdata", fig.directory))
p_kmercombined <- ggplot() +
  ggrastr::rasterize(geom_point(data = combined_mer_df, aes(x = frequency, y = multiplicity, color = method, group = method), size=0.5, alpha= 0.1),dpi = 300) +
  # xlab("k-mer frequency x") +
  # ylab("# unique k-mers with frequency x") +
  labs(x=expression(paste("k-mer frequency ",italic("x"))),
       y=expression(paste("# k-mers with frequency ",italic("x")))) +
  facet_grid(group ~ case, scales="free_x") +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw() + 
  theme(legend.position="none",
        legend.title = element_blank(), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        strip.text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(plot=p_kmercombined, sprintf("%s/kmer_combined_pointplot.pdf", fig.directory), width=10, height=5)


# p_kmercombined <- ggplot() +
#   geom_bar(data = combined_mer_df, aes(x = frequency, y = multiplicity, color = method, group = method), stat="identity", alpha= 0.5) +
#   xlab("k-mer frequency x") +
#   ylab("# unique k-mer with frequency x") +
#   facet_grid(group~case) +
#   scale_y_log10() +
#   theme_bw() + 
#   theme(legend.position="none",
#         legend.title = element_blank(), 
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 15),
#         legend.text=element_text(size=10),
#         strip.text = element_text(size = 15),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# ggsave(plot=p_kmercombined, sprintf("%s/kmer_combined_barplot.pdf", fig.directory))


############# Plot 2 fragment length ############# 
fragment_size_real <- read.table(sprintf("%s/VerifyCoverage/deepTools_scReadSim/fragmentSize_Real.txt", out.directory), header = TRUE)
fragment_size_scReadSim <- read.table(sprintf("%s/VerifyCoverage/deepTools_scReadSim/fragmentSize_scReadSim.txt", out.directory), header = TRUE)

# fragment_size_real <- fragment_size %>% filter(Sample == "Real") 
fragment_size_real_vec <- rep(fragment_size_real$Size, fragment_size_real$Occurrences)
fragment_size_scReadSim_vec <- rep(fragment_size_scReadSim$Size, fragment_size_scReadSim$Occurrences)
fragment_size_vector_df <- data.frame(Size=c(fragment_size_real_vec, fragment_size_scReadSim_vec), method=c(rep("real", length(fragment_size_real_vec)), rep("scReadSim", length(fragment_size_scReadSim_vec))))
save(fragment_size_vector_df, file=sprintf("%s/Plot_FragSize.Rdata", fig.directory))

# p <- ggplot(fragment_size_vector_df, aes(x=Size, fill=method))+
#   geom_histogram(data=subset(fragment_size_vector_df,fragment_size_vector_df$method == 'Real'),aes(y=..density..,fill="Real"),color="#e9ecef", alpha=0.5,binwidth=10)+
#   geom_histogram(data=subset(fragment_size_vector_df,fragment_size_vector_df$method == 'scReadSim'),aes(y=..density.. ,fill="scReadSim"),color="#e9ecef", alpha=0.5,binwidth=10)+
#   scale_fill_manual(values = c("#F8766D", "#00BFC4"), name="") +
#   labs(x="Fragment Size", y = "Frequency")+   
#   theme_light() +
#   scale_x_continuous(limits = c(0, 1000)) +
#   theme(legend.text=element_text(size=14),
#         axis.text=element_text(size=14),
#         axis.title=element_text(size=14),
#         legend.position="bottom",
#         legend.title = element_text(""))

load(file=sprintf("%s/Plot_FragSize.Rdata", fig.directory))
p_FragLen <- ggplot(fragment_size_vector_df, aes(x=Size))+
  geom_histogram(data=subset(fragment_size_vector_df,fragment_size_vector_df$method == 'real'),aes(y=..density..,fill="real"),color="#e9ecef", alpha=1,binwidth=15)+
  geom_histogram(data=subset(fragment_size_vector_df,fragment_size_vector_df$method == 'scReadSim'),aes(y=..density.. ,fill="scReadSim"),color="#e9ecef", alpha=1,binwidth=15)+
  geom_density(aes(y=..density..), size=0.5) +
  scale_fill_manual(values = c("#F8766D", "#00BFC4"), name="") +
  facet_grid(~method) +
  labs(x="fragment size", y = "frequency")+   
  theme_light() +
  scale_x_continuous(limits = c(0, 1000)) +
  theme_bw() +
  theme(legend.position="bottom",
        legend.title = element_text(""),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        strip.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(plot=p_FragLen, sprintf("%s/FragmentSize.pdf", fig.directory))


#############  Plot 3 Error rate ############# 
error_rate_real <- read.table(sprintf("%s/VerifyQuality/fgbio/Real.error_rate_by_read_position.txt", out.directory), header = TRUE)
error_rate_scReadSim <- read.table(sprintf("%s/VerifyQuality/fgbio/scReadSim_manualError.error_rate_by_read_position.txt", out.directory), header = TRUE)

error_rate_real_df <- error_rate_real %>% select(c(read_number, position, error_rate))
error_rate_real_df$group <- factor(error_rate_real_df$read_number, levels=c(1,2), labels=c("read 1", "read 2"))
error_rate_scReadSim_df <- error_rate_scReadSim %>% select(c(read_number, position, error_rate))
error_rate_scReadSim_df$group <- factor(error_rate_scReadSim_df$read_number, levels=c(1,2), labels=c("read 1", "read 2"))

error_rate_df <- rbind(error_rate_real_df, error_rate_scReadSim_df)
error_rate_df$method <- c(rep("real", nrow(error_rate_real_df)), c(rep("scReadSim", nrow(error_rate_scReadSim_df))))

save(error_rate_df, file=sprintf("%s/Errorrate_plot.Rdata", fig.directory))

load(file=sprintf("%s/Errorrate_plot.Rdata", fig.directory))
p_errorrate <- ggplot() +
  geom_line(data = error_rate_df %>% filter(method=="real"), aes(x = position, y = error_rate, color = method, group = method), size=1.5, alpha=0.5) +
  geom_line(data = error_rate_df %>% filter(method=="scReadSim"), aes(x = position, y = error_rate, color = method, group = method), size=0.7, alpha=0.5) +  
  scale_x_continuous(limits = c(0, 49)) +
  xlab(sprintf("position in read")) +
  ylab("error rate") +
  facet_wrap(~group, scales = "free_x") +
  scale_y_log10() +
  theme_bw() + 
  theme(legend.position="none",
        legend.title = element_blank(), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        strip.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(plot=p_errorrate, sprintf("%s/Errorate_plot.pdf", fig.directory), width=10, height = 5)

#############  Plot 4 RPKM  ############# 
setwd(out.directory)
# Read in real feature set
ref_feature_set_dir <- sprintf("%s/%s.MACS3.bed_peaks.bed", out.directory, filename)
bg_ref_feature_set_dir <- sprintf("%s/%s.MACS3_thr10_peaks.COMPLE.bed", out.directory, filename)
ref_feature_set <- read.table(ref_feature_set_dir)
bg_ref_feature_set <- read.table(bg_ref_feature_set_dir)

# Read in input feature set
input_feature_set_dir <- sprintf("%s/%s.INPUT.peaks.bed", out.directory, filename)
bg_input_feature_set_dir <- sprintf("%s/%s.INPUT.COMPLE.peaks.bed", out.directory, filename)
input_feature_set <- read.table(input_feature_set_dir)
bg_input_feature_set<- read.table(bg_input_feature_set_dir)

# Read in synthetic MACS3 feature set
MACS3_feature_set <- read.table("BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.MACS3_peaks.narrowPeak.bed")
bg_MACS3_feature_set <- read.table("BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.MACS3_peaks.narrowPeak.COMPLE.bed")

# Read in synthetic SEACR feature set
SEACR_feature_set <- read.table("BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.SEACR_peaks.narrowPeak.bed")
bg_SEACR_feature_set <- read.table("BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.SEACR_peaks.narrowPeak.COMPLE.bed")

# Read in synthetic HOMER feature set
HOMER_feature_set <- read.table("BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HOMER.merged.bed")
bg_HOMER_feature_set <- read.table("BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HOMER.merged.COMPLE.bed")

# Read in synthetic HMMRATAC feature set
HMMRATAC_feature_set <- read.table("BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HMMRATAC.bed")
bg_HMMRATAC_feature_set <- read.table("BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HMMRATAC.COMPLE.bed")

# Read in real read coverage
count_pergene_vec <- read.table(sprintf("%s/BoneMarrow_62016_chr1.countmatrix.MarginalFeatureCount.txt", out.directory))
count_pergene_vec <- as.numeric(count_pergene_vec[,-1])
bg_count_pergene_vec <- read.table(sprintf("%s/BoneMarrow_62016_chr1.COMPLE.countmatrix.MarginalFeatureCount.txt", out.directory))
bg_count_pergene_vec <- as.numeric(bg_count_pergene_vec[,-1])

# Read in simulated read coverage
samplename <- "BoneMarrow_62016_chr1.countmatrix"
samplename_COMPLE <- "BoneMarrow_62016_chr1.COMPLE.countmatrix"
simu_count_pergene_vec <- read.table(sprintf("%s/%s.scDesign2Simulated.GroundTruth.MarginalFeatureCount.txt", out.directory, samplename))
simu_count_pergene_vec <- as.numeric(simu_count_pergene_vec[,-1])
bg_simu_count_pergene_vec  <- read.table(sprintf("%s/%s.scDesign2Simulated.GroundTruth.MarginalFeatureCount.txt",out.directory, samplename_COMPLE))
bg_simu_count_pergene_vec <- as.numeric(bg_simu_count_pergene_vec[,-1])

# Read in MACS3 synthetic read coverage
simu_count_MACS3peak <- read.table(sprintf("%s/%s.scDesign2Simulated.MACS3Peaks.MarginalFeatureCount.txt", out.directory, samplename))
simu_count_MACS3peak <- as.numeric(simu_count_MACS3peak[,-1])
bg_simu_count_MACS3peak <- read.table(sprintf("%s/%s.scDesign2Simulated.MACS3Peaks.MarginalFeatureCount.txt", out.directory, samplename_COMPLE))
bg_simu_count_MACS3peak <- as.numeric(bg_simu_count_MACS3peak[,-1])

# Read in SEACR synthetic read coverage
simu_count_SEACRpeak <- read.table(sprintf("%s/%s.scDesign2Simulated.SEACRPeaks.MarginalFeatureCount.txt", out.directory, samplename))
simu_count_SEACRpeak <- as.numeric(simu_count_SEACRpeak[,-1])
bg_simu_count_SEACRpeak <- read.table(sprintf("%s/%s.scDesign2Simulated.SEACRPeaks.MarginalFeatureCount.txt", out.directory, samplename_COMPLE))
bg_simu_count_SEACRpeak <- as.numeric(bg_simu_count_SEACRpeak[,-1])

# Read in HOMER synthetic read coverage
simu_count_HOMERpeak <- read.table(sprintf("%s/%s.scDesign2Simulated.HOMERPeaks.MarginalFeatureCount.txt", out.directory, samplename))
simu_count_HOMERpeak <- as.numeric(simu_count_HOMERpeak[,-1])
bg_simu_count_HOMERpeak <- read.table(sprintf("%s/%s.scDesign2Simulated.HOMERPeaks.MarginalFeatureCount.txt", out.directory, samplename_COMPLE))
bg_simu_count_HOMERpeak <- as.numeric(bg_simu_count_HOMERpeak[,-1])

# Read in HMMRATAC synthetic read coverage
simu_count_HMMRATACpeak <- read.table(sprintf("%s/%s.scDesign2Simulated.HMMRATACPeaks.MarginalFeatureCount.txt", out.directory, samplename))
simu_count_HMMRATACpeak <- as.numeric(simu_count_HMMRATACpeak[,-1])
bg_simu_count_HMMRATACpeak <- read.table(sprintf("%s/%s.scDesign2Simulated.HMMRATACPeaks.MarginalFeatureCount.txt", out.directory, samplename_COMPLE))
bg_simu_count_HMMRATACpeak <- as.numeric(bg_simu_count_HMMRATACpeak[,-1])

trun_95 <- function(vec){
  return(vec[vec < quantile(vec, 0.9)])
}

## Calculate read coverage density of reference feature for real data
read_density_real <- 1000 * (1000000 / sum(count_pergene_vec)) * count_pergene_vec / (ref_feature_set[,3] - ref_feature_set[,2])
bg_read_density_real <- 1000 * (1000000 / sum(bg_count_pergene_vec)) * bg_count_pergene_vec / (bg_ref_feature_set[,3] - bg_ref_feature_set[,2])
read_density_real_trun <- trun_95(read_density_real)
bg_read_density_real_trun <- trun_95(bg_read_density_real)

## Calculate read coverage density of input feature for simulated data
read_density_input <- 1000 * (1000000 / sum(simu_count_pergene_vec)) *  simu_count_pergene_vec / (input_feature_set[,3] - input_feature_set[,2])
bg_read_density_input <- 1000 * (1000000 / sum(bg_simu_count_pergene_vec)) * bg_simu_count_pergene_vec / (bg_input_feature_set[,3] - bg_input_feature_set[,2])
read_density_input_trun <- trun_95(read_density_input)
bg_read_density_input[is.na(bg_read_density_input)] <- 0
bg_read_density_input_trun <- trun_95(bg_read_density_input)

## Calculate read coverage density of MACS3 synthetic feature for synthetic data
read_density_MACS3 <- 1000 * (1000000 / sum(simu_count_MACS3peak)) *  simu_count_MACS3peak / (MACS3_feature_set[,3] - MACS3_feature_set[,2])
bg_read_density_MACS3 <- 1000 * (1000000 / sum(bg_simu_count_MACS3peak)) *  bg_simu_count_MACS3peak / (bg_MACS3_feature_set[,3] - bg_MACS3_feature_set[,2])
read_density_MACS3_trun <- trun_95(read_density_MACS3)
bg_read_density_MACS3_trun <- trun_95(bg_read_density_MACS3)

## Calculate read coverage density of SEACR synthetic feature for synthetic data
read_density_SEACR <- 1000 * (1000000 / sum(simu_count_SEACRpeak)) * simu_count_SEACRpeak / (SEACR_feature_set[,3] - SEACR_feature_set[,2])
bg_read_density_SEACR <- 1000 * (1000000 / sum(bg_simu_count_SEACRpeak)) * bg_simu_count_SEACRpeak / (bg_SEACR_feature_set[,3] - bg_SEACR_feature_set[,2])
read_density_SEACR_trun <- trun_95(read_density_SEACR)
bg_read_density_SEACR_trun <- trun_95(bg_read_density_SEACR)

## Calculate read coverage density of HOMER synthetic feature for synthetic data
read_density_HOMER <- 1000 * (1000000 / sum(simu_count_HOMERpeak)) * simu_count_HOMERpeak / (HOMER_feature_set[,3] - HOMER_feature_set[,2])
bg_read_density_HOMER <- 1000 * (1000000 / sum(bg_simu_count_HOMERpeak)) * bg_simu_count_HOMERpeak / (bg_HOMER_feature_set[,3] - bg_HOMER_feature_set[,2])
read_density_HOMER_trun <- trun_95(read_density_HOMER)
bg_read_density_HOMER_trun <- trun_95(bg_read_density_HOMER)

## Calculate read coverage density of HMMRATAC synthetic feature for synthetic data
read_density_HMMRATAC <- 1000 * (1000000 / sum(simu_count_HMMRATACpeak)) * simu_count_HMMRATACpeak / (HMMRATAC_feature_set[,3] - HMMRATAC_feature_set[,2])
bg_read_density_HMMRATAC <- 1000 * (1000000 / sum(bg_simu_count_HMMRATACpeak)) * bg_simu_count_HMMRATACpeak / (bg_HMMRATAC_feature_set[,3] - bg_HMMRATAC_feature_set[,2])
bg_read_density_HMMRATAC[is.na(bg_read_density_HMMRATAC)] <- 0
read_density_HMMRATAC_trun <- trun_95(read_density_HMMRATAC)
bg_read_density_HMMRATAC_trun <- trun_95(bg_read_density_HMMRATAC)


## Make boxplots
save(read_density_real, bg_read_density_real, file=sprintf("%s/BoneMarrow_62016_chr1.VeryRead.ReadDensity.Real.Rdata", fig.directory))
save(read_density_input, bg_read_density_input, file=sprintf("%s/BoneMarrow_62016_chr1.VeryRead.ReadDensity.Input.Rdata", fig.directory))
save(read_density_MACS3, bg_read_density_MACS3, file=sprintf("%s/BoneMarrow_62016_chr1.VeryRead.ReadDensity.Input.MACS3.Rdata", fig.directory))
save(read_density_SEACR, bg_read_density_SEACR, file=sprintf("%s/BoneMarrow_62016_chr1.VeryRead.ReadDensity.Input.SEACR.Rdata", fig.directory))
save(read_density_HOMER, bg_read_density_HOMER, file=sprintf("%s/BoneMarrow_62016_chr1.VeryRead.ReadDensity.Input.HOMER.Rdata", fig.directory))
save(read_density_HMMRATAC, bg_read_density_HMMRATAC, file=sprintf("%s/BoneMarrow_62016_chr1.VeryRead.ReadDensity.Input.HMMRATAC.Rdata", fig.directory))

create_dataftame <- function(inform_RPKM, bg_RPKM, group_name){
  df_tmp <- data.frame(read_density=c(inform_RPKM, bg_RPKM), 
                       group=factor(rep(group_name, length(inform_RPKM) + length(bg_RPKM))),
                       category=factor(c(rep("Foreground Features", length(inform_RPKM)), rep("Background Features", length(bg_RPKM))), levels = c("Foreground Features", "Background Features"), ordered=TRUE))
  return(df_tmp)
}
df_INPUT_trun <- create_dataftame(read_density_input_trun, bg_read_density_input_trun, "User Input")
df_INPUT_Real <- create_dataftame(read_density_real_trun, bg_read_density_real_trun, "Real")
df_INPUT_MACS3 <- create_dataftame(read_density_MACS3_trun, bg_read_density_MACS3_trun, "MACS3")
df_INPUT_SEACR <- create_dataftame(read_density_SEACR_trun, bg_read_density_SEACR_trun, "SEACR")
df_INPUT_HOMER <- create_dataftame(read_density_HOMER_trun, bg_read_density_HOMER_trun, "HOMER")
df_INPUT_HMMRATAC <- create_dataftame(read_density_HMMRATAC_trun, bg_read_density_HMMRATAC_trun, "HMMRATAC")
list_trun <- list(df_INPUT_trun, df_INPUT_Real, df_INPUT_MACS3, df_INPUT_SEACR, df_INPUT_HOMER, df_INPUT_HMMRATAC)
df_trun <- Reduce(rbind, list_trun)

save(df_trun, file=sprintf("%s/RKPM_input.Rdata", fig.directory))


load(sprintf("%s/RKPM_input.Rdata", fig.directory))
measures1 <-  c("Foreground Features", "Background Features")
measures2 <-  c("Peak", "Non-Peak")
df_trun$category <- mapvalues(df_trun$category, from = measures1, to = measures2)
read_des_p_MACS3 <- df_trun %>% filter(group %in% c("Real", "User Input") ) %>% 
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


read_des_p <- df_trun %>% filter(group != "Real") %>%
  ggplot( aes(x=category, y=read_density, fill=category, color=category)) + 
  stat_boxplot(geom = "errorbar",width=0.3,size=1)+
  geom_boxplot(size=1, outlier.size=0.3) +
  scale_fill_manual(values=c("#D9D9D9", "#D9D9D9")) +
  scale_color_manual(values=c("black", "black")) +
  theme_bw() +
  theme(legend.position="none",
        legend.title=element_blank(),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text=element_text(size=25),
        strip.text = element_text(size = 25),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 25, angle = 45,hjust = 1),
        panel.grid.minor = element_blank()) +
  # theme(axis.line = element_line(colour = "black"),
  #   panel.grid.major = element_blank(),
  #   panel.grid.minor = element_blank(),
  #   panel.background = element_blank(),
  #   panel.border = element_blank(),
  #   axis.text = element_text(size = 12),
  #   axis.text.x=element_blank(), axis.title.x=element_blank(), axis.ticks.x=element_blank(),
  #   legend.position="bottom",
  #   legend.title=element_blank(),
  #   legend.text=element_text(size = 12),
  #   strip.text = element_text(size = 12)) +
xlab("") + ylab("RPKM") +
  facet_wrap(~group, scales = "free_x", ncol=3)

ggsave(plot=read_des_p, file=sprintf("%s/NBT_BoneMarrow_62016_chr1.VeryRead.RPKM.tworow.20230316.pdf",out.directory), width=9, height = 8.4)


#############  Plot 5 ROC curve  ############# 
input_peakfile <- sprintf("%s/BoneMarrow_62016_chr1.INPUT.peaks.bed", out.directory)
synthetic_MACS3_peakfile <- sprintf("%s/BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.MACS3_peaks.narrowPeak", out.directory)
synthetic_SEACR_peakfile <- sprintf("%s/BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.SEACR_peaks.narrowPeak.bed", out.directory)
synthetic_HOMER_peakfile <- sprintf("%s/BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HOMER", out.directory)
synthetic_HMMRATAC_peakfile <- sprintf("%s/BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HMMRATAC.bed_peaks.gappedPeak", out.directory)

input_peak <- read.table(input_peakfile, sep="\t")
colnames(input_peak) <- c("chr", "start","end")
synthetic_peak_MACS3 <- read.table(synthetic_MACS3_peakfile, sep="\t")
colnames(synthetic_peak_MACS3)[1:3] <- c("chr", "start","end")
synthetic_peak_SEACR <- read.table(synthetic_SEACR_peakfile, sep="\t")
colnames(synthetic_peak_SEACR)[1:3] <- c("chr", "start","end")
synthetic_peak_HOMER <- read.table(synthetic_HOMER_peakfile, sep="\t")[,-1]
colnames(synthetic_peak_HOMER)[1:3] <- c("chr", "start","end")
synthetic_peak_HMMRATAC <- read.table(synthetic_HMMRATAC_peakfile, sep="\t")
colnames(synthetic_peak_HMMRATAC)[1:3] <- c("chr", "start","end")
# Sort
synthetic_peak_MACS3_sorted <- synthetic_peak_MACS3[order(synthetic_peak_MACS3[,8], decreasing = TRUE),] # sort by p-value
# overlap_size_MACS3 <- min(synthetic_peak_MACS3_sorted[,3] - synthetic_peak_MACS3_sorted[,2])/2
synthetic_peak_SEACR_sorted <- synthetic_peak_SEACR[order(synthetic_peak_SEACR[,4], decreasing = TRUE),] # sort by read sum
# overlap_size_SEACR <- min(synthetic_peak_SEACR_sorted[,3] - synthetic_peak_SEACR_sorted[,2])/2
# synthetic_peak_SEACR_sorted <- synthetic_peak_SEACR[order((synthetic_peak_SEACR[,4]/ (synthetic_peak_SEACR[,3] - synthetic_peak_SEACR[,2] + 1)), decreasing = TRUE),] # sort by average read density
synthetic_peak_HOMER_sorted <- synthetic_peak_HOMER[order(synthetic_peak_HOMER[,9], decreasing = FALSE),] # sort by p-value, increasing!!
# overlap_size_HOMER <- min(synthetic_peak_HOMER_sorted[,3] - synthetic_peak_HOMER_sorted[,2])/2
synthetic_peak_HMMRATAC_sorted <- synthetic_peak_HMMRATAC[order(synthetic_peak_HMMRATAC[,13], decreasing = TRUE),] # sort by HMMRATAC peak score
# overlap_size_HMMRATAC <- min(synthetic_peak_HMMRATAC_sorted[,3] - synthetic_peak_HMMRATAC_sorted[,2])/2
## Use half of minimum ground truth peak size
overlap_size <- min(input_peak[,3] - input_peak[,2])/2

### Function: peak set comparison
PeakCompare <- function(synthetic_peak_trun, real_peak, TN, overlap_size){
  real_peak_GR <- makeGRangesFromDataFrame(as.data.frame(real_peak))
  synthetic_peak_GR_trun <- makeGRangesFromDataFrame(as.data.frame(synthetic_peak_trun))
  
  # TP_vec <- rep(0, nrow(synthetic_peak_trun))
  # for (i in seq(nrow(synthetic_peak_trun))){
  #   TP_vec[i] <- sum(countOverlaps(makeGRangesFromDataFrame(as.data.frame(synthetic_peak_trun[i,])), real_peak_GR, minoverlap=(synthetic_peak_trun[i,3]-synthetic_peak_trun[i,2])/2) >0) > 0
  # }
  # TP <- sum(TP_vec)
  TP <- sum(countOverlaps(synthetic_peak_GR_trun, real_peak_GR, minoverlap = overlap_size) >0)
  TPR <- TP / nrow(real_peak)
  FP <- nrow(synthetic_peak_trun) - TP
  FPR <- FP/ (FP + TN)
  FDR <- FP/nrow(synthetic_peak_trun)
  return(list(TPR=TPR, FPR=FPR, FDR=FDR))
}

rocPeakCalculator <- function(synthetic_peak_sorted, real_peak, overlap_size){
  # synthetic_peak_sorted: synthetic peak list with power decreasing order
  # Count total regions number
  real_peak_GR <- makeGRangesFromDataFrame(as.data.frame(real_peak))
  synthetic_peak_GR <- makeGRangesFromDataFrame(as.data.frame(synthetic_peak_sorted))
  TP_max <- sum(countOverlaps(synthetic_peak_GR, real_peak_GR) >0) # Largest number of TP
  n_total <- nrow(real_peak) + nrow(synthetic_peak_sorted) - TP_max # Some of the positive regions will become true negative because peaks identified decreases
  
  threshold_quantile_list <- seq(0.001,1,by=0.01)
  FPR_vec <- rep(0, length(threshold_quantile_list))
  TPR_vec <- rep(0, length(threshold_quantile_list))
  for (threshold_quantile in threshold_quantile_list){
    TN <- round(threshold_quantile * nrow(synthetic_peak_sorted))
    synthetic_peak_trun <- synthetic_peak_sorted[1:(nrow(synthetic_peak_sorted) - TN),]
    roc_list <- PeakCompare(synthetic_peak_trun, real_peak, TN, overlap_size)
    FPR_vec[which(threshold_quantile_list==threshold_quantile)] <- roc_list[[2]]
    TPR_vec[which(threshold_quantile_list==threshold_quantile)] <- roc_list[[1]]
    # print(which(threshold_quantile_list==threshold_quantile)/length(seq(0.001,1,by=0.01)))
  }
  return(list(TPR_vec=TPR_vec, FPR_vec=FPR_vec))
}

fdrPeakCalculator <- function(synthetic_peak_sorted, real_peak, overlap_size){
  # synthetic_peak_sorted: synthetic peak list with power decreasing order
  # Count total regions number
  real_peak_GR <- makeGRangesFromDataFrame(as.data.frame(real_peak))
  synthetic_peak_GR <- makeGRangesFromDataFrame(as.data.frame(synthetic_peak_sorted))
  TP_max <- sum(countOverlaps(synthetic_peak_GR, real_peak_GR) >0) # Largest number of TP
  n_total <- nrow(real_peak) + nrow(synthetic_peak_sorted) - TP_max # Some of the positive regions will become true negative because peaks identified decreases
  
  threshold_quantile_list <- seq(0.001,1,by=0.01)
  FDR_vec <- rep(0, length(threshold_quantile_list))
  TPR_vec <- rep(0, length(threshold_quantile_list))
  for (threshold_quantile in threshold_quantile_list){
    TN <- round(threshold_quantile * nrow(synthetic_peak_sorted))
    synthetic_peak_trun <- synthetic_peak_sorted[1:(nrow(synthetic_peak_sorted) - TN),]
    roc_list <- PeakCompare(synthetic_peak_trun, real_peak, TN, overlap_size)
    FDR_vec[which(threshold_quantile_list==threshold_quantile)] <- roc_list[[3]]
    TPR_vec[which(threshold_quantile_list==threshold_quantile)] <- roc_list[[1]]
  }
  return(list(TPR_vec=TPR_vec, FDR_vec=FDR_vec))
}

# Default finding TPR
MACS3_default_TPR <- PeakCompare(synthetic_peak_MACS3_sorted, input_peak, TN=0, overlap_size)[[1]]
SEACR_default_TPR <- PeakCompare(synthetic_peak_SEACR_sorted, input_peak, TN=0, overlap_size)[[1]]
HOMER_default_TPR <- PeakCompare(synthetic_peak_HOMER_sorted, input_peak, TN=0, overlap_size)[[1]]
HMMRATAC_default_TPR <- PeakCompare(synthetic_peak_HMMRATAC_sorted, input_peak, TN=0, overlap_size)[[1]]
TPR_vec <- c(MACS3_default_TPR, SEACR_default_TPR, HOMER_default_TPR, HMMRATAC_default_TPR)
MACS3_default_FDR <- PeakCompare(synthetic_peak_MACS3_sorted, input_peak, TN=0, overlap_size)[[3]]
SEACR_default_FDR <- PeakCompare(synthetic_peak_SEACR_sorted, input_peak, TN=0, overlap_size)[[3]]
HOMER_default_FDR <- PeakCompare(synthetic_peak_HOMER_sorted, input_peak, TN=0, overlap_size)[[3]]
HMMRATAC_default_FDR <- PeakCompare(synthetic_peak_HMMRATAC_sorted, input_peak, TN=0, overlap_size)[[3]]
FDR_vec <- c(MACS3_default_FDR, SEACR_default_FDR, HOMER_default_FDR, HMMRATAC_default_FDR)
F1_vec <- 2 / (1/(1-FDR_vec) + 1/TPR_vec)

# ROC list
# MACS3_roc_list <- rocPeakCalculator(synthetic_peak_MACS3_sorted, input_peak, overlap_size_MACS3)
# SEACR_roc_list <- rocPeakCalculator(synthetic_peak_SEACR_sorted, input_peak, overlap_size_SEACR)
# HOMER_roc_list <- rocPeakCalculator(synthetic_peak_HOMER_sorted, input_peak, overlap_size_HOMER)
# HMMRATAC_roc_list <- rocPeakCalculator(synthetic_peak_HMMRATAC_sorted, input_peak, overlap_size_HMMRATAC)
MACS3_roc_list <- rocPeakCalculator(synthetic_peak_MACS3_sorted, input_peak, overlap_size)
SEACR_roc_list <- rocPeakCalculator(synthetic_peak_SEACR_sorted, input_peak, overlap_size)
HOMER_roc_list <- rocPeakCalculator(synthetic_peak_HOMER_sorted, input_peak, overlap_size)
HMMRATAC_roc_list <- rocPeakCalculator(synthetic_peak_HMMRATAC_sorted, input_peak, overlap_size)

# # FPR list
# MACS3_FDR_list <- fdrPeakCalculator(synthetic_peak_MACS3_sorted, input_peak, overlap_size_MACS3)
# SEACR_FDR_list <- fdrPeakCalculator(synthetic_peak_SEACR_sorted, input_peak, overlap_size_SEACR)
# HOMER_FDR_list <- fdrPeakCalculator(synthetic_peak_HOMER_sorted, input_peak, overlap_size_HOMER)
# HMMRATAC_FDR_list <- fdrPeakCalculator(synthetic_peak_HMMRATAC_sorted, input_peak, overlap_size_HMMRATAC)
MACS3_FDR_list <- fdrPeakCalculator(synthetic_peak_MACS3_sorted, input_peak, overlap_size)
SEACR_FDR_list <- fdrPeakCalculator(synthetic_peak_SEACR_sorted, input_peak, overlap_size)
HOMER_FDR_list <- fdrPeakCalculator(synthetic_peak_HOMER_sorted, input_peak, overlap_size)
HMMRATAC_FDR_list <- fdrPeakCalculator(synthetic_peak_HMMRATAC_sorted, input_peak, overlap_size)

# # Data frame for ggplot2
roc_df <- data.frame(FPR=c(MACS3_roc_list$FPR_vec, SEACR_roc_list$FPR_vec, HOMER_roc_list$FPR_vec, HMMRATAC_roc_list$FPR_vec),
                     TPR=c(MACS3_roc_list$TPR_vec, SEACR_roc_list$TPR_vec, HOMER_roc_list$TPR_vec, HMMRATAC_roc_list$TPR_vec),
                     group=c(rep("MACS3", length(MACS3_roc_list$FPR_vec)), rep("SEACR", length(SEACR_roc_list$FPR_vec)), rep("HOMER", length(HOMER_roc_list$FPR_vec)), rep("HMMRATAC", length(HMMRATAC_roc_list$FPR_vec))))
FDR_df <- data.frame(FDR=c(MACS3_FDR_list$FDR_vec, SEACR_FDR_list$FDR_vec, HOMER_FDR_list$FDR_vec, HMMRATAC_FDR_list$FDR_vec),
                     TPR=c(MACS3_FDR_list$TPR_vec, SEACR_FDR_list$TPR_vec, HOMER_FDR_list$TPR_vec, HMMRATAC_FDR_list$TPR_vec),
                     group=c(rep("MACS3", length(MACS3_FDR_list$FDR_vec)), rep("SEACR", length(SEACR_FDR_list$FDR_vec)), rep("HOMER", length(HOMER_FDR_list$FDR_vec)), rep("HMMRATAC", length(HMMRATAC_FDR_list$FDR_vec))))
roc_df$group <- factor(roc_df$group, levels=c("MACS3", "HOMER", "HMMRATAC", "SEACR"), ordered=TRUE)
FDR_df$group <- factor(FDR_df$group, levels=c("MACS3", "HOMER", "HMMRATAC", "SEACR"), ordered=TRUE)

save(roc_df, FDR_df, file=sprintf("%s/FDR_ROC_PeakCalling_plot_20221130.Rdata", fig.directory))
load(sprintf("%s/FDR_ROC_PeakCalling_plot_20221130.Rdata", fig.directory))

# thinned <- floor(seq(from=1,to=nrow(roc_df),length=20))

# Plot
color_vec <-c("#1b9e77",  "#d95f02", "#7570b3", "#e7298a")
roc_plot <- roc_df %>%
  ggplot( aes(x=FPR, y=TPR, group=group, color=group)) +
  geom_line(size=1.5, alpha=1) +
  #   geom_point(size=0.5) +
  # geom_line(size=0.5) +
  # scale_color_brewer(palette="Set3") +
  scale_color_manual(values=color_vec) +
  xlab("FPR") +
  ylab("TPR") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0, 1)) +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_blank(), 
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text=element_text(size=25),
        strip.text = element_text(size = 25))
# panel.grid.major = element_blank(),
# panel.grid.minor = element_blank())
# ggsave(roc_plot, file=sprintf("%s/ROC_benchmarkPeakCalling.pdf", fig.directory))
# Plot
PrecisionRecall_plot <- FDR_df %>%
  ggplot( aes(x=TPR, y=1-FDR, group=group, color=group)) +
  geom_line(size=1.5, alpha=1) +
  #   geom_point(size=0.5) +
  # geom_line(size=0.5) +
  # scale_color_brewer(palette="Set3") +
  scale_color_manual(values=color_vec) +
  ylab("Precision") +
  xlab("Recall") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0, 1)) +
  theme_bw() +
  theme(legend.position="bottom",
        legend.title = element_blank(), 
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text=element_text(size=25),
        strip.text = element_text(size = 25)) +
  guides(colour = guide_legend(override.aes = list(size=5)))
PrecisionRecall_plot_nolegend <- FDR_df %>%
  ggplot( aes(x=TPR, y=1-FDR, group=group, color=group)) +
  geom_line(size=1.5, alpha=1) +
  #   geom_point(size=0.5) +
  # geom_line(size=0.5) +
  # scale_color_brewer(palette="Set3") +
  scale_color_manual(values=color_vec) +
  ylab("Precision") +
  xlab("Recall") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0, 1)) +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_blank(), 
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        legend.text=element_text(size=25),
        strip.text = element_text(size = 25)) +
  guides(colour = guide_legend(override.aes = list(size=5)))

p2_legend <- get_legend(PrecisionRecall_plot)
# panel.grid.major = element_blank(),
# panel.grid.minor = element_blank())
# ggsave(PrecisionRecall_plot, file=sprintf("%s/FDRTPR_benchmarkPeakCalling.pdf", fig.directory))
library(ggpubr)
curve_plot <- ggarrange(roc_plot, PrecisionRecall_plot_nolegend, nrow=2,
                        common.legend = TRUE, legend="bottom", legend.grob=p2_legend)
ggsave(file=sprintf("%s/NBT_BoneMarrow_62016_chr1.VeryRead.ROC.20230406.pdf",fig.directory), curve_plot, height = 10, width = 8.4)



#############  Combine Figures  ############# 
gs <- list(arrangeGrob(p_kmercombined, left = textGrob("A", x = unit(1, "npc"), 
                                                       y = unit(.95, "npc"))),
           arrangeGrob(p_FragLen + theme(legend.position="none"), left = textGrob("B", x = unit(1, "npc"), 
                                                                                  y = unit(.95, "npc"))),
           arrangeGrob(p_errorrate, left = textGrob("C", x = unit(1, "npc"), 
                                                    y = unit(.95, "npc")))
           # arrangeGrob(read_des_p, left = textGrob("D", x = unit(1, "npc"), 
           #                                         y = unit(.95, "npc")))
)
lay <- rbind(c(1,1,1,1),
             c(1,1,1,1),
             c(1,1,1,1),
             c(2,2,3,3),
             c(2,2,3,3))

library(gtable)
legend <- gtable_filter(ggplot_gtable(ggplot_build(p_FragLen)), "guide-box")
p_final <- grid.arrange(arrangeGrob(grobs = gs, layout_matrix = lay,common),
                        legend,
                        heights=c(1.1, 0.1),
                        nrow = 2)
ggsave(file=sprintf("%s/%s.INPUT.VerifyRead.Combined_20230128.pdf",fig.directory, filename), p_final, height = 8, width = 9.6, dpi = 600)



# gs_peak <- list(arrangeGrob(read_des_p, left = textGrob("A", x = unit(1, "npc"), 
#                                                          y = unit(.95, "npc"))),
#             arrangeGrob(roc_plot, left = textGrob("B", x = unit(1, "npc"), 
#                                                          y = unit(.95, "npc"))),
#             arrangeGrob(FDR_plot, left = textGrob("C", x = unit(1, "npc"), 
#                                                          y = unit(.95, "npc"))))
gs_peak <- list(arrangeGrob(read_des_p, left = textGrob("A", x = unit(1, "npc"), 
                                                        y = unit(.95, "npc"))),
                arrangeGrob(curve_plot, left = textGrob("B", x = unit(1, "npc"), 
                                                        y = unit(.95, "npc"))))
lay_peak <- rbind(c(1,1,1,1),
                  c(2,2,3,3))
p_final_peak <-arrangeGrob(grobs = gs_peak, layout_matrix = lay_peak)
ggsave(file=sprintf("%s/%s.INPUT.VerifyRead.Peak.Combined_202201130.pdf",fig.directory, filename), p_final_peak, height = 4.8, width = 9.6, dpi = 600)

