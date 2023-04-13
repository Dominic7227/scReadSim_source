library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(plyr)

## Set directory
out.directory <- "/home/gayan/Projects/scATAC_Simulator/results/20221011_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT"
# out.directory <- "/Users/gayan/Dropbox/Research/ATACseq_peak_calling/Fragment Simulator/Results/20220407_e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1_NONINPUT"
fig.directory <- paste0(out.directory, "/", "VerifyReadsPlots")
filename <- "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1"
dir.create(fig.directory)

############# Plot 1 kmer frequencies ############# 
for (kmer_length in c(11,21,31,41)){
  kmer_count_real_read <- read.table(sprintf("%s/VerifyCoverage/kmer_real/mer%s_counts_read_hist.txt", out.directory, kmer_length), header = FALSE)
  kmer_count_scReadSim_read <- read.table(sprintf("%s/VerifyCoverage/kmer_scReadSim/mer%s_counts_read2_hist.txt", out.directory, kmer_length), header = FALSE)
  colnames(kmer_count_real_read) <- c("frequency", "multiplicity")
  colnames(kmer_count_scReadSim_read) <- c("frequency", "multiplicity")
  
  kmer_count_real_df <- kmer_count_real_read
  # kmer_count_real_df$group <- c(rep("Read 1", nrow(kmer_count_real_read1)), rep("Read 2", nrow(kmer_count_real_read2)))
  kmer_count_scReadSim_df <- kmer_count_scReadSim_read
  # kmer_count_scReadSim_df$group <- c(rep("Read 1", nrow(kmer_count_scReadSim_read1)), rep("Read 2", nrow(kmer_count_scReadSim_read2)))
  
  kmer_count_df <- rbind(kmer_count_real_df, kmer_count_scReadSim_df)
  kmer_count_df$method <- c(rep("Real", nrow(kmer_count_real_df)), c(rep("scReadSim", nrow(kmer_count_scReadSim_df))))
  assign(sprintf("mer%s_count_df", kmer_length), kmer_count_df)
  save(kmer_count_df, file=sprintf("%s/Plot_%smer.Rdata", fig.directory, kmer_length))
  
  p_kmer <- ggplot() +
    geom_line(data = kmer_count_df, aes(x = frequency, y = multiplicity, color = method, group = method), size=1.2, alpha= 0.5) +
    xlab(sprintf("%s-mer frequency x", kmer_length)) +
    ylab(sprintf("# unique %s-mer with frequency x", kmer_length)) +
    # facet_wrap(~group) +
    scale_y_log10() +
    theme_bw() + 
    theme(legend.position="right",
          legend.title = element_blank(), 
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          legend.text=element_text(size=12),
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
  ggrastr::rasterize(geom_point(data = combined_mer_df, position=position_jitter(h=0.15,w=0.15), aes(x = frequency, y = multiplicity, color = method, group = method), size=0.5, alpha= 1),dpi = 300) +
  # xlab("k-mer frequency x") +
  # ylab("# unique k-mers with frequency x") +
  labs(x=expression(paste("k-mer frequency ",italic("x"))),
       y=expression(paste("# unique k-mers with frequency ",italic("x")))) +
  facet_wrap(~ case, scales="free_x",ncol = 2) +
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



#############  Plot 2 Error rate ############# 
error_rate_real <- read.table(sprintf("%s/VerifyQuality/fgbio/Real.error_rate_by_read_position.txt", out.directory), header = TRUE)
# error_rate_real <- read.table(sprintf("%s/VerifyQuality/fgbio/Real.error_rate_by_read_position.txt", out.directory), header = TRUE)
# error_rate_scReadSim <- read.table(sprintf("%s/VerifyQuality/fgbio/scReadSim_manualError.error_rate_by_read_position.txt", out.directory), header = TRUE)
error_rate_scReadSim <- read.table(sprintf("%s/VerifyQuality/fgbio/scReadSim_manualError.error_rate_by_read_position.txt", out.directory), header = TRUE)

error_rate_real_df <- error_rate_real %>% dplyr::select(c(position, error_rate))
error_rate_scReadSim_df <- error_rate_scReadSim %>% dplyr::filter(read_number == 0) %>% dplyr::select(c(position, error_rate))

error_rate_df <- rbind(error_rate_real_df, error_rate_scReadSim_df)
error_rate_df$method <- c(rep("Real", nrow(error_rate_real_df)), c(rep("scReadSim", nrow(error_rate_scReadSim_df))))

save(error_rate_df, file=sprintf("%s/Errorrate_plot.Rdata", fig.directory))

load(file=sprintf("%s/Errorrate_plot.Rdata", fig.directory))
p_errorrate <- ggplot() +
  geom_line(data = error_rate_df %>% filter(method=="Real"), aes(x = position, y = error_rate, color = method, group = method), size=1.5, alpha=0.5) +
  geom_line(data = error_rate_df %>% filter(method=="scReadSim"), aes(x = position, y = error_rate, color = method, group = method), size=0.7, alpha=0.5) +
  scale_x_continuous(limits = c(0, 90)) +
  xlab(sprintf("Position in Read")) +
  ylab("Error rate") +
  # facet_wrap(~group, scales = "free_x") +
  scale_y_log10() +
  theme_bw() + 
  theme(legend.position="bottom",
        legend.title = element_blank(), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        strip.text = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))
ggsave(plot=p_errorrate, sprintf("%s/Errorate_plot_new.pdf", fig.directory), width=10, height = 5)


#############  Combine Figures  ############# 
gs <- list(arrangeGrob(p_kmercombined, left = textGrob("A", x = unit(1, "npc"), 
                                                       y = unit(.95, "npc"))),
           arrangeGrob(p_errorrate, left = textGrob("B", x = unit(1, "npc"), 
                                                    y = unit(.95, "npc"))))
# lay <- rbind(c(1,1,1,1),
#              c(1,1,1,1),
#              c(2,2,2,2),
#              c(2,2,2,2))
lay <- rbind(c(1,1,2,2),
             c(1,1,2,2))
p_final <-arrangeGrob(grobs = gs, layout_matrix = lay)
ggsave(file=sprintf("%s/%s.VerifyRead.Combined_20221021.pdf",fig.directory, filename), p_final, height = 4.8, width = 9.6, dpi = 600)
