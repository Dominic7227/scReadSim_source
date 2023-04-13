library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(plyr)

## Set directory
out.directory <- "/home/gayan/Projects/scATAC_Simulator/results/20221130_BoneMarrow_62016_chr1_NONINPUT"
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
    geom_line(data = kmer_count_df, aes(x = frequency, y = multiplicity, color = method, group = method), size=1.2, alpha= 0.1) +
    xlab(sprintf("%s-mer frequency x", kmer_length)) +
    ylab(sprintf("# %s-mer with frequency x", kmer_length)) +
    facet_wrap(~group) +
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
ggsave(file=sprintf("%s/%s.VeryRead.Combined_20230128.pdf",fig.directory, filename), p_final, height = 8, width = 9.6, dpi = 600)

