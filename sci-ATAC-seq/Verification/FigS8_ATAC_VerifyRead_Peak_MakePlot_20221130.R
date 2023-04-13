library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(plyr)
library(ggvenn)
library(GenomicRanges)
library(VennDiagram) 

## Set directory
out.directory <- "/home/gayan/Projects/scATAC_Simulator/results/20221130_BoneMarrow_62016_chr1_NONINPUT"
fig.directory <- paste0(out.directory, "/", "VerifyReadsPlots")
filename <- "BoneMarrow_62016_chr1"
dir.create(fig.directory)

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


############# Plot 1 Venn Plot ############# 
real_peakfile <- sprintf("%s/BoneMarrow_62016_chr1.MACS3_peaks.bed", out.directory)
# synthetic_peakfile <- sprintf("%s/BoneMarrow_62016_chr1.syntheticBAM.CBincluded.synthetic.chr1.sorted.bam.MACS3.lessStringent_peaks.narrowPeak", out.directory)
# synthetic_peakfile <- sprintf("%s/BoneMarrow_62016_chr1.syntheticBAM.CBincluded.synthetic.chr1.sorted.bam.MACS3_peaks.narrowPeak", out.directory)
synthetic_peakfile <- sprintf("%s/BoneMarrow_62016_chr1.syntheticBAM.CBincluded.synthetic.chr1.sorted.bam.MACS3.stringent01_peaks.narrowPeak", out.directory)

real_peak <- read.table(real_peakfile, sep="\t")
colnames(real_peak) <- c("chr", "start","end")
synthetic_peak_MACS3 <- read.table(synthetic_peakfile, sep="\t")
colnames(synthetic_peak_MACS3)[1:3] <- c("chr", "start","end")
# discovery_length <- synthetic_peak_MACS3[,3] - synthetic_peak_MACS3[,2]
# overlap_size <- min(discovery_length)/2
discovery_length <- real_peak[,3] - real_peak[,2]
overlap_size <- min(discovery_length)/2
# overlap_size <- 0

# Compare
real_peak_GR <- makeGRangesFromDataFrame(as.data.frame(real_peak))
synthetic_peak_GR <- makeGRangesFromDataFrame(as.data.frame(synthetic_peak_MACS3))
TP <- sum(countOverlaps(synthetic_peak_GR, real_peak_GR, minoverlap = overlap_size) >0)
FP <- nrow(synthetic_peak_MACS3) - TP
FN <- nrow(real_peak) - TP

# Prepare data
save(synthetic_peak_MACS3, real_peak, TP, file=sprintf("%s/VennPlot_PeakCalling_plot_20230102.Rdata", fig.directory))

load(sprintf("%s/VennPlot_PeakCalling_plot_20230102.Rdata", fig.directory))

method_list <- c("scReadSim", "Real")
venn_plot <- draw.pairwise.venn( area1 = nrow(synthetic_peak_MACS3), area2 = nrow(real_peak), cross.area = TP, category = c("scReadSim", "Real"), fill = c("#00BFC4", "#F8766D"), cat.col= c("#00BFC4", "#F8766D"), fontfamily = rep("Helvetica",3), cat.fontfamily = rep("Helvetica",2) ,lty = "blank", cex = 1, cat.cex = 1, 
                                 # cat.dist = 0.09, cat.just = list(c(-1,-1), c(1, 1)), 
                                 cat.pos = c(330, 30),ext.pos = 30, ext.dist =-0.05, ext.length = 0.85, ext.line.lwd = 2, ext.line.lty = "dashed" )

# venn_list_2component <- list()
# venn_list_2component$scReadSim <- c(seq(1,FP), -seq(1,TP))
# venn_list_2component$Real <- c(FP + seq(1,FN), -seq(1,TP))
# save(venn_list_2component, file=sprintf("%s/VennPlot_PeakCalling_plot.Rdata", fig.directory))
# venn_plot <- ggvenn(venn_list_2component,  stroke_size = 0.5,text_size = 3,
#                     set_name_color = c("#00BFC4", "#F8766D"), set_name_size = 6,
#                     fill_color = c("#00BFC4", "#F8766D"))

ggsave(venn_plot, file=sprintf("%s/VennPlot_PeakCalling_plot_20230102.pdf", fig.directory))


############# Plot 2 FPR & FDR Plot ############# 
real_peakfile <- sprintf("%s/BoneMarrow_62016_chr1.MACS3_peaks.bed", out.directory)
synthetic_peakfile <- sprintf("%s/BoneMarrow_62016_chr1.syntheticBAM.CBincluded.synthetic.chr1.sorted.bam.MACS3.lessStringent_peaks.narrowPeak", out.directory)
# synthetic_peakfile <- sprintf("%s/BoneMarrow_62016_chr1.syntheticBAM.CBincluded.synthetic.chr1.sorted.bam.MACS3_peaks.narrowPeak", out.directory)

real_peak <- read.table(real_peakfile, sep="\t")
colnames(real_peak) <- c("chr", "start","end")
synthetic_peak_MACS3 <- read.table(synthetic_peakfile, sep="\t")
colnames(synthetic_peak_MACS3)[1:3] <- c("chr", "start","end")
# Sort by p-val
synthetic_peak_MACS3_sorted <- synthetic_peak_MACS3[order(synthetic_peak_MACS3[,8], decreasing = TRUE),]
# discovery_length <- synthetic_peak_MACS3_sorted[,3] - synthetic_peak_MACS3_sorted[,2]
# overlap_size <- min(discovery_length)/2
discovery_length <- real_peak[,3] - real_peak[,2]
overlap_size <- min(discovery_length)/2

MACS3_roc_list <- rocPeakCalculator(synthetic_peak_MACS3_sorted, real_peak, overlap_size)
MACS3_FDR_list <- fdrPeakCalculator(synthetic_peak_MACS3_sorted, real_peak, overlap_size)
save(MACS3_roc_list, MACS3_FDR_list, file=sprintf("%s/FDR_ROC_PeakCalling_plot_20221130.Rdata", fig.directory))

load(sprintf("%s/FDR_ROC_PeakCalling_plot_20221130.Rdata", fig.directory))
PrecisionRecall_plot <- as.data.frame(MACS3_FDR_list) %>%
  ggplot( aes(x=TPR_vec, y=1-FDR_vec)) +
  geom_line(size=0.5, color="#FDB462") +
  ylab("Precision") +
  xlab("Recall") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0, 1.0)) +
  scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0, 1.0)) +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_blank(), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        strip.text = element_text(size = 15))
ggsave(PrecisionRecall_plot, file= sprintf("%s/PrecisionRecall_lessStringent_20221130.pdf", fig.directory))

roc_plot <- as.data.frame(MACS3_roc_list) %>%
  ggplot( aes(x=FPR_vec, y=TPR_vec)) +
  geom_line(size=0.5, color="#FDB462") +
  xlab("False Positive Rate (FPR)") +
  ylab("True Positive Rate (TPR)") +
  scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0, 1.0)) +
  scale_x_continuous(breaks = seq(0,1,0.2), limits = c(0, 1.0)) +
  theme_bw() +
  theme(legend.position="none",
        legend.title = element_blank(), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.text=element_text(size=12),
        strip.text = element_text(size = 15))
ggsave(roc_plot, file= sprintf("%s/ROC_lessStringent_20221130.pdf", fig.directory))
library(ggpubr)
curve_plot <- ggarrange(roc_plot, PrecisionRecall_plot, ncol=2, common.legend = TRUE, legend="bottom")
curve_plot <- annotate_figure(curve_plot, top = text_grob("Consistency of MACS3's peak calling results\n on real and synthetic data", 
                                                          size = 15))                                                          
#############  Combine Figures  ############# 
gs <- list(arrangeGrob(grobTree(venn_plot), left = textGrob("A", x = unit(1, "npc"), 
                                                            y = unit(.95, "npc"))),
           arrangeGrob(curve_plot, left = textGrob("B", x = unit(1, "npc"), 
                                                   y = unit(.95, "npc"))))
lay <- rbind(c(1,1,2,2,2,2),
             c(1,1,2,2,2,2))
p_final <-arrangeGrob(grobs = gs, layout_matrix = lay)
ggsave(file=sprintf("%s/%s.VeryRead.Peak.Combined_20230102.pdf",fig.directory, filename), p_final, height =3.2 , width = 9.6, dpi = 600)


