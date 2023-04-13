#!/usr/bin/env Rscript
# if (suppressMessages(!require("UpSetR"))) suppressMessages(install.packages("UpSetR", repos="http://cran.us.r-project.org"))
# Changed the original make_main_bar funciton: comment out the log-transformation for the frequency. Remove the ylim option for ggplot.
# devtools::install_local("/Users/gayan/Dropbox/Mac/Downloads/UpSetR-master.zip", force = TRUE)
library("UpSetR")
pdf("/Users/gayan/Dropbox/Research/ATACseq_peak_calling/Fragment Simulator/Results/20221130_BoneMarrow_62016_chr1_INPUT_withCluster_mm9TSS/INTERVENE_20221130_noSEACR/Intervene_upset_log_noSEACR.pdf", width=9.6, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('HMMRATAC'=8,'HOMER'=1962,'HOMER&HMMRATAC'=5,'MACS3'=495,'MACS3&HMMRATAC'=2,'MACS3&HOMER'=264,'MACS3&HOMER&HMMRATAC'=2,'INPUT'=70,'INPUT&HMMRATAC'=40,'INPUT&HOMER'=77,'INPUT&HOMER&HMMRATAC'=96,'INPUT&MACS3'=21,'INPUT&MACS3&HMMRATAC'=6,'INPUT&MACS3&HOMER'=234,'INPUT&MACS3&HOMER&HMMRATAC'=2696)
upset(fromExpression(expressionInput), nsets=5, nintersects=30, show.numbers="yes", main.bar.color="#ea5d4e", sets.bar.color="#317eab", empty.intersections=NULL, order.by = "freq", number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size", scale.intersections = "log2", scale.sets = "log2")
invisible(dev.off())


pdf("/Users/gayan/Dropbox/Research/ATACseq_peak_calling/Fragment Simulator/Results/20221130_BoneMarrow_62016_chr1_INPUT_withCluster_mm9TSS/INTERVENE_20221130_withSEACR/Intervene_upset_log_withSEACR.pdf", width=9.6, onefile=FALSE, useDingbats=FALSE)
expressionInput <- c('HMMRATAC'=0,'HOMER'=8,'HOMER&HMMRATAC'=0,'SEACR'=29703,'SEACR&HMMRATAC'=28,'SEACR&HOMER'=1743,'SEACR&HOMER&HMMRATAC'=0,'MACS3'=56,'MACS3&HMMRATAC'=0,'MACS3&HOMER'=0,'MACS3&HOMER&HMMRATAC'=0,'MACS3&SEACR'=439,'MACS3&SEACR&HMMRATAC'=2,'MACS3&SEACR&HOMER'=264,'MACS3&SEACR&HOMER&HMMRATAC'=2,'INPUT'=0,'INPUT&HMMRATAC'=0,'INPUT&HOMER'=0,'INPUT&HOMER&HMMRATAC'=0,'INPUT&SEACR'=70,'INPUT&SEACR&HMMRATAC'=40,'INPUT&SEACR&HOMER'=77,'INPUT&SEACR&HOMER&HMMRATAC'=96,'INPUT&MACS3'=0,'INPUT&MACS3&HMMRATAC'=0,'INPUT&MACS3&HOMER'=0,'INPUT&MACS3&HOMER&HMMRATAC'=0,'INPUT&MACS3&SEACR'=21,'INPUT&MACS3&SEACR&HMMRATAC'=6,'INPUT&MACS3&SEACR&HOMER'=234,'INPUT&MACS3&SEACR&HOMER&HMMRATAC'=2696)
upset(fromExpression(expressionInput), nsets=5, nintersects=30, show.numbers="yes", main.bar.color="#ea5d4e", sets.bar.color="#317eab", empty.intersections=NULL, order.by = "freq", number.angles = 0, mainbar.y.label ="No. of Intersections", sets.x.label ="Set size", scale.intersections = "log2", scale.sets = "log2")
invisible(dev.off())
