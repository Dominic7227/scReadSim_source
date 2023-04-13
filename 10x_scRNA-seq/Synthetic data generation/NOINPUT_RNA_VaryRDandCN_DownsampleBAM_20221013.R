library(dplyr)

outdirectory <- "/home/gayan/Projects/scATAC_Simulator/results/20221013_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_VaryNewCell"


n_downsample <- 1219
adj_factor_vec <- c(0.5, 1, 2, 4)
for (adj_factor in adj_factor_vec) {
	set.seed(2022 + adj_factor*100)
	cellbarcode_input <- sprintf("%s/synthetic_cell_barcode.VaryCellNumber%s.txt", outdirectory, adj_factor)
	cell_barcode_list <- as.data.frame(read.table(cellbarcode_input, sep="\t", header=FALSE))
	sampled_ind <- sample(seq(nrow(cell_barcode_list)), size = n_downsample, replace = FALSE)
	sample_cb <- cell_barcode_list[sampled_ind,]
	write.table(sample_cb, sprintf("%s/synthetic_cell_barcode.VaryCellNumber%s.Downsample1219.txt", outdirectory, adj_factor), sep="\t", row.names = FALSE,col.names = FALSE,quote = FALSE)
}
