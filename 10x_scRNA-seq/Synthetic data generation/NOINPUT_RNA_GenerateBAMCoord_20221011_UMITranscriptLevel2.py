# Create fastq files as simulated reads output
# Start from half-sampled reads coordinates
#BAMfile_halfsampled_coordinates.txt

## JOBS:
## For each referenced peak region, find shift length from true peak
## Write mates for each read, and fill in read end position
## Separate read1 and read2

## Input true and ref peak pos
## 1. Shift the read positions with shift_number
## 2. Find all reads mate
### Load libraries 
import pandas as pd
# import pickle
import numpy as np
import csv
import collections
import time
import sys
pd.options.mode.chained_assignment = None  # default='warn'
import string
import random
from tqdm import tqdm
import pysam
from collections import defaultdict

def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

def cellbarcode_generator(length, size=10):
	chars = 'ACGT'
	cb_list = [''.join(random.choice(chars) for _ in range(size)) for cell in range(length)]
	return cb_list


def PerTruePeakEdition(peak_record, count_vec, read_lines, read_len, jitter_size, random_cellbarcode_list):
	# peak_record = peaks_assignments.loc[1,] # Input
	true_peak_concat = peak_record[0] + ":" + str(peak_record[1]) + "-" + str(peak_record[2])
	reads_cur = read_lines[read_lines['peak_name'] == true_peak_concat] # Input
	nread_cur= np.sum(count_vec).astype(int)
	# Add cell information
	nonempty_cell_ind = np.where(count_vec != 0)[0]
	random_umi_list = cellbarcode_generator(nread_cur, size=10)
	read_code_simu_cur = [random_cellbarcode_list[nonempty_cell_ind[ind]] + ":CellType1" + "CellNo" + str(nonempty_cell_ind[ind] + 1) + ":" + str(true_peak_concat) + "#" + str(count).zfill(4) for ind in range(len(nonempty_cell_ind)) for count in range(count_vec[nonempty_cell_ind[ind]])]
	read_code_simu_cur_withUMI = [read_code_simu_cur[read_id].split(":",1)[0] + random_umi_list[read_id] + ":" + read_code_simu_cur[read_id].split(":",1)[1] for read_id in range(len(read_code_simu_cur))]
# read_code_simu_cur = ["CellType1" + "CellNo" + str(nonempty_cell_ind[ind] + 1) + ":" + str(true_peak_concat) + "#" + str(count).zfill(4) for ind in range(len(nonempty_cell_ind)) for count in range(count_vec[nonempty_cell_ind[ind]])]
	# start = time.time()
	jitter_value_vec = np.random.random_integers(-jitter_size,jitter_size,size=np.shape(reads_cur)[0])  # nrow(reads_cur) should equal to nfrag_cur
	if sum(reads_cur['r1_start'].isna()) == 0:
		reads_cur['r1_start_shifted'] = reads_cur['r1_start'].astype(int) + jitter_value_vec
		reads_cur['r1_end_shifted'] = reads_cur['r1_start'].astype(int) + read_len + jitter_value_vec
		read_1_df = reads_cur[['chr','r1_start_shifted', 'r1_end_shifted']]
		read_1_df['read_name'] = read_code_simu_cur_withUMI
		read_1_df['read_length'] = read_len
		read_1_df['strand'] = '+'
	else:
		read_1_df = pd.DataFrame({'A' : []}) # Create an empty dataframe
	return read_1_df

def PerTruePeakEdition_UMIandRead(peak_record, count_vec, UMI_count_vec, read_lines, read_len, jitter_size, random_cellbarcode_list):
	# peak_record = peaks_assignments.loc[1,] # Input
	true_peak_concat = peak_record[0] + ":" + str(peak_record[1]) + "-" + str(peak_record[2])
	reads_cur = read_lines[read_lines['peak_name'] == true_peak_concat] # Input
	nread_cur= np.sum(count_vec).astype(int)
	# Add cell information
	nonempty_cell_ind = np.where(count_vec != 0)[0]
	# Use UMI_count_vec to guide UMI generation
	random_umi_list = [random.sample(cellbarcode_generator(UMI_count_vec[nonempty_cell_ind_cur], size=10)*count_vec[nonempty_cell_ind_cur], count_vec[nonempty_cell_ind_cur]) if UMI_count_vec[nonempty_cell_ind_cur] > 0 else cellbarcode_generator(1, size=10)*count_vec[nonempty_cell_ind_cur] for nonempty_cell_ind_cur in nonempty_cell_ind]
	random_umi_list = [x for xs in random_umi_list for x in xs]
	# random_umi_list = cellbarcode_generator(nread_cur, size=10)
	read_code_simu_cur = [random_cellbarcode_list[nonempty_cell_ind[ind]] + ":" + "CellNo" + str(nonempty_cell_ind[ind] + 1) + ":" + str(true_peak_concat) + "#" + str(count).zfill(4) for ind in range(len(nonempty_cell_ind)) for count in range(count_vec[nonempty_cell_ind[ind]])]
	read_code_simu_cur_withUMI = [read_code_simu_cur[read_id].split(":",1)[0]  + random_umi_list[read_id] + ":" + read_code_simu_cur[read_id].split(":",1)[1] for read_id in range(len(read_code_simu_cur))]
	# read_code_simu_cur = ["CellType1" + "CellNo" + str(nonempty_cell_ind[ind] + 1) + ":" + str(true_peak_concat) + "#" + str(count).zfill(4) for ind in range(len(nonempty_cell_ind)) for count in range(count_vec[nonempty_cell_ind[ind]])]
	# start = time.time()
	jitter_value_vec = np.random.random_integers(-jitter_size,jitter_size,size=np.shape(reads_cur)[0])  # nrow(reads_cur) should equal to nfrag_cur
	if sum(reads_cur['r1_start'].isna()) == 0:
		reads_cur['r1_start_shifted'] = reads_cur['r1_start'].astype(int) + jitter_value_vec
		reads_cur['r1_end_shifted'] = reads_cur['r1_start'].astype(int) + read_len + jitter_value_vec
		read_1_df = reads_cur[['chr','r1_start_shifted', 'r1_end_shifted']]
		read_1_df['read_name'] = read_code_simu_cur_withUMI
		read_1_df['read_length'] = read_len
		read_1_df['strand'] = '+'
	else:
		read_1_df = pd.DataFrame({'A' : []}) # Create an empty dataframe
	return read_1_df


def UMI_read_dict(bed_file,  UMI_count_mat_file, read1_bedfile, INPUT_bamfile, outdirectory, OUTPUT_cells_barcode_file, jitter_size=5, read_len=90):
	UMI_count_mat_df = pd.read_csv("%s/%s" % (outdirectory, UMI_count_mat_file), header=0, delimiter="\t")
	UMI_count_mat = UMI_count_mat_df.to_numpy()
	UMI_count_mat_cluster = UMI_count_mat_df.columns.to_numpy()
	n_cell = np.shape(UMI_count_mat)[1]
	samfile = pysam.AlignmentFile(INPUT_bamfile, "rb")
	with open(bed_file) as open_peak:
	    reader = csv.reader(open_peak, delimiter="\t")
	    open_peak = np.asarray(list(reader))
	peak_nonzero_id = np.nonzero(UMI_count_mat.sum(axis=1))[0]
	random.seed(2022)
	random_cellbarcode_list = cellbarcode_generator(n_cell, size=16)
	with open(OUTPUT_cells_barcode_file, 'w') as f:
			for item in random_cellbarcode_list:
			    f.write(item + "\n")
	cellbarcode_list_withclusters = np.vstack((random_cellbarcode_list, UMI_count_mat_cluster)).transpose()
	with open(OUTPUT_cells_barcode_file + ".withSynthCluster", 'w') as f:
			for item in cellbarcode_list_withclusters:
			    f.write("\t".join(item) + "\n")
	with open("%s/%s" % (outdirectory, read1_bedfile), 'w') as fp:
		pass
	start = time.time()
	for relative_peak_ind in tqdm(range(len(peak_nonzero_id))):
		peak_ind = peak_nonzero_id[relative_peak_ind]
		# peak_ind = 1
		rec = open_peak[peak_ind]
		rec_name = '_'.join(rec)
		UMI_read_dict_pergene = defaultdict(lambda: [])
		reads = samfile.fetch(rec[0], int(rec[1]), int(rec[2]))
		# TODO: create a dictionary recording UMI(key) to reads starting positions(values)
		for read in reads:
			if read.has_tag('UB:Z'):
				UMI = read.get_tag('UB:Z')
				start = read.reference_start
				if read.is_reverse==1:
					strand = "-" 
				else:
					strand = "+"
				# strand = read.is_reverse
				UMI_read_dict_pergene[UMI].append(str(start) + ":" + str(strand))
		## Real data: Obtain the empirical distribution of # reads of UMI in the gene 
		UMI_real_vec = list(UMI_read_dict_pergene.keys())
		nread_perUMI = [len(x) for x in UMI_read_dict_pergene.values()] # return a vector with the length of # UMIs in the gene
		# nread_perUMI_prob = nread_perUMI / np.sum(nread_perUMI)
		## Synthetic data
		# Sample syntheitc UMI's read number and assign real UMI to synthetic UMI (only for sampling reads)
		UMI_count_vec = UMI_count_mat[peak_ind,:] # Synthetic umi count
		nonzer_UMI_zero_id = np.nonzero(UMI_count_vec)[0]
		synthetic_nread_perUMI_percell = [np.random.choice(nread_perUMI, size=UMI_count_vec[i], replace=True) for i in nonzer_UMI_zero_id] # return a vector with the length of non-zero count cells, each entry idnicates the number of reads for each UMI
		nread_synthetic = sum(np.sum(x) for x in synthetic_nread_perUMI_percell) # total number of synthetic reads
		# Assign real UMI to synthetic UMI based on read number per UMI
		# Construct a dic with keys as number of reasd per UMI and values as the corresponding UMI strings
		nread_UMI_dic = defaultdict(lambda: [])
		for i,item in enumerate(nread_perUMI):
			nread_UMI_dic[item].append(UMI_real_vec[i])
		synthetic_nread_perUMI_percell_unlist = [xs for x in synthetic_nread_perUMI_percell for xs in list(x)]
		synthetic_realUMI_assign_array = np.array(synthetic_nread_perUMI_percell_unlist).astype(str)
		synthetic_nread_perUMI_percell_unlist_array = np.array(synthetic_nread_perUMI_percell_unlist)
		for i in set(synthetic_nread_perUMI_percell_unlist):
			tmp_ind = np.where(synthetic_nread_perUMI_percell_unlist_array == i)[0]
			synthetic_realUMI_assign_array[tmp_ind] = np.random.choice(nread_UMI_dic[i], size=len(tmp_ind), replace=True)
		synthetic_realUMI_assign_vec = list(synthetic_realUMI_assign_array) # a vectir sane length to synthetic_nread_perUMI_percell_unlist, synthetic UMI's number 
		# synthetic_realUMI_assign_vec = [np.random.choice(nread_UMI_dic[nread], size=1, replace=True) for nread in synthetic_nread_perUMI_percell_unlist]
		# 
		# Sample reads for each synthetic UMI
		synthetic_read_list = [np.random.choice(UMI_read_dict_pergene[synthetic_realUMI_assign_vec[id]], size=synthetic_nread_perUMI_percell_unlist[id], replace=True) for id in range(len(synthetic_realUMI_assign_vec))]
		synthetic_read_unlist = [xs for x in synthetic_read_list for xs in list(x)]
		synthetic_read_start_list = [int(xs.split(':', 1)[0]) - 1  for xs in synthetic_read_unlist]  # -1 accounts for the BAM coordinates is +1 compared with getfasta bedtools
		synthetic_read_strand_list = [str(xs.split(':', 1)[1]) for xs in synthetic_read_unlist]
		# Creat synthetic UMI string
		UMI_synthetic_uniq = cellbarcode_generator(np.sum(UMI_count_vec), size=10)
		UMI_synthetic_list = np.repeat(np.array(UMI_synthetic_uniq), synthetic_nread_perUMI_percell_unlist)
		# Create syntheitc CB string
		CB_synthetic_repeat = [np.sum(x) for x in synthetic_nread_perUMI_percell] # same length to nozero cells
		CB_synthetic_list = np.repeat(np.array(random_cellbarcode_list)[nonzer_UMI_zero_id], CB_synthetic_repeat)
		# Create syntheitc read name string
		read_name_remaining_list = ["CellNo" + str(nonzer_UMI_zero_id[ind] + 1) + ":" + str(rec_name) + "#" + str(count).zfill(4) for ind in range(len(CB_synthetic_repeat)) for count in range(CB_synthetic_repeat[ind])]
		read_name_list = [CB_synthetic_list[id] + UMI_synthetic_list[id] + ":" + read_name_remaining_list[id] for id in range(nread_synthetic)]
		# Create output table
		jitter_value_vec = np.random.randint(-jitter_size,jitter_size,size=nread_synthetic)  # nrow(reads_cur) should equal to nfrag_cur
		reads_cur = pd.DataFrame({
			'chr': rec[0],
			'r1_start_shifted': synthetic_read_start_list  + jitter_value_vec,
			'r1_end_shifted': synthetic_read_start_list + jitter_value_vec + read_len,
			'read_name': read_name_list,
			'read_length': read_len,
			'strand': synthetic_read_strand_list
			})
		reads_cur.to_csv("%s/%s" % (outdirectory, read1_bedfile), header=None, index=None, sep='\t', mode='a')
	end = time.time()

def main():
	# random.seed(2022)
	user_args = sys.argv[1:]
	bed_file,  UMI_count_mat_file, read1_bedfile, INPUT_bamfile, outdirectory, OUTPUT_cells_barcode_file = user_args
	jitter_size = 5
	read_len = 90
	UMI_read_dict(bed_file,  UMI_count_mat_file, read1_bedfile, INPUT_bamfile, outdirectory, OUTPUT_cells_barcode_file, jitter_size, read_len)



if __name__ == '__main__':
  main()     


# # Analysis of intergenetic regions
# outdirectory = '/home/gayan/Projects/scATAC_Simulator/results/20221011_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT'
# INPUT_bamfile = "/home/gayan/Projects/scATAC_Simulator/data/10X_MOUSE_BRAIN_ATACandRNA/split.e18_mouse_brain_fresh_5k_gex_possorted_bam/e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.bam"
# OUTPUT_cells_barcode_file = outdirectory + "/synthetic_cell_barcode.txt"
# bed_file = "/home/gayan/Projects/scATAC_Simulator/results/20221011_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_peaks.COMPLE.bed"
# UMI_count_mat_file = "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.UMI.COMPLE.countmatrix.scDesign2Simulated.txt"
# read1_bedfile = "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.syntheticBAM.COMPLE.20221011_UMITransLevel2" + ".read.bed"

# # Analysis of genetic regions
# # outdirectory = '/home/gayan/Projects/scATAC_Simulator/results/20221011_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT'
# INPUT_bamfile = "/home/gayan/Projects/scATAC_Simulator/data/10X_MOUSE_BRAIN_ATACandRNA/split.e18_mouse_brain_fresh_5k_gex_possorted_bam/e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.bam"
# OUTPUT_cells_barcode_file = outdirectory + "/synthetic_cell_barcode.txt"
# bed_file = "/home/gayan/Projects/scATAC_Simulator/results/20221011_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_peaks.bed"
# UMI_count_mat_file = "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.UMI.countmatrix.scDesign2Simulated.txt"
# read1_bedfile = "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.syntheticBAM.20221011_UMITransLevel2" + ".read.bed"
