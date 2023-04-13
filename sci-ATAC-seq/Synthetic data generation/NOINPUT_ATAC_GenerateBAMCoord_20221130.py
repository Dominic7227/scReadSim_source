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
import pickle
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

def flatten(x):
    if isinstance(x, collections.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

def cellbarcode_generator(length, size=10):
	chars = 'ACGT'
	cb_list = [''.join(random.choice(chars) for _ in range(size)) for cell in range(length)]
	return cb_list

def find_leftnearest_nonpeak(non_peak_list, grey_peak):
	dist1 = np.abs(non_peak_list[:,2].astype(int) - int(grey_peak[1]))
	id = dist1.argmin()
	return id

# 20221130 TODO: 
# 1. Include parameter grey_area_modeling: when generate non-peak read, also generate grey-area
# 2. Distinguish peak and non-peaks for generating coordinates
def ATAC_GenerateBAMCoord(bed_file, count_mat_file, read_bedfile_prename, INPUT_bamfile, outdirectory, OUTPUT_cells_barcode_file, jitter_size=5, read_len=50, random_noise_mode = False, grey_area_modeling = True):
	count_mat_df = pd.read_csv("%s/%s" % (outdirectory, count_mat_file), header=0, delimiter="\t")
	count_mat = count_mat_df.to_numpy()
	count_mat_cluster = count_mat_df.columns.to_numpy()
	n_cell = np.shape(count_mat)[1]
	samfile = pysam.AlignmentFile(INPUT_bamfile, "rb")
	with open(bed_file) as file:
		reader = csv.reader(file, delimiter="\t")
		open_peak = np.asarray(list(reader))
	peak_nonzero_id = np.nonzero(count_mat.sum(axis=1))[0]
	random.seed(2022)
	random_cellbarcode_list = cellbarcode_generator(n_cell, size=16)
	with open(OUTPUT_cells_barcode_file, 'w') as f:
		for item in random_cellbarcode_list:
			f.write(item + "\n")
	cellbarcode_list_withclusters = np.vstack((random_cellbarcode_list, count_mat_cluster)).transpose()
	with open(OUTPUT_cells_barcode_file + ".withSynthCluster", 'w') as f:
		for item in cellbarcode_list_withclusters:
			f.write("\t".join(item) + "\n")
	# Create read 1 and read 2 files
	with open("%s/%s.read1.bed" % (outdirectory, read_bedfile_prename), 'w') as f_1:
		pass
	with open("%s/%s.read2.bed" % (outdirectory, read_bedfile_prename), 'w') as f_2:
		pass
	for relative_peak_ind in tqdm(range(len(peak_nonzero_id))):
		peak_ind = peak_nonzero_id[relative_peak_ind]
		rec = open_peak[peak_ind]
		rec_name = '_'.join(rec)
		reads = samfile.fetch(rec[0], int(rec[1]), int(rec[2]))
		# tic = time.time()
		# Extract read information
		reads_str = []
		for read in reads:
			if read.is_reverse==1:
				strand = -1
			else:
				strand = 1
			if read.is_read1==1:
				read_order = 1
			else:
				read_order = 2
			start = read.reference_start
			mate_start = read.next_reference_start
			read_len_cur = read.query_alignment_length
			read_info = [start, mate_start, read_len_cur, read_order, strand]
			reads_str.append(read_info)
		# toc = time.time()
		# toc - tic # 0.23621153831481934
		# print(read_info)
		# tic = time.time()
		# reads_str = [str(read).split("\t") for read in reads]
		# toc = time.time()
		# toc - tic # 0.19811201095581055
		# Sample npair_read_synthetic read 1 from reads and reserve fragment length 
		if len(reads_str) > 0: # If no real reads exist in the peak, skip
			count_vec = count_mat[peak_ind,:] # Synthetic umi count
			count_frag_vec = np.ceil(count_vec/2).astype(int)
			npair_read_synthetic = np.sum(count_frag_vec).astype(int) # nrow(reads_cur) should equal to nfrag_cur
			# npair_read_synthetic = np.ceil(np.sum(count_vec)/2).astype(int) # total number of synthetic reads
			read_sampled_unsplit = np.array(reads_str)[np.random.choice(len(reads_str), size=npair_read_synthetic, replace=True),:]
			# Sample starting position if random noise mode is on, or use real read starting position
			if random_noise_mode == True:
				read_synthetic_start = np.random.randint(int(rec[1]), int(rec[2]), size=npair_read_synthetic)
			else:
				read_synthetic_start = read_sampled_unsplit[:,0].astype(int)
			# Generate Read Name
			nonempty_cell_ind = np.where(count_frag_vec != 0)[0]
			target_peak_concat = rec[0] + ":" + str(rec[1]) + "-" + str(rec[2])
			read_name_list = [random_cellbarcode_list[nonempty_cell_ind[ind]] + ":" + "CellNo" + str(nonempty_cell_ind[ind] + 1) + ":" + str(target_peak_concat) + "#" + str(count).zfill(4) for ind in range(len(nonempty_cell_ind)) for count in range(count_frag_vec[nonempty_cell_ind[ind]])]
			# Create dataframe for unspiltted sampled reads
			reads_cur = pd.DataFrame({
				'chr': rec[0],
				'r_start': read_synthetic_start,
				'mate_start': read_sampled_unsplit[:,1].astype(int) + read_synthetic_start - read_sampled_unsplit[:,0].astype(int),
				'read_name': read_name_list,
				'length': read_sampled_unsplit[:,2],
				'read_order': read_sampled_unsplit[:,3].astype(int),
				'strand': read_sampled_unsplit[:,4].astype(int),
				'mate_strand': -read_sampled_unsplit[:,4].astype(int)
				})
			contain_read_indicator = read_sampled_unsplit[:,0] == read_sampled_unsplit[:,1]
			reads_cur['read_length'] = read_len
			reads_cur['read_length'][contain_read_indicator] = abs(reads_cur['length'].astype(int)[contain_read_indicator])
			# Add jitter size to read positions
			jitter_value_vec = np.random.randint(-jitter_size,jitter_size,size=npair_read_synthetic).astype(int)  # nrow(reads_cur) should equal to nfrag_cur
			reads_cur['r_start_shifted'] = reads_cur['r_start'].astype(int)  + jitter_value_vec
			reads_cur['mate_start_shifted'] = reads_cur['mate_start'].astype(int)  + jitter_value_vec
			reads_cur['r_end_shifted'] = reads_cur['r_start'].astype(int) + reads_cur['read_length'].astype(int) + jitter_value_vec
			reads_cur['mate_end_shifted'] = reads_cur['mate_start'].astype(int) + reads_cur['read_length'].astype(int) + jitter_value_vec
			# Split read 1 and read 2
			read_1_df = pd.concat([reads_cur.loc[reads_cur['read_order'] == 1, ['chr','r_start_shifted', 'r_end_shifted', 'read_length', 'strand']].rename(columns={'r_start_shifted':'r1_start_shifted', 'r_end_shifted':'r1_end_shifted'}), reads_cur.loc[reads_cur['read_order'] == 2, ['chr','mate_start_shifted', 'mate_end_shifted', 'read_length', 'mate_strand']].rename(columns={'mate_start_shifted':'r1_start_shifted', 'mate_end_shifted':'r1_end_shifted', 'mate_strand': 'strand'})], ignore_index=True)
			read_2_df = pd.concat([reads_cur.loc[reads_cur['read_order'] == 1, ['chr','mate_start_shifted', 'mate_end_shifted', 'read_length', 'mate_strand']].rename(columns={'mate_start_shifted':'r2_start_shifted', 'mate_end_shifted':'r2_end_shifted', 'mate_strand': 'strand'}), reads_cur.loc[reads_cur['read_order'] == 2, ['chr','r_start_shifted', 'r_end_shifted', 'read_length', 'strand']].rename(columns={'r_start_shifted':'r2_start_shifted', 'r_end_shifted':'r2_end_shifted'}), ], ignore_index=True)
			read_1_df['read_name'] = read_name_list
			read_2_df['read_name'] = read_name_list
			# read_1_df['read_length'] = read_len
			# read_2_df['read_length'] = read_len
			read_1_df['strand'] = ['+' if x == 1  else '-' for x in read_1_df['strand']]
			read_2_df['strand'] = ['+' if x == 1  else '-' for x in read_2_df['strand']]
			read_1_df_order = read_1_df[['chr','r1_start_shifted', 'r1_end_shifted', 'read_name', 'read_length', 'strand']]
			read_2_df_order = read_2_df[['chr','r2_start_shifted', 'r2_end_shifted', 'read_name', 'read_length', 'strand']]
			if read_1_df_order.shape[0] != read_2_df_order.shape[0]:
				print("[Warning] Peak %s read 1 and read 2 not identical!", relative_peak_ind)
			if np.sum(np.array(read_1_df_order[['r1_start_shifted']] < 0)) + np.sum(np.array(read_2_df_order[['r2_start_shifted']] < 0)) > 0:
				print("[Warning] Synthetic read pair for Peak %s %s has read 1 or read 2 start position negative: synthetic read pair removed!" % (relative_peak_ind, rec_name))
				ind_preserve = np.array(read_1_df_order[['r1_start_shifted']] >= 0) * np.array(read_2_df_order[['r2_start_shifted']] >= 0) # Remove rows with read 1 or read 2's start position negative
				read_1_df_order_removeNegRead = read_1_df_order.loc[ind_preserve]
				read_2_df_order_removeNegRead = read_2_df_order.loc[ind_preserve]
			else:
				read_1_df_order_removeNegRead = read_1_df_order
				read_2_df_order_removeNegRead = read_2_df_order
			read_1_df_order_removeNegRead.to_csv("%s/%s.read1.bed" % (outdirectory, read_bedfile_prename), header=None, index=None, sep='\t', mode='a')
			read_2_df_order_removeNegRead.to_csv("%s/%s.read2.bed" % (outdirectory, read_bedfile_prename), header=None, index=None, sep='\t', mode='a')
			if npair_read_synthetic != read_1_df_order_removeNegRead.shape[0]:
				print("Target read pair %s | Sample synthetic read pair %s" % (npair_read_synthetic, read_1_df_order_removeNegRead.shape[0]))
	## TODO Modeling grey_area
	if grey_area_modeling == True:
		print("Generating reads for Grey Area...")
		with open(outdirectory + "/" + "BoneMarrow_62016_chr1.MACS3_GreyAreas.bed") as file:
			reader = csv.reader(file, delimiter="\t")
			GreyArea_set = np.asarray(list(reader))
		with open("%s/GreyArea_Assigned_Synthetic_CountMatrix.txt" % outdirectory, 'w') as f_2:
			pass
		# Create read 1 and read 2 files
		with open("%s/%s.GreyArea.read1.bed" % (outdirectory, read_bedfile_prename), 'w') as f_1:
			pass
		with open("%s/%s.GreyArea.read2.bed" % (outdirectory, read_bedfile_prename), 'w') as f_2:
			pass
		random.seed(2022)
		for peak_id in tqdm(range(len(GreyArea_set))):
			grey_area = GreyArea_set[peak_id]
			rec_name = '_'.join(grey_area)
			grey_length = int(grey_area[2]) - int(grey_area[1])
			idx = find_leftnearest_nonpeak(open_peak, grey_area)
			nonpeak_cur_count = count_mat[idx,:]
			nonpeak_cur_length = int(open_peak[idx,2]) - int(open_peak[idx,1])
			if np.sum(nonpeak_cur_count) == 0:
				grey_count_vec = pd.DataFrame(np.zeros(n_cell, dtype=int)).T
				grey_count_vec.to_csv("%s/GreyArea_Assigned_Synthetic_CountMatrix.txt" % outdirectory, header=None, index=None, sep='\t', mode='a')
				continue # print zero counts for grey area if no reads in non-peak count mat
			# scaled_grey_count = np.round(nonpeak_cur_count * grey_length / nonpeak_cur_length).astype(int)
			scaled_grey_count = nonpeak_cur_count * np.random.binomial(1, min(grey_length / nonpeak_cur_length, 1), len(nonpeak_cur_count))
			# Write out grey area synthetic count matrix (optional)
			grey_count_vec = pd.DataFrame(np.append(rec_name, scaled_grey_count)).T
			grey_count_vec.to_csv("%s/GreyArea_Assigned_Synthetic_CountMatrix.txt" % outdirectory, header=None, index=None, sep='\t', mode='a')
			if np.sum(scaled_grey_count) == 0:
				continue # if no synthetic count for grey then skip the peak 
			reads = samfile.fetch(grey_area[0], int(grey_area[1]), int(grey_area[2]))
			# Extract read information
			reads_str = []
			for read in reads:
				if read.is_reverse==1:
					strand = -1
				else:
					strand = 1
				if read.is_read1==1:
					read_order = 1
				else:
					read_order = 2
				start = read.reference_start
				mate_start = read.next_reference_start
				read_len_cur = read.query_alignment_length
				read_info = [start, mate_start, read_len_cur, read_order, strand]
				reads_str.append(read_info)
			if len(reads_str) == 0: # If no real reads exist in the peak, skip
				continue
			count_frag_vec = np.ceil(scaled_grey_count/2).astype(int)
			npair_read_synthetic = np.sum(count_frag_vec).astype(int) # nrow(reads_cur) should equal to nfrag_cur
			read_sampled_unsplit = np.array(reads_str)[np.random.choice(len(reads_str), size=npair_read_synthetic, replace=True),:]
			# Sample starting position if random noise mode is on, or use real read starting position
			if random_noise_mode == True:
				read_synthetic_start = np.random.randint(int(grey_area[1]), int(grey_area[2]), size=npair_read_synthetic)
			else:
				read_synthetic_start = read_sampled_unsplit[:,0].astype(int)
			# Generate Read Name
			nonempty_cell_ind = np.where(count_frag_vec != 0)[0]
			target_peak_concat = grey_area[0] + ":" + str(grey_area[1]) + "-" + str(grey_area[2])
			read_name_list = [random_cellbarcode_list[nonempty_cell_ind[ind]] + "CellNo" + str(nonempty_cell_ind[ind] + 1) + ":" + str(target_peak_concat) + "#" + str(count).zfill(4) for ind in range(len(nonempty_cell_ind)) for count in range(count_frag_vec[nonempty_cell_ind[ind]])]
			# Create dataframe for unspiltted sampled reads
			reads_cur = pd.DataFrame({
				'chr': grey_area[0],
				'r_start': read_synthetic_start,
				'mate_start': read_sampled_unsplit[:,1].astype(int) + read_synthetic_start - read_sampled_unsplit[:,0].astype(int),
				'read_name': read_name_list,
				'length': read_sampled_unsplit[:,2],
				'read_order': read_sampled_unsplit[:,3].astype(int),
				'strand': read_sampled_unsplit[:,4].astype(int),
				'mate_strand': -read_sampled_unsplit[:,4].astype(int)
				})
			contain_read_indicator = read_sampled_unsplit[:,0] == read_sampled_unsplit[:,1]
			reads_cur['read_length'] = read_len
			reads_cur['read_length'][contain_read_indicator] = abs(reads_cur['length'].astype(int)[contain_read_indicator])
			# Add jitter size to read positions
			jitter_value_vec = np.random.randint(-jitter_size,jitter_size,size=npair_read_synthetic).astype(int)  # nrow(reads_cur) should equal to nfrag_cur
			reads_cur['r_start_shifted'] = reads_cur['r_start'].astype(int)  + jitter_value_vec
			reads_cur['mate_start_shifted'] = reads_cur['mate_start'].astype(int)  + jitter_value_vec
			reads_cur['r_end_shifted'] = reads_cur['r_start'].astype(int) + reads_cur['read_length'].astype(int) + jitter_value_vec
			reads_cur['mate_end_shifted'] = reads_cur['mate_start'].astype(int) + reads_cur['read_length'].astype(int) + jitter_value_vec
			# Split read 1 and read 2
			read_1_df = pd.concat([reads_cur.loc[reads_cur['read_order'] == 1, ['chr','r_start_shifted', 'r_end_shifted', 'read_length', 'strand']].rename(columns={'r_start_shifted':'r1_start_shifted', 'r_end_shifted':'r1_end_shifted'}), reads_cur.loc[reads_cur['read_order'] == 2, ['chr','mate_start_shifted', 'mate_end_shifted', 'read_length', 'mate_strand']].rename(columns={'mate_start_shifted':'r1_start_shifted', 'mate_end_shifted':'r1_end_shifted', 'mate_strand': 'strand'})], ignore_index=True)
			read_2_df = pd.concat([reads_cur.loc[reads_cur['read_order'] == 1, ['chr','mate_start_shifted', 'mate_end_shifted', 'read_length', 'mate_strand']].rename(columns={'mate_start_shifted':'r2_start_shifted', 'mate_end_shifted':'r2_end_shifted', 'mate_strand': 'strand'}), reads_cur.loc[reads_cur['read_order'] == 2, ['chr','r_start_shifted', 'r_end_shifted', 'read_length', 'strand']].rename(columns={'r_start_shifted':'r2_start_shifted', 'r_end_shifted':'r2_end_shifted'}), ], ignore_index=True)
			read_1_df['read_name'] = read_name_list
			read_2_df['read_name'] = read_name_list
			# read_1_df['read_length'] = read_len
			# read_2_df['read_length'] = read_len
			read_1_df['strand'] = ['+' if x == 1  else '-' for x in read_1_df['strand']]
			read_2_df['strand'] = ['+' if x == 1  else '-' for x in read_2_df['strand']]
			read_1_df_order = read_1_df[['chr','r1_start_shifted', 'r1_end_shifted', 'read_name', 'read_length', 'strand']]
			read_2_df_order = read_2_df[['chr','r2_start_shifted', 'r2_end_shifted', 'read_name', 'read_length', 'strand']]
			if read_1_df_order.shape[0] != read_2_df_order.shape[0]:
				print("[Warning] Grey Area %s read 1 and read 2 not identical!", peak_id)
			if np.sum(np.array(read_1_df_order[['r1_start_shifted']] < 0)) + np.sum(np.array(read_2_df_order[['r2_start_shifted']] < 0)) > 0:
				print("[Warning] Synthetic read pair for Peak %s %s has read 1 or read 2 start position negative: synthetic read pair removed!" % (peak_id, rec_name))
				ind_preserve = np.array(read_1_df_order[['r1_start_shifted']] >= 0) * np.array(read_2_df_order[['r2_start_shifted']] >= 0) # Remove rows with read 1 or read 2's start position negative
				read_1_df_order_removeNegRead = read_1_df_order.loc[ind_preserve]
				read_2_df_order_removeNegRead = read_2_df_order.loc[ind_preserve]
			else:
				read_1_df_order_removeNegRead = read_1_df_order
				read_2_df_order_removeNegRead = read_2_df_order
			read_1_df_order_removeNegRead.to_csv("%s/%s.GreyArea.read1.bed" % (outdirectory, read_bedfile_prename), header=None, index=None, sep='\t', mode='a')
			read_2_df_order_removeNegRead.to_csv("%s/%s.GreyArea.read2.bed" % (outdirectory, read_bedfile_prename), header=None, index=None, sep='\t', mode='a')
			if npair_read_synthetic != read_1_df_order_removeNegRead.shape[0]:
				print("Target read pair %s | Sample synthetic read pair %s" % (npair_read_synthetic, read_1_df_order_removeNegRead.shape[0]))

# Analysis
# Peaks, set 'random_noise_mode = False'
# 35:52
outdirectory = "/home/gayan/Projects/scATAC_Simulator/results/20221121_BoneMarrow_62016_chr1_NONINPUT"
bed_file = outdirectory + "/" + "BoneMarrow_62016_chr1.MACS3_peaks.bed"
count_mat_file = "BoneMarrow_62016_chr1.countmatrix.scDesign2Simulated.txt"
read_bedfile_prename = "BoneMarrow_62016_chr1.syntheticBAM.CBincluded"
INPUT_bamfile = "/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas/BoneMarrow_62016_chr1.bam"
# cellnumberfile = "/home/gayan/Projects/scATAC_Simulator/results/20220408_e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1_1_5000000_NONINPUT/10X_ATAC_chr1_1_5000000.countmatrix.scDesign2Simulated.nPairsRegionmargional.txt"
# BED_filename = "10X_ATAC_chr1_1_5000000.syntheticBAM.CBincluded"
OUTPUT_cells_barcode_file = outdirectory + "/synthetic_cell_barcode.txt"
random_noise_mode = False
grey_area_modeling = False
read_len = 50
jitter_size=5
# Generate bam coordinates for peaks
ATAC_GenerateBAMCoord(bed_file, count_mat_file, read_bedfile_prename, INPUT_bamfile, outdirectory, OUTPUT_cells_barcode_file, jitter_size, read_len, grey_area_modeling)

# Non-peaks, set 'random_noise_mode = True'
# 35:57
outdirectory = "/home/gayan/Projects/scATAC_Simulator/results/20221130_BoneMarrow_62016_chr1_NONINPUT"
bed_file = outdirectory + "/" + "BoneMarrow_62016_chr1.MACS3_thr10_peaks.COMPLE.bed"
count_mat_file = "BoneMarrow_62016_chr1.COMPLE.countmatrix.scDesign2Simulated.txt"
read_bedfile_prename = "BoneMarrow_62016_chr1.syntheticBAM.COMPLE.CBincluded"
INPUT_bamfile = "/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas/BoneMarrow_62016_chr1.bam"
# cellnumberfile = "/home/gayan/Projects/scATAC_Simulator/results/20220408_e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1_1_5000000_NONINPUT/10X_ATAC_chr1_1_5000000.countmatrix.scDesign2Simulated.nPairsRegionmargional.txt"
# BED_filename = "10X_ATAC_chr1_1_5000000.syntheticBAM.CBincluded"
OUTPUT_cells_barcode_file = outdirectory + "/synthetic_cell_barcode.txt"
random_noise_mode = False
grey_area_modeling = True
read_len = 50
jitter_size=5
# Generate bam coordinates for non-peaks and grayareas
ATAC_GenerateBAMCoord(bed_file, count_mat_file, read_bedfile_prename, INPUT_bamfile, outdirectory, OUTPUT_cells_barcode_file, jitter_size, read_len, grey_area_modeling)
