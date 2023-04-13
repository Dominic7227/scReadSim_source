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


# 20221122 TODO: Update the sampling proceuder with random_noise_mode and target_peak_mode
def ATAC_GenerateBAMCoord_InputTargetPeak(target_peak_assignment_file, count_mat_file, read_bedfile_prename, INPUT_bamfile, outdirectory, OUTPUT_cells_barcode_file, jitter_size=5, read_len=50, random_noise_mode = False):
    """
    true_peak_file: peak list taken as ground truth peaks
    target_peak_assignment_file: assignent file between true peaks and target peaks
    """
    count_mat_df = pd.read_csv("%s/%s" % (outdirectory, count_mat_file), header=0, delimiter="\t")
    count_mat = count_mat_df.to_numpy()
    count_mat_cluster = count_mat_df.columns.to_numpy()
    n_cell = np.shape(count_mat)[1]
    samfile = pysam.AlignmentFile(INPUT_bamfile, "rb")
    with open(target_peak_assignment_file) as open_peak:
        reader = csv.reader(open_peak, delimiter="\t")
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
    with open("%s/%s.read1.test.bed" % (outdirectory, read_bedfile_prename), 'w') as f_1:
        pass

    with open("%s/%s.read2.test.bed" % (outdirectory, read_bedfile_prename), 'w') as f_2:
        pass

    # w/ Target Peak
    for relative_peak_ind in tqdm(range(len(peak_nonzero_id))):
        peak_ind = peak_nonzero_id[relative_peak_ind]
        rec = open_peak[peak_ind]
        rec_name = '_'.join(rec)
        true_peak_concat = rec[3] + ":" + str(rec[4]) + "-" + str(rec[5])
        target_peak_concat = rec[0] + ":" + str(rec[1]) + "-" + str(rec[2])
        if int(rec[2]) - int(rec[1]) == 0:
            print("Peak %s has identical start and end position. Skip.")
            continue
        shift_number = int(rec[1]) - int(rec[4])
        reads = samfile.fetch(rec[3], int(rec[4]), int(rec[5])) # Extract reads from true peaks
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
                read_synthetic_start = read_sampled_unsplit[:,0].astype(int) + shift_number
            # Generate Read Name
            nonempty_cell_ind = np.where(count_frag_vec != 0)[0]
            target_peak_concat = rec[0] + ":" + str(rec[1]) + "-" + str(rec[2])
            read_name_list = [random_cellbarcode_list[nonempty_cell_ind[ind]] + ":"+ "CellNo" + str(nonempty_cell_ind[ind] + 1) + ":" + str(target_peak_concat) + "#" + str(count).zfill(4) for ind in range(len(nonempty_cell_ind)) for count in range(count_frag_vec[nonempty_cell_ind[ind]])]
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
            # read_1_df_order.to_csv("%s/%s.read1.test.bed" % (outdirectory, read_bedfile_prename), header=None, index=None, sep='\t', mode='a')
            # read_2_df_order.to_csv("%s/%s.read2.test.bed" % (outdirectory, read_bedfile_prename), header=None, index=None, sep='\t', mode='a')
            if read_2_df_order.shape[0] != read_2_df_order.shape[0]:
                print("Peak %s read 1 and read 2 not identical!", relative_peak_ind)
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



# Peak 
outdirectory = "/home/gayan/Projects/scATAC_Simulator/results/20221130_BoneMarrow_62016_chr1_INPUT_withCluster_mm9TSS"
# true_peak_file = outdirectory + "/" + "BoneMarrow_62016_chr1.MACS3.bed_peaks.bed"
target_peak_assignment_file = outdirectory + "/" + "BoneMarrow_62016_chr1.assigned.peaks.txt"
count_mat_file = "BoneMarrow_62016_chr1.assigned.countmatrix.scDesign2Simulated.txt"
read_bedfile_prename = "BoneMarrow_62016_chr1.syntheticBAM"
INPUT_bamfile = "/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas/BoneMarrow_62016_chr1.bam"
OUTPUT_cells_barcode_file = outdirectory + "/synthetic_cell_barcode.txt"
random_noise_mode = False
read_len = 50
jitter_size=5
# Generate read coordinate for peaks
ATAC_GenerateBAMCoord_InputTargetPeak(target_peak_assignment_file, count_mat_file, read_bedfile_prename, INPUT_bamfile, outdirectory, OUTPUT_cells_barcode_file, jitter_size, read_len, random_noise_mode)


# Non-Peak 
outdirectory = "/home/gayan/Projects/scATAC_Simulator/results/20221130_BoneMarrow_62016_chr1_INPUT_withCluster_mm9TSS"
# true_peak_file = outdirectory + "/" + "BoneMarrow_62016_chr1.MACS3.bed_peaks.bed"
target_peak_assignment_file = outdirectory + "/" + "BoneMarrow_62016_chr1.COMPLE.assigned.peaks.txt"
count_mat_file = "BoneMarrow_62016_chr1.assigned.COMPLE.countmatrix.scDesign2Simulated.txt"
read_bedfile_prename = "BoneMarrow_62016_chr1.syntheticBAM.COMPLE"
INPUT_bamfile = "/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas/BoneMarrow_62016_chr1.bam"
OUTPUT_cells_barcode_file = outdirectory + "/synthetic_cell_barcode.txt"
random_noise_mode = False
read_len = 50
jitter_size = 5
# Generate read coordinate for non-peaks
ATAC_GenerateBAMCoord_InputTargetPeak(target_peak_assignment_file, count_mat_file, read_bedfile_prename, INPUT_bamfile, outdirectory, OUTPUT_cells_barcode_file, jitter_size, read_len, random_noise_mode)
