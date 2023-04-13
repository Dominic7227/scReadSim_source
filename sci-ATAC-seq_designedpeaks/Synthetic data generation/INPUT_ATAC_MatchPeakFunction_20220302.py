import sys
import csv
import numpy as np
import pysam
import pandas as pd
from collections import defaultdict
from time import process_time
import sys
import subprocess
from tqdm import tqdm

def find_nearest_peak(array, value, k, ref_read_density):
    array = np.asarray(array)
    # idx = (np.abs(array - value)).argmin()   
    idx = (np.abs(array - value)).argpartition(k)
    idx_ReadDensity = ref_read_density[idx[:k]].argmax()
    final_idx = idx[idx_ReadDensity]
    return final_idx

def find_nearest_comple(array, value, k, ref_read_density):
    array = np.asarray(array)
    # idx = (np.abs(array - value)).argmin()   
    idx = (np.abs(array - value)).argpartition(k)
    idx_ReadDensity = ref_read_density[idx[:k]].argmin()
    final_idx = idx[idx_ReadDensity]
    return final_idx

def fragment_length(open_peak):
	output = np.asarray([int(x[2]) - int(x[1]) for x in open_peak])
	return output


def match_peak(ref_marginal_file ,true_peak_directory, true_peakfile, ref_peak_directory, ref_peakfile, outdirectory, assignment_file):
    with open("%s/%s" % (true_peak_directory, true_peakfile)) as true_peak_file:
        reader = csv.reader(true_peak_file, delimiter="\t")
        true_peak_set = np.asarray(list(reader))
    with open("%s/%s" % (ref_peak_directory, ref_peakfile)) as ref_peak_file:
        reader = csv.reader(ref_peak_file, delimiter="\t")
        ref_peak_set = np.asarray(list(reader))
    with open("%s/%s" % (outdirectory, ref_marginal_file)) as ref_marginal_file_dir:
        reader = csv.reader(ref_marginal_file_dir, delimiter="\t")
        ref_marginal_count = np.asarray(list(reader))
    ref_peak_fraglen = fragment_length(ref_peak_set)
    ref_read_density = ref_marginal_count[:,1].astype(int) / ref_peak_fraglen
    # ref_peak_fraglen.view('i8,i8,i8').sort(order=['f0'], axis=0)
    with open("%s/%s" % (outdirectory, assignment_file), 'w') as outsfile:
        for true_peak_id in tqdm(range(len(true_peak_set))):
            true_peak = true_peak_set[true_peak_id]
            true_length = int(true_peak[2]) - int(true_peak[1])
            idx = find_nearest_peak(ref_peak_fraglen, true_length, 50, ref_read_density)
            print("\t".join(true_peak[0:3]) + '\t' + "\t".join(ref_peak_set[idx][0:3]), file=outsfile)

def match_complepeak(ref_marginal_file ,true_peak_directory, true_peakfile, ref_peak_directory, ref_peakfile, outdirectory, assignment_file):
    with open("%s/%s" % (true_peak_directory, true_peakfile)) as true_peak_file:
        reader = csv.reader(true_peak_file, delimiter="\t")
        true_peak_set = np.asarray(list(reader))
    with open("%s/%s" % (ref_peak_directory, ref_peakfile)) as ref_peak_file:
        reader = csv.reader(ref_peak_file, delimiter="\t")
        ref_peak_set = np.asarray(list(reader))
    with open("%s/%s" % (outdirectory, ref_marginal_file)) as ref_marginal_file_dir:
        reader = csv.reader(ref_marginal_file_dir, delimiter="\t")
        ref_marginal_count = np.asarray(list(reader))
    ref_peak_fraglen = fragment_length(ref_peak_set)
    ref_read_density = ref_marginal_count[:,1].astype(int) / ref_peak_fraglen
    # ref_peak_fraglen.view('i8,i8,i8').sort(order=['f0'], axis=0)
    with open("%s/%s" % (outdirectory, assignment_file), 'w') as outsfile:
        for true_peak_id in tqdm(range(len(true_peak_set))):
            true_peak = true_peak_set[true_peak_id]
            true_length = int(true_peak[2]) - int(true_peak[1])
            idx = find_nearest_comple(ref_peak_fraglen, true_length, 50, ref_read_density)
            print("\t".join(true_peak[0:3]) + '\t' + "\t".join(ref_peak_set[idx][0:3]), file=outsfile)

def bam2MarginalCount(bed_directory, bed_file, sam_filename, outdirectory, ref_marginal_file):
    with open("%s/%s" % (outdirectory, ref_marginal_file), 'w') as outsfile:
        samfile = pysam.AlignmentFile(sam_filename, "rb")
        with open("%s/%s" % (bed_directory, bed_file)) as open_peak:
            reader = csv.reader(open_peak, delimiter="\t")
            open_peak = np.asarray(list(reader))
        k = 0
        peaksdic = defaultdict(lambda: [None])
        for rec in open_peak:
            rec_name = '_'.join(rec)
            peaksdic[rec_name] = k
            k += 1
        peaks_n = len(open_peak)
        print("Converting marginal count vector...\n")
        # for rec in open_peak:
        for rec_id in tqdm(range(len(open_peak))):
            rec = open_peak[rec_id]
            rec_name = '_'.join(rec)
            currcounts =  0
            reads = samfile.fetch(rec[0], int(rec[1]), int(rec[2]))
            for read in reads:
                currcounts += 1
            # if sum(currcounts) > 0:
            print(rec_name + "\t" + str(currcounts),file = outsfile)


def main():
    user_args = sys.argv[1:]
    sam_filename, ref_marginal_file, ref_marginal_comple_file, true_peak_directory, true_peakfile, true_comple_peakfile, ref_peak_directory, ref_peakfile, ref_comple_peakfile, outdirectory, assignment_file, assignment_comple_file = user_args
    bam2MarginalCount(ref_peak_directory, ref_peakfile, sam_filename, outdirectory, ref_marginal_file)
    bam2MarginalCount(ref_peak_directory, ref_comple_peakfile, sam_filename, outdirectory, ref_marginal_comple_file)
    match_peak(ref_marginal_file, true_peak_directory, true_peakfile, ref_peak_directory, ref_peakfile, outdirectory, assignment_file)
    match_complepeak(ref_marginal_comple_file, true_peak_directory, true_comple_peakfile, ref_peak_directory, ref_comple_peakfile, outdirectory, assignment_comple_file)

if __name__ == '__main__':
  	main()


