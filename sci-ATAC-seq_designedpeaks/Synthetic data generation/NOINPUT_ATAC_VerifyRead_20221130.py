import csv
import numpy as np
import pysam
import pandas as pd
from collections import defaultdict
from time import process_time
import sys
import subprocess
from tqdm import tqdm

def bam2MarginalCount(outdirectory, bed_file, sam_filename, count_mat_filename):
    with open("%s/%s" % (outdirectory, count_mat_filename), 'w') as outsfile:
        samfile = pysam.AlignmentFile(sam_filename, "rb")
        with open("%s/%s" % (outdirectory, bed_file)) as open_peak:
            reader = csv.reader(open_peak, delimiter="\t")
            open_peak = np.asarray(list(reader))
        k = 0
        peaksdic = defaultdict(lambda: [None])
        for rec in open_peak:
            rec_name = '_'.join(rec)
            peaksdic[rec_name] = k
            k += 1
        peaks_n = len(open_peak)
        print("Converting count matrix...\n")
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
    outdirectory, bed_file, sam_filename, output_count_vec_filename = user_args
    bam2MarginalCount(outdirectory, bed_file, sam_filename, output_count_vec_filename)

if __name__ == '__main__':
  	main()
