## Construct cell bp matrix
import csv
import numpy as np
import pysam
import pandas as pd
from collections import defaultdict
import sys

def regionLocation(open_peak):
    output = np.asmatrix([[int(x[1]),int(x[2]), int(x[2]) - int(x[1])] for x in open_peak])
    return output

def comple_peak_set(chr_title, peak_directory, peakfile, outdirectory, comple_peakfile):
    with open("%s/%s" % (peak_directory, peakfile)) as open_peak:
        reader = csv.reader(open_peak, delimiter="\t")
        open_peak = np.asarray(list(reader))
    peak_coord = open_peak[:,[1,2]]
    # open_peak = pd.read_csv("%s/%s" % (peak_directory, peakfile), delimiter="\t",  names=['chr', 'start', 'end'])
    # peak_coord = open_peak[['start', 'end']].to_numpy()
    # peak_coord = regionLocation(open_peak)
    # peak_coord.view('i8,i8,i8').sort(order=['f0'], axis=0)
    region_list = list()
    for ind in range(np.shape(open_peak)[0]-2):
        peak_cur = np.asarray(peak_coord[ind,:]).flatten()
        peak_next = np.asarray(peak_coord[ind+1,:]).flatten()
        # region_list.append([chr_title, str(peak_cur[0]), str(peak_cur[1])])
        cur_region_left = int(peak_cur[1]) + 1
        add_bin_length = int(peak_next[0]) - int(peak_cur[1]) - 1
        if add_bin_length > 0:
            region_list.append([chr_title, str(cur_region_left), str(int(peak_next[0]) - 1)])
        elif add_bin_length < 0:
            print("WARNING: OVERLAPPED Peak" + "peak_next_start" + str(peak_next[1]) + "peak_cur_end" + str(peak_cur[1]))                
    with open("%s/%s" % (outdirectory, comple_peakfile), 'w') as outsfile:
        for item in region_list:
            print("\t".join(item), file=outsfile)

def main():
    user_args = sys.argv[1:]
    chr_title, peak_directory, peakfile, outdirectory, comple_peakfile = user_args
    comple_peak_set(chr_title, peak_directory, peakfile, outdirectory, comple_peakfile)

if __name__ == '__main__':
  main()















