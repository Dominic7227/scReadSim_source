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

### Synethetic MACS3
# filename = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1"
# outdirectory = "/home/gayan/Projects/scATAC_Simulator/results/20220314_%s_NONINPUT" % filename

# sam_filename = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.syntheticBAM.CBincluded.synthetic.chr1.sorted.bam"
# ## 
# bed_file = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.syntheticBAM.synthetic.chr1.sorted.bam.MACS3_peaks.narrowPeak.bed"
# count_mat_filename = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.syntheticBAM.MACS3Peaks.MarginalFeatureCount.txt"
# bed_comple_file = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.syntheticBAM.synthetic.chr1.sorted.bam.MACS3_peaks.narrowPeak.COMPLE.bed"
# count_comple_mat_filename = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.syntheticBAM.MACS3Peaks.COMPLE.MarginalFeatureCount.txt"
# bam2MarginalCount(outdirectory, bed_file, sam_filename, outdirectory, count_mat_filename)
# bam2MarginalCount(outdirectory, bed_comple_file, sam_filename, outdirectory, count_comple_mat_filename)

# ## Real MACS3
# real_sam_filename = "/home/gayan/Projects/scATAC_Simulator/data/10X_MOUSE_BRAIN_ATACandRNA/split.e18_mouse_brain_fresh_5k_atac_possorted_bam/e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.bam"
# bed_file = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.MACS3_peaks.bed"
# count_mat_filename = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.countmatrix.real.MarginalFeatureCount.txt"
# bed_comple_file = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.MACS3_peaks.COMPLE.bed"
# count_comple_mat_filename = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.COMPLE.countmatrix.real.MarginalFeatureCount.txt"
# bam2MarginalCount(outdirectory, bed_file, real_sam_filename, outdirectory, count_mat_filename)
# bam2MarginalCount(outdirectory, bed_comple_file, real_sam_filename, outdirectory, count_comple_mat_filename)


# ### Synethetic MACS3 with less stringent rules
# filename = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1"
# outdirectory = "/home/gayan/Projects/scATAC_Simulator/results/20220314_%s_NONINPUT" % filename
# sam_filename = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.syntheticBAM.CBincluded.synthetic.chr1.sorted.bam"
# ## 
# bed_file = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.syntheticBAM.synthetic.chr1.sorted.bam.MACS3.lessStringent_peaks.narrowPeak.bed"
# count_mat_filename = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.syntheticBAM.MACS3Peaks.lessStringent.MarginalFeatureCount.txt"
# bed_comple_file = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.syntheticBAM.synthetic.chr1.sorted.bam.MACS3.lessStringent_peaks.narrowPeak.COMPLE.bed"
# count_comple_mat_filename = "e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1.syntheticBAM.MACS3Peaks.lessStringent.COMPLE.MarginalFeatureCount.txt"
# bam2MarginalCount(outdirectory, bed_file, sam_filename, outdirectory, count_mat_filename)
# bam2MarginalCount(outdirectory, bed_comple_file, sam_filename, outdirectory, count_comple_mat_filename)

# Used in main_20231130_GreyArea.sh
def main():
    user_args = sys.argv[1:]
    outdirectory, bed_file, sam_filename, output_count_vec_filename = user_args
    bam2MarginalCount(outdirectory, bed_file, sam_filename, output_count_vec_filename)

if __name__ == '__main__':
  	main()
