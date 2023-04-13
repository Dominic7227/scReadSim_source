'''
Task 2: Generate Per-Cell-Type Ground Truth Peaks
Date: 20230311
'''
import pandas as pd
import csv
import numpy as np
from collections import defaultdict
from tqdm import tqdm
import random
import scReadSim.Utility as Utility
import scReadSim.GenerateSyntheticCount as GenerateSyntheticCount
import scReadSim.scATAC_GenerateBAM as scATAC_GenerateBAM
import os

outdirectory = "/home/gayan/Projects/scATAC_Simulator/results/20230310_SplitCellType/InputGroundTruths" # use absolute path
outdirectory_celltype1 = "/home/gayan/Projects/scATAC_Simulator/results/20230310_SplitCellType/InputGroundTruths/Hematopoieticprogenitors" # use absolute path
outdirectory_celltype2 = "/home/gayan/Projects/scATAC_Simulator/results/20230310_SplitCellType/InputGroundTruths/Erythroblasts" # use absolute path
os.mkdir(outdirectory)
os.mkdir(outdirectory_celltype1)
os.mkdir(outdirectory_celltype2)

# Specify software paths
samtools_directory = "/home/gayan/Tools/samtools/bin"
macs3_directory = "/home/gayan/.local/bin"
bedtools_directory = "/home/gayan/Tools/bedtools/bedtools2/bin"
seqtk_directory = "/home/gayan/Tools/seqtk/seqtk"
bowtie2_directory = "/home/gayan/Tools/bowtie2/bowtie2-2.3.4.1-linux-x86_64"

# Step 1: Prepare Input Data
# Generate input cell-type-wise cell barcode files
cells = pd.read_csv("/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas/cell_metadata_bonemarrow.csv", delimiter="\t")
cells["cell"] = cells["cell"].astype(str)
# cells["cell_label"].value_counts() 
cell_type_selection_list = ["Hematopoietic progenitors", "Erythroblasts", "Monocytes", "Immature B cells"]
## Write out cell barcode file for four celltypes
for cell_type in cell_type_selection_list:
    cells.loc[cells["cell_label"]==cell_type][['cell']].to_csv(outdirectory + "/CB.noCBZ.%s.Real.txt" % cell_type.replace(" ",""), sep="\t", header=False, index=False)
# Specify input files
INPUT_cells_barcode_file_1 = outdirectory + "/CB.noCBZ.Hematopoieticprogenitors.Real.txt"
INPUT_cells_barcode_file_2 = outdirectory + "/CB.noCBZ.Erythroblasts.Real.txt"
filename = "BoneMarrow_62016_chr1"
INPUT_bamfile = "/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas/BoneMarrow_62016_chr1.bam"
# Linux: get mm9 chrm size from ref genome fa file
# pip install pyfaidx
# faidx /home/gayan/Projects/scATAC_Simulator/data/mm9_ref_genome_GECODE/NCBIM37.genome.fa -i chromsizes > /home/gayan/Projects/scATAC_Simulator/results/20230310_SplitCellType/InputGroundTruths/mm9.chrom.sizes
INPUT_genome_size_file = "/home/gayan/Projects/scATAC_Simulator/results/20230310_SplitCellType/InputGroundTruths/mm9.chrom.sizes" # for mm9
cell_type_list = ["Hematopoietic progenitors", "Erythroblasts"]

# Step 2: Prepare Features
# Specify the path to input peak and non-peak bed files
INPUT_peakfile = "/home/gayan/Projects/scATAC_Simulator/results/20221130_BoneMarrow_62016_chr1_INPUT_withCluster_mm9TSS/" + filename + ".MACS3.bed_peaks.bed"
INPUT_nonpeakfile= "/home/gayan/Projects/scATAC_Simulator/results/20221130_BoneMarrow_62016_chr1_INPUT_withCluster_mm9TSS/" + filename + ".MACS3_thr10_peaks.COMPLE.bed"
# Generate cell-type-wise output peak bed files
bulk_peaks = pd.read_csv(INPUT_peakfile, header=None, delimiter="\t")
random.seed(2022)
n_homogeneous_peaks = int(np.floor(bulk_peaks.shape[0]/2))
homogeneous_peaks_id = np.random.choice(bulk_peaks.shape[0], size=n_homogeneous_peaks, replace=False)
heterogeneous_peaks_id = np.setxor1d(np.arange(bulk_peaks.shape[0]), homogeneous_peaks_id)
heterogeneous_peaks_id_1 = np.random.choice(heterogeneous_peaks_id, size=int(np.round(len(heterogeneous_peaks_id)/2)), replace=False)
heterogeneous_peaks_id_2 = np.setxor1d(heterogeneous_peaks_id, heterogeneous_peaks_id_1)
peaks_id_1 = np.concatenate((heterogeneous_peaks_id_1, homogeneous_peaks_id))
peaks_id_1.sort(kind='mergesort')
peaks_id_2 = np.concatenate((heterogeneous_peaks_id_2, homogeneous_peaks_id))
peaks_id_2.sort(kind='mergesort')
# Write out output peak set files
bulk_peaks.iloc[peaks_id_1,:].to_csv(outdirectory + "/GroundTruth.peaks.Hematopoieticprogenitors.txt", sep="\t", header=False, index=False)
bulk_peaks.iloc[peaks_id_2,:].to_csv(outdirectory + "/GroundTruth.peaks.Erythroblasts.txt", sep="\t", header=False, index=False)
# Specify the path to output peak bed files
OUTPUT_peakfile_1 = outdirectory + "/GroundTruth.peaks.Hematopoieticprogenitors.txt"
OUTPUT_peakfile_2 = outdirectory + "/GroundTruth.peaks.Erythroblasts.txt"
# Prepare features with user-specified peaks and non-peaks
Utility.scATAC_CreateFeatureSets(INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory_celltype1, INPUT_genome_size_file, macs3_directory, INPUT_peakfile, INPUT_nonpeakfile, OUTPUT_peakfile=OUTPUT_peakfile_1)
Utility.scATAC_CreateFeatureSets(INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory_celltype2, INPUT_genome_size_file, macs3_directory, INPUT_peakfile, INPUT_nonpeakfile, OUTPUT_peakfile=OUTPUT_peakfile_2)

# Step 3: Count Matrix
# Cell type 1: Hematopoieticprogenitors
# Specify the path to bed files generated by Utility.scATAC_CreateFeatureSets
input_peaks = outdirectory_celltype1 + "/" + "scReadSim.UserInput.peak.bed"
input_nonpeaks = outdirectory_celltype1 + "/" + "scReadSim.UserInput.nonpeak.bed"
output_peaks = outdirectory_celltype1 + "/" + "scReadSim.output.peak.bed"
output_nonpeaks = outdirectory_celltype1 + "/" + "scReadSim.output.nonpeak.bed"
# Specify the names of peak mapping files
assignment_peak_file = filename + ".assigned.peaks.txt"
assignment_nonpeak_file = filename + ".assigned.nonpeaks.txt"
# Generate mappings for peaks and nonpeaks
Utility.FeatureMapping(INPUT_bamfile=INPUT_bamfile, input_peaks=input_peaks, input_nonpeaks=input_nonpeaks, output_peaks=output_peaks, output_nonpeaks=output_nonpeaks, outdirectory=outdirectory_celltype1, assignment_peak_file=assignment_peak_file, assignment_nonpeak_file=assignment_nonpeak_file, n_top=1)
# Specify the output count matrices' base names
count_mat_peak_filename = "%s.output.peak.countmatrix" % filename
count_mat_nonpeak_filename = "%s.output.nonpeak.countmatrix" % filename
# Construct count matrix for peaks
Utility.scATAC_bam2countmat_OutputPeak(cells_barcode_file=INPUT_cells_barcode_file_1, assignment_file=outdirectory_celltype1 + "/" + assignment_peak_file, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory_celltype1, count_mat_filename=count_mat_peak_filename)
# Construct count matrix for non-peaks
Utility.scATAC_bam2countmat_OutputPeak(cells_barcode_file=INPUT_cells_barcode_file_1, assignment_file=outdirectory_celltype1 + "/" + assignment_nonpeak_file, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory_celltype1, count_mat_filename=count_mat_nonpeak_filename)
# Cell type 2: Erythroblasts
# Specify the path to bed files generated by Utility.scATAC_CreateFeatureSets
input_peaks = outdirectory_celltype2 + "/" + "scReadSim.UserInput.peak.bed"
input_nonpeaks = outdirectory_celltype2 + "/" + "scReadSim.UserInput.nonpeak.bed"
output_peaks = outdirectory_celltype2 + "/" + "scReadSim.output.peak.bed"
output_nonpeaks = outdirectory_celltype2 + "/" + "scReadSim.output.nonpeak.bed"
# Specify the names of peak mapping files
assignment_peak_file = filename + ".assigned.peaks.txt"
assignment_nonpeak_file = filename + ".assigned.nonpeaks.txt"
# Generate mappings for peaks and nonpeaks
Utility.FeatureMapping(INPUT_bamfile=INPUT_bamfile, input_peaks=input_peaks, input_nonpeaks=input_nonpeaks, output_peaks=output_peaks, output_nonpeaks=output_nonpeaks, outdirectory=outdirectory_celltype2, assignment_peak_file=assignment_peak_file, assignment_nonpeak_file=assignment_nonpeak_file, n_top=1)
# Specify the output count matrices' base names
count_mat_peak_filename = "%s.output.peak.countmatrix" % filename
count_mat_nonpeak_filename = "%s.output.nonpeak.countmatrix" % filename
# Construct count matrix for peaks
Utility.scATAC_bam2countmat_OutputPeak(cells_barcode_file=INPUT_cells_barcode_file_2, assignment_file=outdirectory_celltype2 + "/" + assignment_peak_file, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory_celltype2, count_mat_filename=count_mat_peak_filename)
# Construct count matrix for non-peaks
Utility.scATAC_bam2countmat_OutputPeak(cells_barcode_file=INPUT_cells_barcode_file_2, assignment_file=outdirectory_celltype2 + "/" + assignment_nonpeak_file, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory_celltype2, count_mat_filename=count_mat_nonpeak_filename)

# Step 4: Synthetic Count Matrix
# Cell type 1: Hematopoieticprogenitors
# Generate synthetic count matrix for peak-by-cell count matrix
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_filename=count_mat_peak_filename, directory=outdirectory_celltype1, outdirectory=outdirectory_celltype1, n_cluster=1)
# Generate synthetic count matrix for nonpeak-by-cell count matrix
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_filename=count_mat_nonpeak_filename, directory=outdirectory_celltype1, outdirectory=outdirectory_celltype1, n_cluster=1)
# Cell type 2: Erythroblasts
# Generate synthetic count matrix for peak-by-cell count matrix
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_filename=count_mat_peak_filename, directory=outdirectory_celltype2, outdirectory=outdirectory_celltype2, n_cluster=1)
# Generate synthetic count matrix for nonpeak-by-cell count matrix
GenerateSyntheticCount.scATAC_GenerateSyntheticCount(count_mat_filename=count_mat_nonpeak_filename, directory=outdirectory_celltype2, outdirectory=outdirectory_celltype2, n_cluster=1)

# Step 5: Generate Read BED Files
# Cell type 1: Hematopoieticprogenitors
# Specify the names of synthetic count matrices (generated by GenerateSyntheticCount.scATAC_GenerateSyntheticCount)
synthetic_countmat_peak_file = count_mat_peak_filename + ".scDesign2Simulated.txt"
synthetic_countmat_nonpeak_file = count_mat_nonpeak_filename + ".scDesign2Simulated.txt"
synthetic_cell_label_file = count_mat_peak_filename + ".scDesign2Simulated.CellTypeLabel.txt"
# Specify the base name of bed files containing synthetic reads
OUTPUT_cells_barcode_file = "synthetic_cell_barcode.txt"
peak_read_bedfile_prename = "%s.syntheticBAM.peak" % filename
nonpeak_read_bedfile_prename = "%s.syntheticBAM.nonpeak" % filename
BED_filename_combined_pre_1 = "Hematopoieticprogenitors.syntheticBAM.combined"
# Create synthetic read bed file for peaks
scATAC_GenerateBAM.scATAC_GenerateBAMCoord_OutputPeak(target_peak_assignment_file=outdirectory_celltype1 + "/" + assignment_peak_file, count_mat_file=outdirectory_celltype1 + "/" + synthetic_countmat_peak_file, synthetic_cell_label_file=outdirectory_celltype1 + "/" + synthetic_cell_label_file, read_bedfile_prename=peak_read_bedfile_prename, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory_celltype1, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=50)
# Create synthetic read bed file for non-peaks
scATAC_GenerateBAM.scATAC_GenerateBAMCoord_OutputPeak(target_peak_assignment_file=outdirectory_celltype1 + "/" + assignment_nonpeak_file, count_mat_file=outdirectory_celltype1 + "/" + synthetic_countmat_nonpeak_file, synthetic_cell_label_file=outdirectory_celltype1 + "/" + synthetic_cell_label_file, read_bedfile_prename=nonpeak_read_bedfile_prename, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory_celltype1, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=50)
# Combine bed files
scATAC_GenerateBAM.scATAC_CombineBED(outdirectory=outdirectory_celltype1, peak_read_bedfile_prename=peak_read_bedfile_prename, nonpeak_read_bedfile_prename=nonpeak_read_bedfile_prename, BED_filename_combined_pre=BED_filename_combined_pre_1, GrayAreaModeling=False)
# Cell type 2: Erythroblasts
# Specify the names of synthetic count matrices (generated by GenerateSyntheticCount.scATAC_GenerateSyntheticCount)
synthetic_countmat_peak_file = count_mat_peak_filename + ".scDesign2Simulated.txt"
synthetic_countmat_nonpeak_file = count_mat_nonpeak_filename + ".scDesign2Simulated.txt"
synthetic_cell_label_file = count_mat_peak_filename + ".scDesign2Simulated.CellTypeLabel.txt"
# Specify the base name of bed files containing synthetic reads
OUTPUT_cells_barcode_file = "synthetic_cell_barcode.txt"
peak_read_bedfile_prename = "%s.syntheticBAM.peak" % filename
nonpeak_read_bedfile_prename = "%s.syntheticBAM.nonpeak" % filename
BED_filename_combined_pre_2 = "Erythroblasts.syntheticBAM.combined"
# Create synthetic read bed file for peaks
scATAC_GenerateBAM.scATAC_GenerateBAMCoord_OutputPeak(target_peak_assignment_file=outdirectory_celltype2 + "/" + assignment_peak_file, count_mat_file=outdirectory_celltype2 + "/" + synthetic_countmat_peak_file, synthetic_cell_label_file=outdirectory_celltype2 + "/" + synthetic_cell_label_file, read_bedfile_prename=peak_read_bedfile_prename, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory_celltype2, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=50)
# Create synthetic read bed file for non-peaks
scATAC_GenerateBAM.scATAC_GenerateBAMCoord_OutputPeak(target_peak_assignment_file=outdirectory_celltype2 + "/" + assignment_nonpeak_file, count_mat_file=outdirectory_celltype2 + "/" + synthetic_countmat_nonpeak_file, synthetic_cell_label_file=outdirectory_celltype2 + "/" + synthetic_cell_label_file, read_bedfile_prename=nonpeak_read_bedfile_prename, INPUT_bamfile=INPUT_bamfile, outdirectory=outdirectory_celltype2, OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file, jitter_size=5, read_len=50)
# Combine bed files
scATAC_GenerateBAM.scATAC_CombineBED(outdirectory=outdirectory_celltype2, peak_read_bedfile_prename=peak_read_bedfile_prename, nonpeak_read_bedfile_prename=nonpeak_read_bedfile_prename, BED_filename_combined_pre=BED_filename_combined_pre_2, GrayAreaModeling=False)

# Step 6: Convert BED files to FASTQ files
# Cell type 1: Hematopoieticprogenitors
referenceGenome_name = "NCBIM37.genome"
referenceGenome_dir = "/home/gayan/Projects/scATAC_Simulator/data/mm9_ref_genome_GECODE"
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)
BED_filename_combined_pre_1 = "Hematopoieticprogenitors.syntheticBAM.combined"
synthetic_fastq_prename_1 = "Hematopoieticprogenitors.syntheticBAM.combined"
# Convert combined bed file into FASTQ files
scATAC_GenerateBAM.scATAC_BED2FASTQ(bedtools_directory=bedtools_directory, seqtk_directory=seqtk_directory, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory_celltype1, BED_filename_combined=BED_filename_combined_pre_1, synthetic_fastq_prename=synthetic_fastq_prename_1)
# Cell type 2: Erythroblasts
# Manually remove read AGGTCAAGATGTAATT:CellNo1090:chr1:197146478-197195432#0000 from read 1 and read 2 bed files
# Reason: AGGTCAAGATGTAATT:CellNo1090:chr1:197146478-197195432#0000 has read 1 larger than chr1 genome sizes
referenceGenome_name = "NCBIM37.genome"
referenceGenome_dir = "/home/gayan/Projects/scATAC_Simulator/data/mm9_ref_genome_GECODE"
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)
BED_filename_combined_pre_2 = "Erythroblasts.syntheticBAM.combined"
synthetic_fastq_prename_2 = "Erythroblasts.syntheticBAM.combined"
# Convert combined bed file into FASTQ files
scATAC_GenerateBAM.scATAC_BED2FASTQ(bedtools_directory=bedtools_directory, seqtk_directory=seqtk_directory, referenceGenome_file=referenceGenome_file, outdirectory=outdirectory_celltype2, BED_filename_combined=BED_filename_combined_pre_2, synthetic_fastq_prename=synthetic_fastq_prename_2)

# Step 7: Read alignment
# Cell type 1: Hematopoieticprogenitors
output_BAM_pre_1 =  "Hematopoieticprogenitors.syntheticBAM.combined"
# Convert FASTQ files to BAM file
scATAC_GenerateBAM.AlignSyntheticBam_Pair(bowtie2_directory=bowtie2_directory, samtools_directory=samtools_directory, outdirectory=outdirectory_celltype1, referenceGenome_name=referenceGenome_name, referenceGenome_dir=referenceGenome_dir, synthetic_fastq_prename=synthetic_fastq_prename_1, output_BAM_pre=output_BAM_pre_1)
# Cell type 2: Erythroblasts
output_BAM_pre_2 =  "Erythroblasts.syntheticBAM.combined"
# Convert FASTQ files to BAM file
scATAC_GenerateBAM.AlignSyntheticBam_Pair(bowtie2_directory=bowtie2_directory, samtools_directory=samtools_directory, outdirectory=outdirectory_celltype2, referenceGenome_name=referenceGenome_name, referenceGenome_dir=referenceGenome_dir, synthetic_fastq_prename=synthetic_fastq_prename_2, output_BAM_pre=output_BAM_pre_2)
