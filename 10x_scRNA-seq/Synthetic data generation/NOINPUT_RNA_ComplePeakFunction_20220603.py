## Construct cell bp matrix
import csv
import numpy as np
import pysam
import pandas as pd
from collections import defaultdict
import sys
import subprocess
import os

def regionLocation(open_peak):
    output = np.asmatrix([[int(x[1]),int(x[2]), int(x[2]) - int(x[1])] for x in open_peak])
    return output

def comple_peak_set(chr_title, peak_directory, peakfile, outdirectory, comple_peakfile):
    # with open("%s/%s" % (peak_directory, peakfile)) as open_peak:
    #     reader = csv.reader(open_peak, delimiter="\t")
    #     open_peak = np.asarray(list(reader))
    open_peak = pd.read_csv("%s/%s" % (peak_directory, peakfile), delimiter="\t",  names=['chr', 'start', 'end'])
    peak_coord = open_peak[['start', 'end']].to_numpy()
    # peak_coord = regionLocation(open_peak)
    # peak_coord.view('i8,i8,i8').sort(order=['f0'], axis=0)
    region_list = list()
    for ind in range(np.shape(open_peak)[0]-2):
        peak_cur = np.asarray(peak_coord[ind,:]).flatten()
        peak_next = np.asarray(peak_coord[ind+1,:]).flatten()
        # region_list.append([chr_title, str(peak_cur[0]), str(peak_cur[1])])
        cur_region_left =  peak_cur[1] + 1
        add_bin_length = peak_next[0] - peak_cur[1] - 1
        if add_bin_length > 0:
            region_list.append([chr_title, str(cur_region_left), str(peak_next[0] - 1)])
        elif add_bin_length < 0:
            print("WARNING: OVERLAPPED Peak" + "peak_next_start" + str(peak_next[1]) + "peak_cur_end" + str(peak_cur[1]))                
    with open("%s/%s" % (outdirectory, comple_peakfile), 'w') as outsfile:
        for item in region_list:
            print("\t".join(item), file=outsfile)


def ExtractBAMCoverage(INPUT_bamfile, samtools_directory, outdirectory):
    """Examine the covered chromosome names for the input bam file.

    Parameters
    ----------
    INPUT_bamfile: `str`
        Directory of input BAM file.
    samtools_directory: `str`
        Directory of software samtools.
    outdirectory: `str`
        Output directory.
    
    Return
    ------
    chromosomes_coverd: `list`
        List of chromosome names that the input bam files covers.
    """
    cmd = "%s/samtools idxstats %s > %s/bam.stats.txt" % (samtools_directory, INPUT_bamfile, outdirectory)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
         print('[ERROR] Fail to rerturn the index BAM information:\n', error.decode())    
    bamstats = pd.read_csv("%s/bam.stats.txt" % outdirectory, header=None, delimiter="\t").to_numpy()
    chromosomes_coverd = bamstats[np.nonzero(bamstats[:,2])[0],0].tolist()
    return chromosomes_coverd


def scRNA_CreateFeatureSets(INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory, genome_annotation, genome_size_file, ref_peakfile, ref_comple_peakfile):
    """Create the foreground and background feature set for the input scRNA-seq bam file.
    Parameters
    ----------
    INPUT_bamfile: `str`
        Input BAM file.
    samtools_directory: `str`
        Path to software `samtools`.
    bedtools_directory: `str`
        Path to software `bedtools`.
    outdirectory: `str`
        Specify the output directory of the features files.
    genome_annotation: `str`
        Genome annotation file for the reference genome that the input BAM aligned on or the synthetic BAM should align on.
    genome_size_file: `str`
        Genome sizes file. The file should be a tab delimited text file with two columns: first column for the chromosome name, second column indicates the size.
    ref_peakfile: `str`
        Specify the name of output foreground feature bed file.
    ref_comple_peakfile: `str`
        Specify the name of output background feature bed file.    
    """
    genome_size_df = pd.read_csv(genome_size_file, header=None, delimiter="\t")
    chromosomes_coverd = ExtractBAMCoverage(INPUT_bamfile, samtools_directory, outdirectory)
    search_string_chr = '|'.join(chromosomes_coverd)
    cmd = "cat %s | grep -Ew '%s' > %s/genome_size_selected.txt" % (genome_size_file, search_string_chr, outdirectory)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
        print('[ERROR] Fail to extract corresponding chromosomes from genome size file:\n', error.decode())
    cmd = """awk -F"\t" '$3=="gene"' %s | cut -f1,4,5 > %s/gene_region.bed""" % (genome_annotation, outdirectory)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
        print('[ERROR] Fail to extract gene regions from genome annotation file:\n', error.decode())
    cmd = "%s/bedtools sort -i %s/gene_region.bed | %s/bedtools merge | grep -Ew '%s' > %s/%s" % (bedtools_directory, outdirectory, bedtools_directory, search_string_chr, outdirectory, ref_peakfile)
    output, error = subprocess.Popen(cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
        print('[ERROR] Fail to create feature set:\n', error.decode())
    os.system("rm %s/gene_region.bed" % outdirectory)
    complement_cmd = "%s/bedtools complement -i %s/%s -g %s/genome_size_selected.txt > %s/%s" % (bedtools_directory, outdirectory, ref_peakfile, outdirectory, outdirectory, ref_comple_peakfile)
    output, error = subprocess.Popen(complement_cmd, shell=True, executable="/bin/bash", stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    if error:
        print('[ERROR] Fail to create complementary feature set:\n', error.decode())
    print('Done!\n') 


def main():
    user_args = sys.argv[1:]
    INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory, genome_annotation, genome_size_file, ref_peakfile, ref_comple_peakfile = user_args
    scRNA_CreateFeatureSets(INPUT_bamfile, samtools_directory, bedtools_directory, outdirectory, genome_annotation, genome_size_file, ref_peakfile, ref_comple_peakfile)

if __name__ == '__main__':
  main()























