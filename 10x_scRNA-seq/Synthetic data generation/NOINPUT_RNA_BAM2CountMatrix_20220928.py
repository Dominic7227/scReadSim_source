import csv
import numpy as np
import pysam
import pandas as pd
from collections import defaultdict
from time import process_time
import sys
from tqdm import tqdm

## Construct expression matrix for single paired read
def read_pair_generator(bam, contig_str=None,start_int = None,stop_int=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(contig=contig_str, start=start_int, stop=stop_int):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary or read.is_duplicate:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

def read_pair_generator_supp(bam, contig_str=None,start_int = None,stop_int=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(contig=contig_str, start=start_int, stop=stop_int):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary or read.is_duplicate:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            del read_dict[qname]
    single_reads = []
    for qname in read_dict:
        if read_dict[qname][0] == None:
            single_reads.append(read_dict[qname][1])
        else:
            single_reads.append(read_dict[qname][0])
    return single_reads

def getBarcodesFromBam(sam_filename):
    """Identify unique barcodes from the bam file
    
    Args:
        sam_filename: a bam file
    Returns:
        A dictionary contains all barcodes, otherwise None
    """
    barcode_dict = collections.OrderedDict();
    samfile = pysam.AlignmentFile(sam_filename, "rb");
    i = 1
    for _read in samfile:
        barcode = _read.qname.split(":")[0].upper();
        if barcode not in barcode_dict:
            barcode_dict[barcode] = i
            i = i + 1
    samfile.close();
    return barcode_dict

def bam2countmat(cells_barcode_file, bed_directory, bed_file, sam_filename, outdirectory, count_mat_filename):
    cells = pd.read_csv(cells_barcode_file, sep="\t")
    cells = cells.values.tolist()
    cells_barcode = [item[0] for item in cells]
    with open("%s/%s" % (outdirectory, count_mat_filename), 'w') as outsfile:
        samfile = pysam.AlignmentFile(sam_filename, "rb")
        with open("%s/%s" % (bed_directory, bed_file)) as open_peak:
            reader = csv.reader(open_peak, delimiter="\t")
            open_peak = np.asarray(list(reader))
        k = 0
        cellsdic = defaultdict(lambda: [None])
        for cell in cells_barcode:
            cellsdic[cell] = k
            k += 1
        k = 0
        peaksdic = defaultdict(lambda: [None])
        for rec in open_peak:
            rec_name = '_'.join(rec)
            peaksdic[rec_name] = k
            k += 1
        cells_n = len(cells_barcode)
        peaks_n = len(open_peak)
        print("Converting count matrix...\n")
        for rec in open_peak:
            rec_name = '_'.join(rec)
            reads = samfile.fetch(rec[0], int(rec[1]), int(rec[2]))
            currcounts = [0]*cells_n
            for read in reads:
                if read.has_tag('CB:Z'):
                    cell = read.get_tag('CB:Z')
                    # print(cell)
                    if cell in cells_barcode:
                        try:
                            currcounts[cellsdic[cell]] += 1
                        except KeyError:
                            pass
            # print(sum(currcounts))
            # if sum(currcounts) > 0:
            print(rec_name + "\t" + "\t".join([str(x) for x in currcounts]),file = outsfile)


def bam2countmat(cells_barcode_file, bed_file, INPUT_bamfile, outdirectory, count_mat_filename):
    """Construct count matrix.
    Parameters
    ----------
    cells_barcode_file: `str`
        Cell barcode file corresponding to the input BAM file.
    bed_file: `str`
        Features bed file to generate the count matrix.
    INPUT_bamfile: `str`
        Input BAM file for anlaysis.
    outdirectory: `str`
        Specify the output directory of the count matrix file.
    count_mat_filename: `str`
        Specify the base name of output count matrix.
    """
    cells = pd.read_csv(cells_barcode_file, sep="\t")
    cells = cells.values.tolist()
    cells_barcode = [item[0] for item in cells]
    with open("%s/%s.txt" % (outdirectory, count_mat_filename), 'w') as outsfile:
        samfile = pysam.AlignmentFile(INPUT_bamfile, "rb")
        with open(bed_file) as open_peak:
            reader = csv.reader(open_peak, delimiter="\t")
            open_peak = np.asarray(list(reader))
        k = 0
        cellsdic = defaultdict(lambda: [None])
        for cell in cells_barcode:
            cellsdic[cell] = k
            k += 1
        k = 0
        peaksdic = defaultdict(lambda: [None])
        for rec in open_peak:
            rec_name = '_'.join(rec)
            peaksdic[rec_name] = k
            k += 1
        cells_n = len(cells_barcode)
        peaks_n = len(open_peak)
        # marginal_count_vec = [0] * len(open_peak)
        print("Converting count matrix...\n")
        # for rec in open_peak:
        for rec_id in tqdm(range(len(open_peak))):
            rec = open_peak[rec_id]
            rec_name = '_'.join(rec)
            currcounts = [0]*cells_n
            reads = samfile.fetch(rec[0], int(rec[1]), int(rec[2]))
            for read in reads:
                cell = read.qname.split(":")[0].upper()
                if cell in cells_barcode:
                    try:
                        currcounts[cellsdic[cell]] += 1
                    except KeyError:
                        pass
            # marginal_count_vec[rec_id] = sum(currcounts)
            # if sum(currcounts) > 0:
            print(rec_name + "\t" + "\t".join([str(x) for x in currcounts]),file = outsfile)

def bam2countmat_UMI(cells_barcode_file, bed_file, INPUT_bamfile, outdirectory, UMI_count_mat_filename):
    """Construct count matrix.
    Parameters
    ----------
    cells_barcode_file: `str`
        Cell barcode file corresponding to the input BAM file.
    bed_file: `str`
        Features bed file to generate the count matrix.
    INPUT_bamfile: `str`
        Input BAM file for anlaysis.
    outdirectory: `str`
        Specify the output directory of the count matrix file.
    """
    cells = pd.read_csv(cells_barcode_file, sep="\t")
    cells = cells.values.tolist()
    cells_barcode = [item[0] for item in cells]
    samfile = pysam.AlignmentFile(INPUT_bamfile, "rb")
    with open(bed_file) as open_peak:
        reader = csv.reader(open_peak, delimiter="\t")
        open_peak = np.asarray(list(reader))
    k = 0
    cellsdic = defaultdict(lambda: [None])
    for cell in cells_barcode:
        cellsdic[cell] = k
        k += 1
    k = 0
    peaksdic = defaultdict(lambda: [None])
    for rec in open_peak:
        rec_name = '_'.join(rec)
        peaksdic[rec_name] = k
        k += 1
    cells_n = len(cells_barcode)
    peaks_n = len(open_peak)
    UMI_count_mat = np.zeros((peaks_n,cells_n), dtype="int")
    # marginal_count_vec = [0] * len(open_peak)
    print("Converting UMI count matrix...\n")
    # # for rec in open_peak:
    # with open("%s/%s.txt" % (outdirectory, count_mat_filename), 'w') as outsfile:
    #     for rec_id in tqdm(range(len(open_peak))):
    #         rec = open_peak[rec_id]
    #         rec_name = '_'.join(rec)
    #         read_currcounts = [0]*cells_n
    #         UMI_currlist  = [["empty UMI"] for _ in range(cells_n)]  # Create netsed list to store UMI for each cell within one peak. Remeber to minus 1 for "empty UMI" when count unique UMIs.
    #         reads = samfile.fetch(rec[0], int(rec[1]), int(rec[2]))
    #         for read in reads:
    #             cell = read.qname.split(":")[0].upper()
    #             if cell in cells_barcode:
    #                 try:
    #                     read_currcounts[cellsdic[cell]] += 1
    #                     if read.has_tag('UB:Z'):
    #                         UMI = read.get_tag('UB:Z')
    #                         UMI_currlist[cellsdic[cell]].append(UMI)
    #                 except KeyError:
    #                     pass
    #         # marginal_count_vec[rec_id] = sum(currcounts)
    #         # if sum(currcounts) > 0:
    #         UMI_count_mat[rec_id,:] = [len(set(UMIs_percell))-1 for UMIs_percell in UMI_currlist]
    #         print(rec_name + "\t" + "\t".join([str(x) for x in read_currcounts]),file = outsfile)
    with open("%s/%s.txt" % (outdirectory, UMI_count_mat_filename), 'w') as outsfile:
        for rec_id in tqdm(range(len(open_peak))):
            rec = open_peak[rec_id]
            rec_name = '_'.join(rec)
            UMI_currlist  = [["empty UMI"] for _ in range(cells_n)]  # Create netsed list to store UMI for each cell within one peak. Remeber to minus 1 for "empty UMI" when count unique UMIs.
            reads = samfile.fetch(rec[0], int(rec[1]), int(rec[2]))
            for read in reads:
                cell = read.qname.split(":")[0].upper()
                if cell in cells_barcode:
                    try:
                        if read.has_tag('UB:Z'):
                            UMI = read.get_tag('UB:Z')
                            UMI_currlist[cellsdic[cell]].append(UMI)
                    except KeyError:
                        pass
            UMI_count_mat[rec_id,:] = [len(set(UMIs_percell))-1 for UMIs_percell in UMI_currlist]
            print(rec_name + "\t" + "\t".join([str(x) for x in UMI_count_mat[rec_id,:]]),file = outsfile)

def bam2countmat_UMIandRead(cells_barcode_file, bed_file, INPUT_bamfile, outdirectory, count_mat_filename, UMI_count_mat_filename):
    """Construct count matrix.
    Parameters
    ----------
    cells_barcode_file: `str`
        Cell barcode file corresponding to the input BAM file.
    bed_file: `str`
        Features bed file to generate the count matrix.
    INPUT_bamfile: `str`
        Input BAM file for anlaysis.
    outdirectory: `str`
        Specify the output directory of the count matrix file.
    count_mat_filename: `str`
        Specify the base name of output count matrix.
    """
    cells = pd.read_csv(cells_barcode_file, sep="\t")
    cells = cells.values.tolist()
    cells_barcode = [item[0] for item in cells]
    samfile = pysam.AlignmentFile(INPUT_bamfile, "rb")
    with open(bed_file) as open_peak:
        reader = csv.reader(open_peak, delimiter="\t")
        open_peak = np.asarray(list(reader))
    k = 0
    cellsdic = defaultdict(lambda: [None])
    for cell in cells_barcode:
        cellsdic[cell] = k
        k += 1
    k = 0
    peaksdic = defaultdict(lambda: [None])
    for rec in open_peak:
        rec_name = '_'.join(rec)
        peaksdic[rec_name] = k
        k += 1
    cells_n = len(cells_barcode)
    peaks_n = len(open_peak)
    UMI_count_mat = np.zeros((peaks_n,cells_n), dtype="int")
    # marginal_count_vec = [0] * len(open_peak)
    print("Converting count matrix...\n")
    # for rec in open_peak:
    with open("%s/%s.txt" % (outdirectory, count_mat_filename), 'w') as outsfile:
        for rec_id in tqdm(range(len(open_peak))):
            rec = open_peak[rec_id]
            rec_name = '_'.join(rec)
            read_currcounts = [0]*cells_n
            UMI_currlist  = [["empty UMI"] for _ in range(cells_n)]  # Create netsed list to store UMI for each cell within one peak. Remeber to minus 1 for "empty UMI" when count unique UMIs.
            reads = samfile.fetch(rec[0], int(rec[1]), int(rec[2]))
            for read in reads:
                cell = read.qname.split(":")[0].upper()
                if cell in cells_barcode:
                    try:
                        read_currcounts[cellsdic[cell]] += 1
                        if read.has_tag('UB:Z'):
                            UMI = read.get_tag('UB:Z')
                            UMI_currlist[cellsdic[cell]].append(UMI)
                    except KeyError:
                        pass
            # marginal_count_vec[rec_id] = sum(currcounts)
            # if sum(currcounts) > 0:
            UMI_count_mat[rec_id,:] = [len(set(UMIs_percell))-1 for UMIs_percell in UMI_currlist]
            print(rec_name + "\t" + "\t".join([str(x) for x in read_currcounts]),file = outsfile)
    with open("%s/%s.txt" % (outdirectory, UMI_count_mat_filename), 'w') as outsfile:
        for rec_id in tqdm(range(len(open_peak))):
            rec = open_peak[rec_id]
            rec_name = '_'.join(rec)
            print(rec_name + "\t" + "\t".join([str(x) for x in UMI_count_mat[rec_id,:]]),file = outsfile)

# INPUT_bamfile = "/home/gayan/Projects/scATAC_Simulator/results/20220603_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/filtered_CB.bam"
# outdirectory = "/home/gayan/Projects/scATAC_Simulator/results/20220603_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT"
# bed_file = "/home/gayan/Projects/scATAC_Simulator/results/20220603_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_peaks.bed"
# cells_barcode_file = "/home/gayan/Projects/scATAC_Simulator/data/10X_MOUSE_BRAIN_ATACandRNA/filtered_feature_bc_matrix/barcodes.tsv"
# count_mat_filename = "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.countmatrix"
# UMI_count_mat_filename = "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.UMI.countmatrix"

def main():
    user_args = sys.argv[1:]
    cells_barcode_file, bed_file, sam_filename, outdirectory, UMI_count_mat_filename = user_args
    bam2countmat_UMI(cells_barcode_file, bed_file, sam_filename, outdirectory, UMI_count_mat_filename)

if __name__ == '__main__':
    main()
