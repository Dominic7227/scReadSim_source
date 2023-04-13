import csv
import numpy as np
import pysam
import pandas as pd
from collections import defaultdict
from time import process_time
import sys


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
        print("Converting count matrix...\n")
        for rec in open_peak:
            rec_name = '_'.join(rec)
            peaksdic[rec_name] = k
            k += 1
        cells_n = len(cells_barcode)
        peaks_n = len(open_peak)
        for rec in open_peak:
            rec_name = '_'.join(rec)
            currcounts = [0]*cells_n
            for read1, read2 in read_pair_generator(samfile, rec[3], int(rec[4]), int(rec[5])):
                read1_str = str(read1).split("\t")
                read2_str = str(read2).split("\t")
                cell = read1_str[0].split(':')[0]
                if cell in cells_barcode:
                    try:
                        currcounts[cellsdic[cell]] += 2
                    except KeyError:
                        pass
            for read in read_pair_generator_supp(samfile, rec[3], int(rec[4]), int(rec[5])):
                read_str = str(read).split("\t")
                cell = read_str[0].split(':')[0]
                if cell in cells_barcode:
                    try:
                        currcounts[cellsdic[cell]] += 1
                    except KeyError:
                        pass
            # print(sum(currcounts))
            # if sum(currcounts) > 0:
            print(rec_name + "\t" + "\t".join([str(x) for x in currcounts]),file = outsfile)

def main():
    user_args = sys.argv[1:]
    cells_barcode_file, bed_directory, bed_file, sam_filename, outdirectory, count_mat_filename = user_args
    bam2countmat(cells_barcode_file, bed_directory, bed_file, sam_filename, outdirectory, count_mat_filename)

if __name__ == '__main__':
    main()
