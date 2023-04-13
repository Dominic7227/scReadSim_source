#!/bin/bash

chr_title=chr1
filename=BoneMarrow_62016_${chr_title}
directory=/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas
INPUT_noCB_bamfile=${directory}/${filename}.bam  
INPUT_cells_barcode_file=/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas/cell_metadata_bonemarrow.csv

# Tools directory
script_directory=/home/gayan/Projects/scATAC_Simulator/scripts/03092022/sciATAC/ATAC-seq
samtools_directory=/home/gayan/Tools/samtools/bin
macs3_directory=/home/gayan/.local/bin
bedtools_directory=/home/gayan/Tools/bedtools/bedtools2/bin
export PATH=${macs3_directory}:${samtools_directory}:${bedtools_directory}:${PATH}

# Output directory
outdirectory=/home/gayan/Projects/scATAC_Simulator/results/20230310_SplitCellType
mkdir ${outdirectory}
cd ${outdirectory}

####################################################################################
############################  Examine scReadSim per-cell-type bam file #######################
####################################################################################
# Python start here
##  process cell type file
import pandas as pd
cells = pd.read_csv("/home/gayan/Projects/scATAC_Simulator/results/20221130_BoneMarrow_62016_chr1_NONINPUT/synthetic_cell_barcode.txt.withSynthCluster", delimiter="\t",  names=['cell', 'cell_type'])
cells["cell_type_splitted"] = cells["cell_type"].str.split(".", expand=True)[0]
cells["cell"] = "CB:Z:" + cells["cell"].astype(str)
##  Generate cellbarcode file for four cell types from synthetic cells
outdirectory="/home/gayan/Projects/scATAC_Simulator/results/20230310_SplitCellType"
cell_type_selection_list = ["Hematopoietic progenitors", "Erythroblasts", "Monocytes", "Immature B cells"]
## Write out cell barcode file per cell type
for cell_type in cell_type_selection_list:
    cells.loc[cells["cell_type_splitted"]==cell_type][['cell', 'cell_type_splitted']].to_csv(outdirectory + "/CBwithCellType.%s.scReadSim.txt" % cell_type.replace(" ",""), sep="\t", header=False, index=False)
    cells.loc[cells["cell_type_splitted"]==cell_type][['cell']].to_csv(outdirectory + "/CB.%s.scReadSim.txt" % cell_type.replace(" ",""), sep="\t", header=False, index=False)
# Python end here

# Subset scReadSim synthetic bam into per celltype (f)
scReadSim_bamfile=/home/gayan/Projects/scATAC_Simulator/results/20230304_SCANATACSim/BoneMarrow_62016_chr1.syntheticBAM.CBincluded.synthetic.chr1.ErrorIncluded.bam
for cell_type in Hematopoieticprogenitors Erythroblasts Monocytes ImmatureBcells
do
    celltype_cellbarcode_file=CB.${cell_type}.scReadSim.txt
    samtools view -H $scReadSim_bamfile > ${outdirectory}/SAM_header_${cell_type}
    # Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
    samtools view $scReadSim_bamfile | LC_ALL=C grep -F -f ${outdirectory}/${celltype_cellbarcode_file} > ${outdirectory}/filtered_SAM_body_${cell_type}
    # Combine header and body
    cat ${outdirectory}/SAM_header_${cell_type} ${outdirectory}/filtered_SAM_body_${cell_type} | samtools view -bS - > ${outdirectory}/${cell_type}.scReadSim.bam
    rm ${outdirectory}/SAM_header_${cell_type} ${outdirectory}/filtered_SAM_body_${cell_type}
    # samtools index
    samtools sort ${outdirectory}/${cell_type}.scReadSim.bam > ${outdirectory}/${cell_type}.scReadSim.sorted.bam
    samtools index ${outdirectory}/${cell_type}.scReadSim.sorted.bam
done

####################################################################################
#####################  Task 2: Per-Cell-type Ground Truth Peaks #####################
####################################################################################
# Python script: Fig2b_ATAC_SplitCellType_20230310.py