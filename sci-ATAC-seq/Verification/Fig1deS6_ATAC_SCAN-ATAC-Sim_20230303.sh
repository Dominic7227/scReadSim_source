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
outdirectory=/home/gayan/Projects/scATAC_Simulator/results/20230304_SCANATACSim
mkdir ${outdirectory}
cd ${outdirectory}

####################################################################################
############################ Prepare scReadSim BAM file ############################
####################################################################################
scReadSim_noCB_bamfile=/home/gayan/Projects/scATAC_Simulator/results/20221130_BoneMarrow_62016_chr1_NONINPUT/BoneMarrow_62016_chr1.syntheticBAM.CBincluded.synthetic.chr1.noCB.ErrorIncluded.bam
scReadSim_bamfile=${outdirectory}/BoneMarrow_62016_chr1.syntheticBAM.CBincluded.synthetic.chr1.ErrorIncluded.bam
scReadSim_cells_barcode_file=/home/gayan/Projects/scATAC_Simulator/results/20221130_BoneMarrow_62016_chr1_NONINPUT/synthetic_cell_barcode.txt.withSynthCluster
# Regenerate CB tag for noCB BAM file (20221130 version CB not added)
time(samtools view $scReadSim_noCB_bamfile -H > ${outdirectory}/BoneMarrow_62016_chr1.syntheticBAM.CBincluded.synthetic.chr1.noCB.header.sam)
time(cat <( cat ${outdirectory}/BoneMarrow_62016_chr1.syntheticBAM.CBincluded.synthetic.chr1.noCB.header.sam ) \
<( paste <(samtools view $scReadSim_noCB_bamfile ) <(samtools view $scReadSim_noCB_bamfile | cut -f1 | cut -d':' -f1 | cut -c-16 | sed -e 's/^/CB:Z:/')) | samtools view -bS - > ${outdirectory}/BoneMarrow_62016_chr1.syntheticBAM.CBincluded.synthetic.chr1.ErrorIncluded.bam)
rm ${outdirectory}/BoneMarrow_62016_chr1.syntheticBAM.CBincluded.synthetic.chr1.noCB.header.sam

# Python start here
##  process cell type file
import pandas as pd
cells = pd.read_csv("/home/gayan/Projects/scATAC_Simulator/results/20221130_BoneMarrow_62016_chr1_NONINPUT/synthetic_cell_barcode.txt.withSynthCluster", delimiter="\t",  names=['cell', 'cell_type'])
cells["cell_type_splitted"] = cells["cell_type"].str.split(".", expand=True)[0]
cells["cell"] = "CB:Z:" + cells["cell"].astype(str)
## Cells types 
# Hematopoietic progenitors    1621 *
# Erythroblasts                1153 *
# Monocytes                     531 *
# Immature B cells              296 *
# Collisions                    102
# Unknown                        91
# Macrophages                    69
# B cells                        58
# Dendritic cells                50
# T cells                        22
# Regulatory T cells             17
# NK cells                       17
# Activated B cells               6
##  Generate cellbarcode file for four cell types from synthetic cells
outdirectory="/home/gayan/Projects/scATAC_Simulator/results/20230304_SCANATACSim"
cell_type_selection_list = ["Hematopoietic progenitors", "Erythroblasts", "Monocytes", "Immature B cells"]
## Write out cell barcode file
cells.loc[cells["cell_type_splitted"].isin(cell_type_selection_list)][['cell', 'cell_type_splitted']].to_csv(outdirectory + "/CBwithCellType.HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.scReadSim.txt", sep="\t", header=False, index=False)
cells.loc[cells["cell_type_splitted"].isin(cell_type_selection_list)][['cell']].to_csv(outdirectory + "/CB.HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.scReadSim.txt", sep="\t", header=False, index=False)
# Python end here

# Subset scReadSim bam file into four cell types
# Save the header lines
samtools view -H $scReadSim_bamfile > ${outdirectory}/SAM_header_scReadSim
# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
samtools view $scReadSim_bamfile | LC_ALL=C grep -F -f ${outdirectory}/CB.HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.scReadSim.txt > ${outdirectory}/filtered_SAM_body_scReadSim
# Combine header and body
cat ${outdirectory}/SAM_header_scReadSim ${outdirectory}/filtered_SAM_body_scReadSim | samtools view -bS - > ${outdirectory}/HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.scReadSim.bam
rm ${outdirectory}/SAM_header_scReadSim ${outdirectory}/filtered_SAM_body_scReadSim
# Attach CB barcode to the header of reads. This is for scReadSim.Utility.scATAC_bam2countmat_paral
mkdir tmp
samtools view ${outdirectory}/HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.scReadSim.bam -H > tmp/header.sam
time(cat <( cat tmp/header.sam ) \
 <( samtools view ${outdirectory}/HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.scReadSim.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) \
 | samtools view -bS - > ${outdirectory}/HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.scReadSim.forReadCounting.bam)
rm -dr tmp
# samtools index
samtools sort ${outdirectory}/HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.scReadSim.forReadCounting.bam > ${outdirectory}/HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.scReadSim.forReadCounting.sorted.bam
samtools index ${outdirectory}/HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.scReadSim.forReadCounting.sorted.bam


######################################################################################
################################ Prepare real BAM file ###############################
######################################################################################
# Regenerate CB tag for real BAM file 
time(samtools view $INPUT_noCB_bamfile -H > ${outdirectory}/BoneMarrow_62016_chr1.header.sam)
time(cat <( cat ${outdirectory}/BoneMarrow_62016_chr1.header.sam) \
<( paste <(samtools view $INPUT_noCB_bamfile ) <(samtools view $INPUT_noCB_bamfile | cut -f1 | cut -d':' -f1 | cut -c-36 | sed -e 's/^/CB:Z:/')) | samtools view -bS - > ${outdirectory}/BoneMarrow_62016_chr1.CBincluded.bam)
rm ${outdirectory}/BoneMarrow_62016_chr1.header.sam
INPUT_bamfile=${outdirectory}/BoneMarrow_62016_chr1.CBincluded.bam

# Python start here
##  process cell type file
import pandas as pd
cells = pd.read_csv("/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas/cell_metadata_bonemarrow.csv", delimiter="\t")
cells["cell"] = "CB:Z:" + cells["cell"].astype(str)
# cells["cell_label"].value_counts() 
cell_type_selection_list = ["Hematopoietic progenitors", "Erythroblasts", "Monocytes", "Immature B cells"]
## Write out cell barcode file for four celltypes
cells.loc[cells["cell_label"].isin(cell_type_selection_list)][['cell', 'cell_label']].to_csv(outdirectory + "/CBwithCellType.HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.Real.txt", sep="\t", header=False, index=False)
cells.loc[cells["cell_label"].isin(cell_type_selection_list)][['cell']].to_csv(outdirectory + "/CB.HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.Real.txt", sep="\t", header=False, index=False)
## Write out cell barcode file per cell type
for cell_type in cell_type_selection_list:
    cells.loc[cells["cell_label"]==cell_type][['cell', 'cell_label']].to_csv(outdirectory + "/CBwithCellType.%s.Real.txt" % cell_type.replace(" ",""), sep="\t", header=False, index=False)
    cells.loc[cells["cell_label"]==cell_type][['cell']].to_csv(outdirectory + "/CB.%s.Real.txt" % cell_type.replace(" ",""), sep="\t", header=False, index=False)
# Python end here

# Subset real bam into four celltypes
# Save the header lines
samtools view -H $INPUT_bamfile > ${outdirectory}/SAM_header
# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
samtools view $INPUT_bamfile | LC_ALL=C grep -F -f ${outdirectory}/CB.HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.Real.txt > ${outdirectory}/filtered_SAM_body
# Combine header and body
cat ${outdirectory}/SAM_header ${outdirectory}/filtered_SAM_body | samtools view -bS - > ${outdirectory}/HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.Real.bam
rm ${outdirectory}/SAM_header ${outdirectory}/filtered_SAM_body
# samtools index
samtools sort ${outdirectory}/HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.Real.bam > ${outdirectory}/HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.Real.sorted.bam
samtools index ${outdirectory}/HematopoieticProgenitor_Erythroblasts_Monocytes_ImmatureB.Real.sorted.bam

# Subset real bam into per celltype (for scan-atac-sim)
for cell_type in Hematopoieticprogenitors Erythroblasts Monocytes ImmatureBcells
do
    celltype_cellbarcode_file=CB.${cell_type}.Real.txt
    samtools view -H $INPUT_bamfile > ${outdirectory}/SAM_header_${cell_type}
    # Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
    samtools view $INPUT_bamfile | LC_ALL=C grep -F -f ${outdirectory}/${celltype_cellbarcode_file} > ${outdirectory}/filtered_SAM_body_${cell_type}
    # Combine header and body
    cat ${outdirectory}/SAM_header_${cell_type} ${outdirectory}/filtered_SAM_body_${cell_type} | samtools view -bS - > ${outdirectory}/${cell_type}.Real.bam
    rm ${outdirectory}/SAM_header_${cell_type} ${outdirectory}/filtered_SAM_body_${cell_type}
    # samtools index
    samtools sort ${outdirectory}/${cell_type}.Real.bam > ${outdirectory}/${cell_type}.Real.sorted.bam
    samtools index ${outdirectory}/${cell_type}.Real.sorted.bam
done

###############################################################################
################################ SCAN-ATAC-Sim  ###############################
###############################################################################
# Call peaks
for cell_type in Hematopoieticprogenitors Erythroblasts Monocytes ImmatureBcells
do
    # Peak calling
    ${macs3_directory}/macs3 callpeak -f BAMPE -t ${outdirectory}/${cell_type}.Real.sorted.bam -g mm -n ${outdirectory}/${cell_type}.Real. -B -q 0.05
done

# Run SCAN-ATAC-Sim
# Prepare ./bam and ./peak directory for SCAN-ATAC-Sim
mkdir ${outdirectory}/SCANATACSim_bam
mkdir ${outdirectory}/SCANATACSim_peak
mkdir ${outdirectory}/SCANATACSim_output

for cell_type in Hematopoieticprogenitors Erythroblasts Monocytes ImmatureBcells
do
    cp ${outdirectory}/${cell_type}.Real.sorted.bam ${outdirectory}/SCANATACSim_bam/${cell_type}.bam
    cp ${outdirectory}/${cell_type}.Real._peaks.narrowPeak ${outdirectory}/SCANATACSim_peak/${cell_type}.narrowPeak
done

# Run SCAN-ATAC-Sim
SCANATACSim_dir=/home/gayan/Tools/SCAN-ATAC-Sim/finalized
# Creat conda envrionment
conda create -n SCAN-ATAC-Sim
conda activate SCAN-ATAC-Sim
conda install -c bioconda pysam=0.15.4 pybedtools=0.8.0 pybigwig=0.3.17
conda install pip # Use conda installed pip to install packages 
/home/gayan/miniconda3/envs/SCAN-ATAC-Sim/bin/pip install numpy==1.16.4 pandas==1.0.3 wget==3.2

# Preprocess bam and called peaks
time(python ${SCANATACSim_dir}/preprocess_edited.py \
    -c Hematopoieticprogenitors,Erythroblasts,Monocytes,ImmatureBcells \
    -i ${outdirectory}/SCANATACSim_peak/ \
    -j ${outdirectory}/SCANATACSim_bam/ \
    -e 1000 \
    -b 1000 \
    -o ${outdirectory}/SCANATACSim_output/)

temp_dir=${outdirectory}/SCANATACSim_output
array=( Hematopoieticprogenitors Erythroblasts Monocytes ImmatureBcells )
array2=( 1621 1153 531 296 )
for i in "${!array[@]}"
do
    cell_name=${array[i]}
    cell_number=${array2[i]}
    printf "Cell type %s: %s\n" "$cell_name" "$cell_number" 
    ${SCANATACSim_dir}/weighted_sampling -f ${temp_dir}/${cell_name}.peak_counts.bed -b ${temp_dir}/bg_counts.bed -of ${temp_dir}/${cell_name}.foreground.sampled.bed -ob ${temp_dir}/${cell_name}.background.sampled.bed -n 2000 -nv 0.5 -c $cell_number -s 0.6 -min 1000 -max 20000
    ${SCANATACSim_dir}/uniform_sampling ${temp_dir}/${cell_name}.peak_intersect.bed ${temp_dir}/${cell_name}.foreground.sampled.bed ${temp_dir}/${cell_name}.foreground.bed
    ${SCANATACSim_dir}/uniform_sampling ${temp_dir}/bg_intersect.expanded.bed ${temp_dir}/${cell_name}.background.sampled.bed ${temp_dir}/${cell_name}.background.bed
    cat ${temp_dir}/${cell_name}.foreground.bed ${temp_dir}/${cell_name}.background.bed > ${temp_dir}/${cell_name}.SCAN-ATAC-Sim.bed
done

# Generate Count matrix for SCAN-ATAC-Sim using Python script: Fig1deS6_ATAC_SCAN-ATAC-Sim_CountMat_20230310.py

# Comparison 2: read coverage
cat ${outdirectory}/${BED_filename_pre}.read1.bed ${outdirectory}/${BED_COMPLE_filename_pre}.read1.bed  > ${outdirectory}/${BED_filename_combined_pre}.read1.unsort.bed
