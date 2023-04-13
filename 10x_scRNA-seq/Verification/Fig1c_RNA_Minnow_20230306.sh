#!/bin/bash

full_filename=e18_mouse_brain_fresh_5k_gex_possorted_bam
full_directory=/home/gayan/Projects/scATAC_Simulator/data/10X_MOUSE_BRAIN_ATACandRNA
INPUT_bamfile=${full_directory}/split.${full_filename}/${full_filename}_chr1.bam

referenceGenome_dir=/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE
referenceGenome_name=GRCm38.primary_assembly.genome
referenceGenome_file=${referenceGenome_dir}/${referenceGenome_name}.fa

chr_title=chr1
filename=${full_filename}_${chr_title}
outdirectory=/home/gayan/Projects/scATAC_Simulator/results/20230304_Minnow
mkdir ${outdirectory}

# Run Minnow
export PATH=/home/gayan/Tools/salmon-1.8.0/bin:$PATH
export PATH=/home/gayan/Tools/minnow/minnow-latest_linux_x86_64/bin:$PATH
export LD_LIBRARY_PATH=/home/gayan/Tools/minnow/minnow-latest_linux_x86_64/lib:$LD_LIBRARY_PATH

# Step 0: Salmon index and Alevin generating Bfh file
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M10/gencode.vM10.transcripts.fa.gz
gunzip gencode.vM10.transcripts.fa.gz
grep ">" ${referenceGenome_file} | cut -d ">" -f 2 | cut -d " " -f 1 > ${outdirectory}/GRCm38.primary_assembly.genome.chrnames.txt
salmon index \
-t ${outdirectory}/gencode.vM10.pc_transcripts.fa \
-i ${outdirectory}/gencode.vM10.annotation.sidx --gencode -p 10 
# -d ${outdirectory}/GRCm38.primary_assembly.genome.chrnames.txt

# Alvein count a dataset from the same organism
# Real BAM file to fastq (using bamtofastq from 10x)
/home/gayan/Tools/bam2fastq/bamtofastq-1.4.1/target/release/bamtofastq --reads-per-fastq=500000000 ${INPUT_bamfile} ${outdirectory}/FASTQtmp
# Construct a geneid to transcript id tab-separated file. Col 1 is transcript id, col 2 is gene id
awk -F "\t" '$3 == "transcript" { print $9 }' /home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE/gencode.vM10.annotation.gtf | tr -d ";\"" | awk '{print $4"\t"$2}' > ${outdirectory}/gencode.vM10.annotation.tx2gene.tsv
# Generate BFH file using Alevin
salmon alevin -l ISR -i ${outdirectory}/gencode.vM10.annotation.sidx \
    -1 ${outdirectory}/FASTQtmp/e18_mouse_brain_fresh_5k_0_1_HGWCKDSXY/bamtofastq_S1_L001_R1_001.fastq.gz \
    -2 ${outdirectory}/FASTQtmp/e18_mouse_brain_fresh_5k_0_1_HGWCKDSXY/bamtofastq_S1_L001_R2_001.fastq.gz \
    -o alevin_out -p 36 \
    --tgMap ${outdirectory}/gencode.vM10.annotation.tx2gene.tsv \
    --chromium \
    --dumpFeatures --expectCells 4877 \
    --dumpBfh

# Step 1: Minnow index
minnow index -r ${outdirectory}/gencode.vM10.pc_transcripts.fa -k 101 -f 20 --tmpdir tmp -p 10 -o minnow_index
# Step 2: Minnow count
minnow estimate -o minnow_estimate -r ${outdirectory}/gencode.vM10.pc_transcripts.fa --g2t ${outdirectory}/gencode.vM10.annotation.tx2gene.tsv --bfh ${outdirectory}/alevin_out/alevin/bfh.txt
# Step 3: Prepare count matrix 
# Synthetic count matrix
mkdir ${outdirectory}/ground_truth_scReadSim_countmat
cp /home/gayan/Projects/scATAC_Simulator/results/20221011_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/synthetic_cell_barcode.txt ${outdirectory}/ground_truth_scReadSim_countmat/quants_mat_cols.txt
ground_truth_count_matrix=${outdirectory}/ground_truth_scReadSim_countmat/quants_mat.txt # Generate usign Rscript NOINPUT_RNA_Minnow_PrepareCountMatrix.R
ground_truth_count_matrix_row=${outdirectory}/ground_truth_scReadSim_countmat/quants_mat_rows.txt # Generate usign Rscript NOINPUT_RNA_Minnow_PrepareCountMatrix.R
ground_truth_count_matrix_col=${outdirectory}/ground_truth_scReadSim_countmat/quants_mat_cols.txt

# Step 4: Simulate reads in FASTQs 
minnow simulate --splatter-mode \
    --g2t  ${outdirectory}/gencode.vM10.annotation.tx2gene.tsv \
    --inputdir ${outdirectory}/ground_truth_scReadSim_countmat \
    --PCR 4 \
    -r ${outdirectory}/minnow_index/ref_k101_fixed.fa \
    -e 0.01 \
    -p 16 \
    -o ${outdirectory}/minnow_simulate \
    --dbg \
    --gfa ${outdirectory}/minnow_index/dbg.gfa \
    -w ${outdirectory}/737K-august-2016.txt \
    --countProb ${outdirectory}/minnow_estimate/countProb.txt \
    --custom \
    --gencode 

# Step 5: Read alignment
export PATH=/home/gayan/Tools/bowtie2/bowtie2-2.3.4.1-linux-x86_64:$PATH
script_directory=/home/gayan/Projects/scATAC_Simulator/scripts/03092022/sciATAC/ATAC-seq
samtools_directory=/home/gayan/Tools/samtools/bin
macs3_directory=/home/gayan/.local/bin
bedtools_directory=/home/gayan/Tools/bedtools/bedtools2/bin
export PATH=${macs3_directory}:${samtools_directory}:${bedtools_directory}:${PATH}

time(bowtie2 -x ${referenceGenome_dir}/${referenceGenome_name} -U ${outdirectory}/minnow_simulate/hg_100_S1_L001_R2_001.fastq.gz -S ${outdirectory}/minnow_simulate/hg_100_S1_L001_R2_001.sam) # 151m37.660s
samtools view -S -b ${outdirectory}/minnow_simulate/hg_100_S1_L001_R2_001.sam > ${outdirectory}/minnow_simulate/hg_100_S1_L001_R2_001.bam
# Convert sam to bam
samtools sort ${outdirectory}/minnow_simulate/hg_100_S1_L001_R2_001.bam > ${outdirectory}/minnow_simulate/hg_100_S1_L001_R2_001.sorted.bam
samtools index ${outdirectory}/minnow_simulate/hg_100_S1_L001_R2_001.sorted.bam


################# Plot read coverage #################
# Python script: Fig1c_RNA_Minnow_PlotCoverage_20230316.py