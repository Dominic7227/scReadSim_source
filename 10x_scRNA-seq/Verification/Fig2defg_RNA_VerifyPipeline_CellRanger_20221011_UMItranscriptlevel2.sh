#!/bin/bash
full_filename=e18_mouse_brain_fresh_5k_gex_possorted_bam
full_directory=/home/gayan/Projects/scATAC_Simulator/data/10X_MOUSE_BRAIN_ATACandRNA

referenceGenome_name=GRCm38.primary_assembly.genome
referenceGenome_dir=/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE
referenceGenome_file=${referenceGenome_dir}/${referenceGenome_name}.fa

chr_title=chr1
filename=${full_filename}_${chr_title}
outdirectory=/home/gayan/Projects/scATAC_Simulator/results/20221011_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT
BED_filename_combined_pre=${filename}.syntheticBAM.combined.20221011_UMITransLevel2
synthetic_read1_fq=${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq
synthetic_read2_fq=${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq

verification_dir=${outdirectory}/verification_cellranger
mkdir ${verification_dir}
export PATH=/home/gayan/Tools/cellranger-7.0.0:$PATH

mkdir ${outdirectory}/verification_cellranger_fastqs
cp $synthetic_read1_fq ${outdirectory}/verification_cellranger_fastqs/SampleName_S1_L001_R1_001.fastq
cp $synthetic_read2_fq ${outdirectory}/verification_cellranger_fastqs/SampleName_S1_L001_R2_001.fastq

# Note: 
# Add cell barcode to cellranger whitelist: /home/gayan/Tools/cellranger-7.0.0/lib/python/cellranger/barcodes/737K-august-2016.txt for 3 Prime V2 chemistr, which UMI is 10 bp
## Step1: Cellranger count
# rm -r verification_cellranger
# cellranger count --id=verification_cellranger_20220629 \
# --fastqs=${outdirectory}/verification_cellranger_fastqs \
# --sample=SampleName \
# --expect-cells 4877 \
# --transcriptome=/home/gayan/Tools/cellranger-7.0.0/refdata-gex-mm10-2020-A \
# --chemistry=SC3Pv2 

######################## 20221009 ################################
cd ${verification_dir}
rm -r cellranger_results
cellranger count --id=cellranger_results \
--fastqs=${outdirectory}/verification_cellranger_fastqs \
--sample=SampleName \
--expect-cells 4877 \
--transcriptome=/home/gayan/Projects/scATAC_Simulator/results/20220608_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/verification_cellranger_20220926/ref_data/mm10_cellranger_genome \
--chemistry=SC3Pv2 
