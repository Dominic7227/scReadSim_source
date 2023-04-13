#!/bin/bash
full_filename=e18_mouse_brain_fresh_5k_gex_possorted_bam
full_directory=/home/gayan/Projects/scATAC_Simulator/data/10X_MOUSE_BRAIN_ATACandRNA

referenceGenome_dir=/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE
referenceGenome_name=GRCm38.primary_assembly.genome
referenceGenome_file=${referenceGenome_dir}/${referenceGenome_name}.fa

chr_title=chr1
filename=${full_filename}_${chr_title}
outdirectory=/home/gayan/Projects/scATAC_Simulator/results/20221011_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT
BED_filename_combined_pre=${filename}.syntheticBAM.combined.20221011_UMITransLevel2
synthetic_read1_fq=${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq
synthetic_read2_fq=${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq


######################## 20221009 ################################
# Tutorial: https://combine-lab.github.io/alevin-tutorial/2018/setting-up-resources/
verification_dir=verification_alevin
mkdir ${verification_dir}
ref_dir=/home/gayan/Projects/scATAC_Simulator/results/20220608_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/verification_alevin_20220930/ref_data
export PATH=/home/gayan/Tools/salmon-1.8.0/bin:$PATH
# Step3: Alvein count
salmon alevin -l ISR -1 $synthetic_read1_fq -2 $synthetic_read2_fq --chromium  -p 10 -o ${outdirectory}/${verification_dir}  -i ${ref_dir}/transcripts_index_tutorial --tgMap ${ref_dir}/txp2gene.tsv --whitelist ${outdirectory}/synthetic_cell_barcode.txt --dumpMtx
cd ${outdirectory}/${verification_dir}/alevin
gzip -d quants_mat.mtx.gz


