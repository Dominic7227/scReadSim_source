#!/bin/bash
full_filename=e18_mouse_brain_fresh_5k_gex_possorted_bam
full_directory=/home/gayan/Projects/scATAC_Simulator/data/10X_MOUSE_BRAIN_ATACandRNA

referenceGenome_dir=/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE
referenceGenome_name=GRCm38.primary_assembly.genome
referenceGenome_file=${referenceGenome_dir}/${referenceGenome_name}.fa

chr_title=chr1
filename=${full_filename}_${chr_title}
outdirectory=/home/gayan/Projects/scATAC_Simulator/results/20221013_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_VaryNewCell

alevin_verification_dir=verification_alevin
mkdir ${alevin_verification_dir}
ref_dir=/home/gayan/Projects/scATAC_Simulator/results/20220608_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/verification_alevin_20220930/ref_data
export PATH=/home/gayan/Tools/salmon-1.8.0/bin:$PATH


BED_filename_combined_pre=${filename}.syntheticBAM.combined.20221011_UMITransLevel2

###################################################################
######################## VarySeqDepth ################################
###################################################################
timer_alevin_VarySeqDepth (){
local adj_factor=$1
BED_filename_pre=${filename}.syntheticBAM.VarySeqDepth${adj_factor}
BED_COMPLE_filename_pre=${filename}.syntheticBAM.VarySeqDepth${adj_factor}.COMPLE
BED_filename_combined_pre=${filename}.syntheticBAM.VarySeqDepth${adj_factor}.combined
output_BAM_pre=${BED_filename_pre}
synthetic_cell_barcode_file=synthetic_cell_barcode.VarySeqDepth${adj_factor}.txt

synthetic_read1_fq=${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq
synthetic_read2_fq=${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq

verification_dir=${outdirectory}/${alevin_verification_dir}/Case_VarySeqDepth${adj_factor}
mkdir ${verification_dir}

START=$(date +%s)
# Alvein count
salmon alevin -l ISR -1 $synthetic_read1_fq -2 $synthetic_read2_fq --chromium  -p 10 -o ${verification_dir}  -i ${ref_dir}/transcripts_index_tutorial --tgMap ${ref_dir}/txp2gene.tsv --whitelist ${outdirectory}/${synthetic_cell_barcode_file} --dumpMtx
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "[Task VarySeq ${adj_factor}] $DIFF seconds"
echo "[Task VarySeq ${adj_factor}] Lapsed $DIFF seconds" >> ${outdirectory}/verification_alevin/TIMER_VarySeqDepth.txt
}          
echo $(date) > ${outdirectory}/verification_alevin/TIMER_VarySeqDepth.txt
for adj_factor in 0.25 0.5 1 2 4
do 
     timer_alevin_VarySeqDepth "$adj_factor"
done


###################################################################
######################## VaryCell ################################
###################################################################
timer_alevin_VaryCell (){
local adj_factor=$1
BED_filename_pre=${filename}.syntheticBAM.VaryCellNumber${adj_factor}
BED_COMPLE_filename_pre=${filename}.syntheticBAM.VaryCellNumber${adj_factor}.COMPLE
BED_filename_combined_pre=${filename}.syntheticBAM.VaryCellNumber${adj_factor}.combined
output_BAM_pre=${BED_filename_pre}
synthetic_cell_barcode_file=synthetic_cell_barcode.VaryCellNumber${adj_factor}.txt

synthetic_read1_fq=${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq
synthetic_read2_fq=${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq

verification_dir=${outdirectory}/${alevin_verification_dir}/Case_VaryCell${adj_factor}
mkdir ${verification_dir}

START=$(date +%s)
# Alvein count
salmon alevin -l ISR -1 $synthetic_read1_fq -2 $synthetic_read2_fq --chromium  -p 10 -o ${verification_dir}  -i ${ref_dir}/transcripts_index_tutorial --tgMap ${ref_dir}/txp2gene.tsv --whitelist ${outdirectory}/${synthetic_cell_barcode_file} --dumpMtx
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "[Task VaryCell ${adj_factor}] $DIFF seconds"
echo "[Task VaryCell ${adj_factor}] Lapsed $DIFF seconds" >> ${outdirectory}/verification_alevin/TIMER_VaryCell.txt
}          
echo $(date) > ${outdirectory}/verification_alevin/TIMER_VaryCell.txt
for adj_factor in 0.25 0.5 1 2 4
do 
     timer_alevin_VaryCell "$adj_factor"
done