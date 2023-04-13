#!/bin/bash
full_filename=e18_mouse_brain_fresh_5k_gex_possorted_bam
full_directory=/home/gayan/Projects/scATAC_Simulator/data/10X_MOUSE_BRAIN_ATACandRNA

referenceGenome_name=GRCm38.primary_assembly.genome
referenceGenome_dir=/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE
referenceGenome_file=${referenceGenome_dir}/${referenceGenome_name}.fa

chr_title=chr1
filename=${full_filename}_${chr_title}
outdirectory=/home/gayan/Projects/scATAC_Simulator/results/20221013_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_VaryNewCell

###################################################################
######################## VarySeqDepth ################################
###################################################################
timer_CellRanger_VarySeqDepth (){
local adj_factor=$1
BED_filename_pre=${filename}.syntheticBAM.VarySeqDepth${adj_factor}
BED_COMPLE_filename_pre=${filename}.syntheticBAM.VarySeqDepth${adj_factor}.COMPLE
BED_filename_combined_pre=${filename}.syntheticBAM.VarySeqDepth${adj_factor}.combined
output_BAM_pre=${BED_filename_pre}

synthetic_read1_fq=${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq
synthetic_read2_fq=${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq

verification_dir=${outdirectory}/verification_cellranger/Case_VarySeqDepth${adj_factor}
mkdir ${verification_dir}
export PATH=/home/gayan/Tools/cellranger-7.0.0:$PATH

######################## Copy FQs ################################
fastq_dir=verification_cellranger_fastqs
mkdir ${verification_dir}/${fastq_dir}
# cp $synthetic_read1_fq ${verification_dir}/${fastq_dir}/SampleName_S1_L001_R1_001.fastq
# cp $synthetic_read2_fq ${verification_dir}/${fastq_dir}/SampleName_S1_L001_R2_001.fastq

START=$(date +%s)
######################## Run Cellranger ################################
cd ${verification_dir}
rm -r cellranger_results
cellranger count --id=cellranger_results \
--fastqs=${verification_dir}/verification_cellranger_fastqs \
--sample=SampleName \
--expect-cells 4877 \
--transcriptome=/home/gayan/Projects/scATAC_Simulator/results/20220608_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/verification_cellranger_20220926/ref_data/mm10_cellranger_genome \
--chemistry=SC3Pv2 
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "[Task VarySeq ${adj_factor}] $DIFF seconds"
echo "[Task VarySeq ${adj_factor}] Lapsed $DIFF seconds" >> ${outdirectory}/verification_cellranger/TIMER_VarySeqDepth.txt
}
# for adj_factor in 0.125 0.25 0.5 1 2 4 8; do GeneateBAM "$adj_factor" & done
echo $(date) >> ${outdirectory}/verification_cellranger/TIMER_VarySeqDepth.txt
for adj_factor in 0.25 0.5 1 2 4
do
	timer_CellRanger_VarySeqDepth "$adj_factor"  
done

# Sat Oct 15 19:35:45 UTC 2022
# [Task VarySeq 0.25] Lapsed 12274 seconds
# [Task VarySeq 0.5] Lapsed 3005 seconds
# [Task VarySeq 1] Lapsed 4532 seconds
# [Task VarySeq 2] Lapsed 9284 seconds
# [Task VarySeq 4] Lapsed 9401 seconds
# Sun Oct 16 16:38:01 UTC 2022
# [Task VarySeq 0.25] Lapsed 1671 seconds
# [Task VarySeq 0.5] Lapsed 2235 seconds
# [Task VarySeq 1] Lapsed 5118 seconds
# [Task VarySeq 2] Lapsed 6722 seconds
# [Task VarySeq 4] Lapsed 7007 seconds


###################################################################
######################## VaryCell ################################
###################################################################
timer_CellRanger_VaryCell (){
	local adj_factor=$1
  	BED_filename_pre=${filename}.syntheticBAM.VaryCellNumber${adj_factor}
  	BED_COMPLE_filename_pre=${filename}.syntheticBAM.VaryCellNumber${adj_factor}.COMPLE
  	BED_filename_combined_pre=${filename}.syntheticBAM.VaryCellNumber${adj_factor}.combined
	output_BAM_pre=${BED_filename_pre}

	synthetic_read1_fq=${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq
	synthetic_read2_fq=${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq

	verification_dir=${outdirectory}/verification_cellranger/Case_VaryCell${adj_factor}
	mkdir ${verification_dir}
	export PATH=/home/gayan/Tools/cellranger-7.0.0:$PATH

	######################## Copy FQs ################################
	fastq_dir=verification_cellranger_fastqs
	mkdir ${verification_dir}/${fastq_dir}
	# cp $synthetic_read1_fq ${verification_dir}/${fastq_dir}/SampleName_S1_L001_R1_001.fastq
	# cp $synthetic_read2_fq ${verification_dir}/${fastq_dir}/SampleName_S1_L001_R2_001.fastq

	START=$(date +%s)
	######################## Run Cellranger ################################
	cd ${verification_dir}
	rm -r cellranger_results
	cellranger count --id=cellranger_results \
	--fastqs=${verification_dir}/verification_cellranger_fastqs \
	--sample=SampleName \
	--expect-cells 4877 \
	--transcriptome=/home/gayan/Projects/scATAC_Simulator/results/20220608_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/verification_cellranger_20220926/ref_data/mm10_cellranger_genome \
	--chemistry=SC3Pv2
	END=$(date +%s)
	DIFF=$(( $END - $START ))
	echo "[Task VaryCell ${adj_factor}] $DIFF seconds"
	echo "[Task VaryCell ${adj_factor}] Lapsed $DIFF seconds" >> ${outdirectory}/verification_cellranger/TIMER_VaryCell.txt
}
# for adj_factor in 0.125 0.25 0.5 1 2 4 8; do GeneateBAM "$adj_factor" & done
echo $(date) >> ${outdirectory}/verification_cellranger/TIMER_VaryCell.txt
for adj_factor in 0.25 0.5 1 2 4
do 
	timer_CellRanger_VaryCell "$adj_factor" 
done




# -rw-rw-r-- 1 gayan gayan  6882670516 Oct 14 19:26 e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.syntheticBAM.VaryCellNumber0.25.combined.read2.bed2fa.sorted.fq
# -rw-rw-r-- 1 gayan gayan  7181978174 Oct 14 20:14 e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.syntheticBAM.VaryCellNumber0.5.combined.read2.bed2fa.sorted.fq
# -rw-rw-r-- 1 gayan gayan  7155522338 Oct 14 20:14 e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.syntheticBAM.VaryCellNumber1.combined.read2.bed2fa.sorted.fq
# -rw-rw-r-- 1 gayan gayan  7065327961 Oct 14 20:14 e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.syntheticBAM.VaryCellNumber2.combined.read2.bed2fa.sorted.fq
# -rw-rw-r-- 1 gayan gayan  7104459958 Oct 14 20:14 e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.syntheticBAM.VaryCellNumber4.combined.read2.bed2fa.sorted.fq