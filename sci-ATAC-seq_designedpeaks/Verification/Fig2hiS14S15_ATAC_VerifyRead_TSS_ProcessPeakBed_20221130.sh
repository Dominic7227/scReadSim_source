## TSS
filename=BoneMarrow_62016_chr1
outdirectory=~/Projects/scATAC_Simulator/results/20221130_${filename}_INPUT_withCluster_mm9TSS
chr_title=chr1
scReadSim_bamfile=${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam

script_directory=/home/gayan/Projects/scATAC_Simulator/scripts/03092022/sciATAC/INPUT_ATAC-seq/WithCB_scDesign2forsynthetic
samtools_directory=/home/gayan/Tools/samtools/bin
macs3_directory=/home/gayan/.local/bin
bedtools_directory=/home/gayan/Tools/bedtools/bedtools2/bin
export PATH=${macs3_directory}:${samtools_directory}:${bedtools_directory}:${PATH}

########################## With stringent rules ##########################
# MACS3
synthetic_MACS3_peakcall_filename=${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.MACS3
# Verification purpose: verify the peak calling results
bedtools sort -i ${outdirectory}/${synthetic_MACS3_peakcall_filename}_peaks.narrowPeak | ${bedtools_directory}/bedtools merge  > ${outdirectory}/${synthetic_MACS3_peakcall_filename}_peaks.narrowPeak.bed 
python ${script_directory}/INPUT_ATAC_ComplePeakFunction.py ${chr_title} ${outdirectory} ${synthetic_MACS3_peakcall_filename}_peaks.narrowPeak.bed  ${outdirectory} ${synthetic_MACS3_peakcall_filename}_peaks.narrowPeak.COMPLE.bed
# Synthetic BAM read coverage over MACS3 peaks and non-peaks
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} ${synthetic_MACS3_peakcall_filename}_peaks.narrowPeak.bed ${scReadSim_bamfile} ${count_mat_filename}.scDesign2Simulated.MACS3Peaks.MarginalFeatureCount.txt
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} ${synthetic_MACS3_peakcall_filename}_peaks.narrowPeak.COMPLE.bed ${scReadSim_bamfile} ${count_mat_comple_filename}.scDesign2Simulated.MACS3Peaks.MarginalFeatureCount.txt


# SEACR
mv  ${outdirectory}/BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.searc.stringent.bed ${outdirectory}/BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.SEACR_peaks.narrowPeak.bed
python3 ${script_directory}/INPUT_ATAC_ComplePeakFunction.py chr1 ${outdirectory} BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.SEACR_peaks.narrowPeak.bed ${outdirectory} BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.SEACR_peaks.narrowPeak.COMPLE.bed
# Synthetic BAM read coverage over SEACR peaks and non-peaks
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.SEACR_peaks.narrowPeak.bed ${scReadSim_bamfile} ${count_mat_filename}.scDesign2Simulated.SEACRPeaks.MarginalFeatureCount.txt
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.SEACR_peaks.narrowPeak.COMPLE.bed ${scReadSim_bamfile} ${count_mat_comple_filename}.scDesign2Simulated.SEACRPeaks.MarginalFeatureCount.txt


# HOMER
python3 ${script_directory}/INPUT_ATAC_ComplePeakFunction.py chr1 ${outdirectory} BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HOMER.merged.bed ${outdirectory} BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HOMER.merged.COMPLE.bed
# Synthetic BAM read coverage over HOMER peaks and non-peaks
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HOMER.merged.bed ${scReadSim_bamfile} ${count_mat_filename}.scDesign2Simulated.HOMERPeaks.MarginalFeatureCount.txt
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HOMER.merged.COMPLE.bed ${scReadSim_bamfile} ${count_mat_comple_filename}.scDesign2Simulated.HOMERPeaks.MarginalFeatureCount.txt


# HMMRATAC
python3 ${script_directory}/INPUT_ATAC_ComplePeakFunction.py chr1 ${outdirectory} BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HMMRATAC.bed ${outdirectory} BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HMMRATAC.COMPLE.bed
# Synthetic BAM read coverage over HMMRATAC peaks and non-peaks
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HMMRATAC.bed ${scReadSim_bamfile} ${count_mat_filename}.scDesign2Simulated.HMMRATACPeaks.MarginalFeatureCount.txt
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HMMRATAC.COMPLE.bed ${scReadSim_bamfile} ${count_mat_comple_filename}.scDesign2Simulated.HMMRATACPeaks.MarginalFeatureCount.txt

################## Make plots based on peak calling sets ##################
mkdir ${outdirectory}/INTERVENE_20221130_noSEACR
mkdir ${outdirectory}/INTERVENE_20221130_withSEACR
intervene venn -i BoneMarrow_62016_chr1.INPUT.peaks.bed BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.MACS3_peaks.narrowPeak.bed BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.SEACR_peaks.narrowPeak.bed BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HOMER.merged.bed BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HMMRATAC.bed -o ${outdirectory}/INTERVENE_20221130_withSEACR --names INPUT,MACS3,SEACR,HOMER,HMMRATAC --fontsize 20
intervene upset -i BoneMarrow_62016_chr1.INPUT.peaks.bed BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.MACS3_peaks.narrowPeak.bed BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.SEACR_peaks.narrowPeak.bed BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HOMER.merged.bed BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HMMRATAC.bed -o ${outdirectory}/INTERVENE_20221130_withSEACR --names INPUT,MACS3,SEACR,HOMER,HMMRATAC 
# no SEACR
intervene venn -i BoneMarrow_62016_chr1.INPUT.peaks.bed BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.MACS3_peaks.narrowPeak.bed BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HOMER.merged.bed BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HMMRATAC.bed -o ${outdirectory}/INTERVENE_20221130_noSEACR --names INPUT,MACS3,HOMER,HMMRATAC --fontsize 20
intervene upset -i BoneMarrow_62016_chr1.INPUT.peaks.bed BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.MACS3_peaks.narrowPeak.bed BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HOMER.merged.bed BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HMMRATAC.bed -o ${outdirectory}/INTERVENE_20221130_noSEACR --names INPUT,MACS3,HOMER,HMMRATAC 
