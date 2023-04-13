#!/bin/bash

# Task: mouse sci-ATAC-seq synthetic reads simulation
# Date: 20221130
# Notes
# 1. Generate ground-truth peaks (around 200bp) from Transcription Starting Site (TSS)
#    Use generated peaks and non-peaks as ground-truth peaks and non-peaks


#####################################################################
############################ Prepare the data #######################
#####################################################################
full_directory=/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas
full_filename=BoneMarrow_62016
chr_title=chr1
filename=BoneMarrow_62016_${chr_title}
directory=/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas
INPUT_bamfile=${directory}/${filename}.bam 
INPUT_cells_barcode_file=/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas/cell_metadata_bonemarrow.csv
OUTPUT_cells_barcode_file=${outdirectory}/synthetic_cell_barcode.txt

outdirectory=/home/gayan/Projects/scATAC_Simulator/results/20221130_BoneMarrow_62016_chr1_INPUT_withCluster_mm9TSS
# outdirectory=~/Projects/scATAC_Simulator/results/20220302_${filename}_INPUT_withCluster_mm9TSS
mkdir ${outdirectory}
cd ${outdirectory}

script_directory=/home/gayan/Projects/scATAC_Simulator/scripts/03092022/sciATAC/INPUT_ATAC-seq/WithCB_scDesign2forsynthetic
samtools_directory=/home/gayan/Tools/samtools/bin
macs3_directory=/home/gayan/.local/bin
bedtools_directory=/home/gayan/Tools/bedtools/bedtools2/bin
export PATH=${macs3_directory}:${samtools_directory}:${bedtools_directory}:${PATH}



######################## Prestep Attach barcodes to BAM file ########################
# extract the header file
mkdir ${full_directory}/tmp
${samtools_directory}/samtools view ${full_directory}/${full_filename}.bam -H > ${full_directory}/tmp/${full_filename}.header.sam

# create a bam file with the barcode embedded into the read name
time(cat <( cat ${full_directory}/tmp/${full_filename}.header.sam ) \
<( ${samtools_directory}/samtools view ${full_directory}/${full_filename}.bam | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) \
| ${samtools_directory}/samtools view -bS - > ${full_directory}/${full_filename}.scReadSim.bam) 
rm -d ${full_directory}/tmp/

######################## Split, Index and Sort BAM file ######################## 
${samtools_directory}/samtools sort ${full_directory}/${full_filename}.scReadSim.bam -o ${full_directory}/${full_filename}.scReadSim.sorted.bam
${samtools_directory}/samtools index ${full_directory}/${full_filename}.scReadSim.sorted.bam ${full_directory}/${full_filename}.scReadSim.sorted.bam.bai 

## Separate bam files
for (( counter=1; counter<20; counter++ ))
do
time(${samtools_directory}/samtools view -b ${full_directory}/${full_filename}.scReadSim.sorted.bam chr$counter > ${directory}/${full_filename}_chr${counter}.bam)
time(${samtools_directory}/samtools index ${directory}/${full_filename}_chr${counter}.bam)
printf "chr$counter\tDone\n"
done
printf "Done\n"

time(${samtools_directory}/samtools view -b ${full_directory}/${full_filename}.scReadSim.sorted.bam chrX > "${directory}/${full_filename}_chrX.bam")
time(${samtools_directory}/samtools index "${directory}/${full_filename}_chrX.bam")
printf "chrX\tDone\n"

time(${samtools_directory}/samtools view -b ${full_directory}/${full_filename}.scReadSim.sorted.bam chrY > "${directory}/${full_filename}_chrY.bam")
time(${samtools_directory}/samtools index "${directory}/${full_filename}_chrY.bam")
printf "chrY\tDone\n"

time(${samtools_directory}/samtools view -b ${full_directory}/${full_filename}.scReadSim.sorted.bam chrM > "${directory}/${full_filename}_chrM.bam")
time(${samtools_directory}/samtools index "${directory}/${full_filename}_chrM.bam")
printf "chrM\tDone\n"


#####################################################################
############################ Main Analysis #########################
#####################################################################
######################## User Specified Peaks Set ######################## 
## User input bed file
# Store as true-peak.bed
# Take complementary regions and stored as true-non-peak.bed 
cat /home/gayan/Projects/scATAC_Simulator/data/mm9_ref_genome_GECODE/gencode.vM1.annotation.gtf | grep -v '^#' | awk '$1=="chr1" && $3=="transcript"' | cut -f1,4,5 | sort -k1,1 -k2,2n -k3,3n > ${outdirectory}/mm9.transcripts.${chr_title}.bed
# true_peakfile=${filename}.INPUT.peaks.bed
# true_peak_directory=${outdirectory}
# bedtools sort -i ${true_peak_directory}/mm10.transcripts.${chr_title}.bed | bedtools merge  > ${true_peak_directory}/${true_peakfile}
true_peak_directory=${outdirectory}
true_peakfile=${filename}.INPUT.peaks.bed
true_comple_peak_directory=${outdirectory}
true_comple_peakfile=${filename}.INPUT.COMPLE.peaks.bed
## Generate true peak using gene transcripts area
python ${script_directory}/INPUT_ATAC_demoInputPeak.py ${true_peak_directory} mm9.transcripts.${chr_title}.bed unsorted.${true_peakfile}
bedtools sort -i ${true_peak_directory}/unsorted.${true_peakfile} | bedtools merge  > ${true_peak_directory}/${true_peakfile}
python ${script_directory}/INPUT_ATAC_ComplePeakFunction.py ${chr_title} ${true_peak_directory} ${true_peakfile} ${true_comple_peak_directory} ${true_comple_peakfile}
# cp ~/Projects/scATAC_Simulator/results/20221121_${filename}_INPUT_withCluster_mm9TSS/${true_peakfile} ./
# cp ~/Projects/scATAC_Simulator/results/20221121_${filename}_INPUT_withCluster_mm9TSS/${true_comple_peakfile} ./

######################## Ground Truth Peature Set ######################## 
## MACS3 peak calling results
# Store as referenced-peak.bed
# Take complementary regions and stored as referenced-non-peak.bed 
# MACS3 call peaks
# macs3 callpeak -f BAMPE -t $INPUT_bamfile -g mm -n ${outdirectory}/${filename}.MACS3 -B -q 0.01
ref_peak_directory=${outdirectory}
ref_peakfile=${filename}.MACS3.bed_peaks.bed
ref_comple_peak_directory=${outdirectory}
ref_comple_peakfile=${filename}.MACS3_thr10_peaks.COMPLE.bed

# Generate reference peak set
# reference peak and nonpeaks are copied from noninput results
bedtools sort -i ${directory}/${filename}.MACS3.bed_peaks.narrowPeak | bedtools merge  > ${ref_peak_directory}/${ref_peakfile}
# # Generate reference non peak set
python ${script_directory}/INPUT_ATAC_ComplePeakFunction.py ${chr_title} ${ref_peak_directory} ${ref_peakfile} ${ref_comple_peak_directory} ${ref_comple_peakfile}
# cp ~/Projects/scATAC_Simulator/results/20221121_${filename}_INPUT_withCluster_mm9TSS/${ref_peakfile} ./
# cp ~/Projects/scATAC_Simulator/results/20221130_${filename}_NONINPUT/${ref_comple_peakfile} ./

######################## Generate Count matrix ######################## 
## Construct assignment between true-feature set and referenced-feature set
# Assignment between true peak set and reference peak set
ref_marginal_file=BoneMarrow_62016_chr1.ref.peak.marginalcount.txt
ref_marginal_comple_file=BoneMarrow_62016_chr1.ref.complepeak.marginalcount.txt
assignment_file=${filename}.assigned.peaks.txt
assignment_comple_file=${filename}.COMPLE.assigned.peaks.txt
# Match ref and input peaks and nonpeaks
python ${script_directory}/INPUT_ATAC_MatchPeakFunction_20220302.py ${INPUT_bamfile} $ref_marginal_file $ref_marginal_comple_file $true_peak_directory $true_peakfile $true_comple_peakfile $ref_peak_directory $ref_peakfile $ref_comple_peakfile $outdirectory $assignment_file $assignment_comple_file
# # Peak (Run after non-peaks as nonpeak step produces both peak and nonpeak assignment, copy the original one to overwrite the newly produced)
# cp ~/Projects/scATAC_Simulator/results/20220407_BoneMarrow_62016_chr1_INPUT_withCluster_mm9TSS/${assignment_file} ./

## Construct assigned count matrix
INPUT_cells_barcode_file=/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas/cell_metadata_bonemarrow.csv
count_mat_filename=${filename}.assigned.countmatrix
count_mat_format=txt
#Extract referenced peaks set and convert to count matrices
# cp ~/Projects/scATAC_Simulator/results/20221121_BoneMarrow_62016_chr1_INPUT_withCluster_mm9TSS/${count_mat_filename}.txt ./ 
time(python ${script_directory}/INPUT_ATAC_BAM2CountMatrix.py ${INPUT_cells_barcode_file} ${outdirectory} ${assignment_file} ${INPUT_bamfile} ${outdirectory} ${count_mat_filename}.${count_mat_format}) # 17mins
#Extract referenced non-peaks set and convert to count matrices
count_mat_comple_filename=${filename}.assigned.COMPLE.countmatrix
time(python ${script_directory}/INPUT_ATAC_BAM2CountMatrix.py ${INPUT_cells_barcode_file} ${outdirectory} ${assignment_comple_file} ${INPUT_bamfile} ${outdirectory} ${count_mat_comple_filename}.${count_mat_format}) # 50mins
# cp ~/Projects/scATAC_Simulator/results/20220407_BoneMarrow_62016_chr1_INPUT_withCluster_mm9TSS/${count_mat_filename}.txt ~/Projects/scATAC_Simulator/results/20220407_BoneMarrow_62016_chr1_INPUT_withCluster_mm9TSS/${count_mat_comple_filename}.txt ./ 

## Train scDesign2 on reconstructed count matrix
# Synthetic matrix file named as filename.scDesign2Simulated.txt, note the first column with peak name is omitted 
synthetic_countmat_file=${count_mat_filename}.scDesign2Simulated.${count_mat_format}
synthetic_countmat_file_comple=${count_mat_comple_filename}.scDesign2Simulated.${count_mat_format}
# Peak
# cp ~/Projects/scATAC_Simulator/results/20221121_BoneMarrow_62016_chr1_INPUT_withCluster_mm9TSS/${synthetic_countmat_file} ./ 
time(Rscript ${script_directory}/INPUT_ATAC_SyntheticMat_withCluster_SelectCellType.R ${count_mat_filename} ${count_mat_format} ${outdirectory} ${outdirectory} ${INPUT_cells_barcode_file}) # 173m30.792s
# Non-Peak
time(Rscript ${script_directory}/INPUT_ATAC_SyntheticMat_withCluster_SelectCellType.R ${count_mat_comple_filename} ${count_mat_format} ${outdirectory} ${outdirectory} ${INPUT_cells_barcode_file}) # 86m8.349s


######################## Generate BAM file ######################## 
BED_filename_pre=${filename}.syntheticBAM
BED_COMPLE_filename_pre=${filename}.syntheticBAM.COMPLE
BED_filename_combined_pre=${filename}.syntheticBAM.combined
## Parsing bam files according to referenced features, modify the position according to true features
cd ${outdirectory}

# Peaks and Non-peaks, set 'random_noise_mode = False'
time(python3 ${script_directory}/INPUT_ATAC_GenerateBAMCoord_20221130.py # 6 min and 7 min

## Combine peak and comple.peak
# cat ${outdirectory}/${BED_filename_pre}.read1.bed ${outdirectory}/${BED_COMPLE_filename_pre}.read1.bed | sort -k1,1 -k2,2n > ${outdirectory}/${BED_filename_combined_pre}.read1.bed
# cat ${outdirectory}/${BED_filename_pre}.read2.bed ${outdirectory}/${BED_COMPLE_filename_pre}.read2.bed | sort -k1,1 -k2,2n > ${outdirectory}/${BED_filename_combined_pre}.read2.bed
cat ${outdirectory}/${BED_filename_pre}.read1.bed ${outdirectory}/${BED_COMPLE_filename_pre}.read1.bed  > ${outdirectory}/${BED_filename_combined_pre}.read1.unsort.bed
cat ${outdirectory}/${BED_filename_pre}.read2.bed ${outdirectory}/${BED_COMPLE_filename_pre}.read2.bed  > ${outdirectory}/${BED_filename_combined_pre}.read2.unsort.bed


## Convert bed files to FASTQ files
# Double check each fasta read header contains single cell barcode, read 1or2 index, single cell's read pseudo order
# Option: -s Forcing the extracted sequence to reflect the requested strand
# Option: -name Using the BED “name” column as a FASTA header.
# bedtools getfasta -fi ~/Projects/scATAC_Simulator/data/mm9_genome/chr1.fa -bed  ~/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas/BoneMarrow_62016_chr1.bed -fo ~/Projects/scATAC_Simulator/data/bedTobam_output/BoneMarrow_62016_chr1.bed2fa.fa  -s -name
referenceGenome_name=NCBIM37.genome
referenceGenome_dir=/home/gayan/Projects/scATAC_Simulator/data/mm9_ref_genome_GECODE
referenceGenome_file=${referenceGenome_dir}/${referenceGenome_name}.fa
output_BAM_pre=${filename}.syntheticBAM

time(bedtools getfasta  -s -nameOnly -fi ${referenceGenome_file} -bed ${outdirectory}/${BED_filename_combined_pre}.read1.unsort.bed -fo ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.strand.fa) # 4 mins
time(bedtools getfasta  -s -nameOnly -fi ${referenceGenome_file} -bed ${outdirectory}/${BED_filename_combined_pre}.read2.unsort.bed -fo ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.strand.fa) # 4 mins
sed '/^>/s/.\{3\}$//' ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.strand.fa > ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.fa # remove (-)
sed '/^>/s/.\{3\}$//' ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.strand.fa > ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fa # remove (-)

# time(${bedtools_directory}/bedtools getfasta -fi ${referenceGenome_file} -bed ${outdirectory}/${BED_filename_combined_pre}.read1.bed -fo ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.fa -nameOnly) # 4 mins
# awk 'OFS="\t" {print $0, "-"}'  ${outdirectory}/${BED_filename_combined_pre}.read2.bed > ${outdirectory}/${BED_filename_combined_pre}.read2.strand.bed
# time(bedtools getfasta -s -fi ${referenceGenome_file} -bed ${outdirectory}/${BED_filename_combined_pre}.read2.strand.bed -fo ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.strand.fa -nameOnly)
# sed '/^>/s/.\{3\}$//' ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.strand.fa > ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fa # remove (-)
# time(${bedtools_directory}/bedtools getfasta -fi ${referenceGenome_file} -bed ${outdirectory}/${BED_filename_combined_pre}.read2.bed -fo ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fa -nameOnly)

# Convert fasta to fastq
# ~/Tools/seqtk/seqtk/seqtk seq -F '#' ~/Projects/scATAC_Simulator/data/bedTobam_output/BoneMarrow_62016_chr1.bed2fa.fa > ~/Projects/scATAC_Simulator/data/bedTobam_output/BoneMarrow_62016_chr1.bed2fa.fq
time(~/Tools/seqtk/seqtk/seqtk seq -F 'F' ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.fa > ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.fq)
time(~/Tools/seqtk/seqtk/seqtk seq -F 'F' ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fa > ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fq)

# Sort fastq files
time(cat ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq)
time(cat ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq)

# Bowtie align sam
bowtie2_directory=~/Tools/bowtie2
# cd ${referenceGenome_dir}
# bowtie2-build ${referenceGenome_file} ${referenceGenome_name}
# time(bowtie2 -x ${referenceGenome_dir}/${referenceGenome_name} -1 ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq -2 ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq | samtools view -bS - > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.quickversion.bam) # 151m37.660s
# samtools sort ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.quickversion.bam > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.quickversion.sorted.bam
# samtools index ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.quickversion.sorted.bam
time(bowtie2 --minins 0 --maxins 1200 -x ${referenceGenome_dir}/${referenceGenome_name} -1 ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq -2 ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq | samtools view -bS - > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.bam) # 99 minutes
# add cell barcode to bam file
time(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.bam -H > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.header.sam)
time(cat <( cat ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.header.sam ) \
<( paste <(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.bam ) <(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.bam | cut -f1 | cut -d':' -f1 | sed -e 's/^/CB:Z:/')) | ${samtools_directory}/samtools view -bS - > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.bam)
rm ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.header.sam 

# Sort bam
samtools sort ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.bam > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam
samtools index ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam

 
######################## 20221121 Peak Calling analysis  ######################## 
outdirectory=/home/gayan/Projects/scATAC_Simulator/results/20221130_BoneMarrow_62016_chr1_INPUT_withCluster_mm9TSS
output_BAM_pre=BoneMarrow_62016_chr1.syntheticBAM
chr_title=chr1
filename=BoneMarrow_62016_${chr_title}
script_directory=/home/gayan/Projects/scATAC_Simulator/scripts/03092022/sciATAC/INPUT_ATAC-seq/WithCB_scDesign2forsynthetic
samtools_directory=/home/gayan/Tools/samtools/bin
macs3_directory=/home/gayan/.local/bin
bedtools_directory=/home/gayan/Tools/bedtools/bedtools2/bin
export PATH=${macs3_directory}:${samtools_directory}:${bedtools_directory}:${PATH}

## Calculate marginal count for RPKM
ref_gt_peak=${filename}.MACS3.bed_peaks.bed
ref_gt_nonpeak=${filename}.MACS3_thr10_peaks.COMPLE.bed
true_peakfile=${filename}.INPUT.peaks.bed
true_comple_peakfile=${filename}.INPUT.COMPLE.peaks.bed
scReadSim_bamfile=${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam
INPUT_bamfile=${directory}/${filename}.bam  
count_mat_filename=${filename}.countmatrix
count_mat_comple_filename=${filename}.COMPLE.countmatrix
# Synthetic BAM read coverage over groundtruth peaks and non-peaks
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} ${true_peakfile} ${scReadSim_bamfile} ${count_mat_filename}.scDesign2Simulated.GroundTruth.MarginalFeatureCount.txt
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} ${true_comple_peakfile} ${scReadSim_bamfile} ${count_mat_comple_filename}.scDesign2Simulated.GroundTruth.MarginalFeatureCount.txt
# Real BAM read coverage over real peaks and non-peaks
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} ${ref_gt_peak} ${INPUT_bamfile} ${count_mat_filename}.MarginalFeatureCount.txt
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} ${ref_gt_nonpeak} ${INPUT_bamfile} ${count_mat_comple_filename}.MarginalFeatureCount.txt


## MACS3 call peak 
synthetic_MACS3_peakcall_filename=${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.MACS3
macs3 callpeak -f BAMPE -t ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam -g mm -n ${synthetic_MACS3_peakcall_filename} -B -q 0.05

## SEARC call peak
cp  ~/Projects/scATAC_Simulator/results/20221121_${filename}_INPUT_withCluster_mm9TSS/mm9.chrom.sizes ./
synthetic_searc_peakcall_filename=${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.searc
scReadSim_bam=${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam
time(samtools sort -n -o ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.namesorted.bam ${scReadSim_bam} )
time(bedtools bamtobed -bedpe -i ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.namesorted.bam > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.bed)
awk '$1==$4 && $6-$2 < 1000 {print $0}' ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.bed > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.clean.bed
cut -f 1,2,6  ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.clean.bed | sort -k1,1 -k2,2n -k3,3n  | grep chr1 >  ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.fragments.bed
bedtools genomecov -bg -i ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.fragments.bed -g ${outdirectory}/mm9.chrom.sizes > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.fragments.bedgraph
searc_path=/home/gayan/Tools/searc/SEACR
export PATH=$searc_path:$PATH
cd $searc_path
bash SEACR_1.3.sh ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.fragments.bedgraph 0.5 norm stringent ${outdirectory}/${synthetic_searc_peakcall_filename}

## SEARC call peak 
# Test with fewer discoveries 20230227
synthetic_searc_peakcall_filename=${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.searc.test20230227
searc_path=/home/gayan/Tools/searc/SEACR
export PATH=$searc_path:$PATH
cd $searc_path
bash SEACR_1.3.sh ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.fragments.bedgraph 0.01 norm stringent ${synthetic_searc_peakcall_filename}


## HOMER call peak
scReadSim_bam=${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam
PATH=$PATH:/home/gayan/Tools/HOMER/bin/
mkdir ${outdirectory}/HOMER_tags
makeTagDirectory ${outdirectory}/HOMER_tags/ ${scReadSim_bam}
findPeaks ${outdirectory}/HOMER_tags/ -region -o ${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.HOMER -gsize 2.5e9 minDist 150
pos2bed.pl ${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.HOMER | bedtools sort | bedtools merge > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.HOMER.merged.bed

## HMMRATAC call peak
scReadSim_bam=${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam
samtools view -H ${scReadSim_bam}| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > ${outdirectory}/HMMRATAC.genome.info
export PATH=/home/gayan/Tools/HMMRATAC:$PATH
java -jar /home/gayan/Tools/HMMRATAC/HMMRATAC_V1.2.10_exe.jar -b ${scReadSim_bam} -i ${scReadSim_bam}.bai -g ${outdirectory}/HMMRATAC.genome.info -o ${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.HMMRATAC.bed
mv BoneMarrow_62016_chr1.syntheticBAM.synthetic.chr1.sorted.bam.HMMRATAC.bed_peaks.gappedPeak ${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.HMMRATAC.bed

######################## 20220421 Introduce Error using NONINPUT_ATAC_ReadError.py  ######################## 
### Plot 1: fragment length
mkdir -p ${outdirectory}/VerifyCoverage/deepTools_scReadSim
scReadSim_bamfile=${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam
# bamPEFragmentSize \
# -hist ${outdirectory}/VerifyCoverage/deepTools_scReadSim/fragmentSize.png \
# -T "Fragment length" \
# -b $scReadSim_bamfile \
# $INPUT_bamfile \
# --maxFragmentLength 1000 \
# --outRawFragmentLengths ${outdirectory}/VerifyCoverage/deepTools_scReadSim/fragmentSize.txt \
# --samplesLabel scReadSim Real

bamPEFragmentSize \
-hist ${outdirectory}/VerifyCoverage/deepTools_scReadSim/fragmentSize_scReadSim.png \
-T "Fragment size" \
-b $scReadSim_bamfile \
--maxFragmentLength 1000 \
--outRawFragmentLengths ${outdirectory}/VerifyCoverage/deepTools_scReadSim/fragmentSize_scReadSim.txt \
--samplesLabel scReadSim

bamPEFragmentSize \
-hist ${outdirectory}/VerifyCoverage/deepTools_scReadSim/fragmentSizer_Real.png \
-T "Fragment size" \
-b $INPUT_bamfile \
--maxFragmentLength 1000 \
--outRawFragmentLengths ${outdirectory}/VerifyCoverage/deepTools_scReadSim/fragmentSize_Real.txt \
--samplesLabel Real


### Plot 2: Kmer spectrum
## Convert real bam to fasta
samtools fasta -n $INPUT_bamfile -1 ${outdirectory}/real_bam2fasta_read1.fa -2 ${outdirectory}/real_bam2fasta_read2.fa
## Kmer on real fa
mkdir -p ${outdirectory}/VerifyCoverage/kmer_real
mkdir -p ${outdirectory}/VerifyCoverage/kmer_scReadSim
for kmer_len in 11 21 31 41 51
do 
  jellyfish count -m ${kmer_len} -o ${outdirectory}/VerifyCoverage/kmer_real/mer${kmer_len}_counts_read1.jf -s 100M -t 10 -C ${outdirectory}/real_bam2fasta_read1.fa
  jellyfish count -m ${kmer_len} -o ${outdirectory}/VerifyCoverage/kmer_real/mer${kmer_len}_counts_read2.jf -s 100M -t 10 -C ${outdirectory}/real_bam2fasta_read1.fa
  jellyfish histo ${outdirectory}/VerifyCoverage/kmer_real/mer${kmer_len}_counts_read1.jf > ${outdirectory}/VerifyCoverage/kmer_real/mer${kmer_len}_counts_read1_hist.txt
  jellyfish histo ${outdirectory}/VerifyCoverage/kmer_real/mer${kmer_len}_counts_read2.jf > ${outdirectory}/VerifyCoverage/kmer_real/mer${kmer_len}_counts_read2_hist.txt
  ## Kmer on simulated fa
  jellyfish count -m ${kmer_len} -o ${outdirectory}/VerifyCoverage/kmer_scReadSim/mer${kmer_len}_counts_read1.jf -s 100M -t 10 -C ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.fa
  jellyfish count -m ${kmer_len} -o ${outdirectory}/VerifyCoverage/kmer_scReadSim/mer${kmer_len}_counts_read2.jf -s 100M -t 10 -C ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fa
  jellyfish histo ${outdirectory}/VerifyCoverage/kmer_scReadSim/mer${kmer_len}_counts_read1.jf > ${outdirectory}/VerifyCoverage/kmer_scReadSim/mer${kmer_len}_counts_read1_hist.txt
  jellyfish histo ${outdirectory}/VerifyCoverage/kmer_scReadSim/mer${kmer_len}_counts_read2.jf > ${outdirectory}/VerifyCoverage/kmer_scReadSim/mer${kmer_len}_counts_read2_hist.txt
done


### Plot 3: error rates
## fgbio 
mkdir -p ${outdirectory}/VerifyQuality/fgbio
# java -jar /home/gayan/Tools/picard/build/libs/picard.jar CreateSequenceDictionary \
#       -R $referenceGenome_file \
#       -O ${referenceGenome_file}.dict
java -jar /home/gayan/Tools/fgbio/target/scala-2.13/fgbio-2.0.1-e884860-SNAPSHOT.jar ErrorRateByReadPosition -i $INPUT_bamfile -r $referenceGenome_file -o ${outdirectory}/VerifyQuality/fgbio/Real --collapse false

## Introduce error rate bam
python3 ${script_directory}/NOINPUT_ATAC_ReadError_20220421.py ${outdirectory}/VerifyQuality/fgbio/Real.error_rate_by_read_position.txt ${outdirectory} ${BED_filename_combined_pre}
time(cat ${outdirectory}/${BED_filename_combined_pre}.ErrorIncluded.read1.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > ${outdirectory}/${BED_filename_combined_pre}.ErrorIncluded.read1.bed2fa.sorted.fq)
time(cat ${outdirectory}/${BED_filename_combined_pre}.ErrorIncluded.read2.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > ${outdirectory}/${BED_filename_combined_pre}.ErrorIncluded.read2.bed2fa.sorted.fq)
## Bowtie align sam
# cd ${referenceGenome_dir}
bowtie2_directory=~/Tools/bowtie2
# bowtie2-build ${referenceGenome_file} ${referenceGenome_name}
time(bowtie2 --minins 0 --maxins 1200 -x ${referenceGenome_dir}/${referenceGenome_name} -1 ${outdirectory}/${BED_filename_combined_pre}.ErrorIncluded.read1.bed2fa.sorted.fq -2 ${outdirectory}/${BED_filename_combined_pre}.ErrorIncluded.read2.bed2fa.sorted.fq -S ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.sam) # 151m37.660s
${samtools_directory}/samtools view -S -b ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.sam > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.bam
## add cell barcode to bam file
time(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.bam -H > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.header.sam)
time(cat <( cat ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.header.sam ) \
<( paste <(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.bam ) <(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.bam | cut -f1 | cut -d':' -f1 | sed -e 's/^/CB:Z:/')) | ${samtools_directory}/samtools view -bS - > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.ErrorIncluded.bam)
rm ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.header.sam 
## Sort and index bam 
${samtools_directory}/samtools sort ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.ErrorIncluded.bam > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.ErrorIncluded.sorted.bam
${samtools_directory}/samtools index ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.ErrorIncluded.sorted.bam

## fgbio 
java -jar /home/gayan/Tools/fgbio/target/scala-2.13/fgbio-2.0.1-e884860-SNAPSHOT.jar ErrorRateByReadPosition -i $scReadSim_bamfile -r $referenceGenome_fissle -o ${outdirectory}/VerifyQuality/fgbio/scReadSim --collapse false
java -jar /home/gayan/Tools/fgbio/target/scala-2.13/fgbio-2.0.1-e884860-SNAPSHOT.jar ErrorRateByReadPosition -i ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.ErrorIncluded.sorted.bam -r $referenceGenome_file -o ${outdirectory}/VerifyQuality/fgbio/scReadSim_manualError --collapse false



































