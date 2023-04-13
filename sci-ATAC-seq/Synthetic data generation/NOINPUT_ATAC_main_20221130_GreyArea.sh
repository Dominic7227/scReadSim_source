#!/bin/bash

# Task: mouse sci-ATAC-seq synthetic reads simulation
# Date: 20221130
# Notes
# 1. No designed ground-truth peaks, use trustworthy peaks and non-peaks as ground-truth peaks and non-peaks


#####################################################################
############################ Prepare the data #######################
#####################################################################
full_directory=/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas
full_filename=BoneMarrow_62016
chr_title=chr1
filename=BoneMarrow_62016_${chr_title}
directory=/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas
INPUT_bamfile=${directory}/${filename}.bam  
# Tools directory
script_directory=/home/gayan/Projects/scATAC_Simulator/scripts/03092022/sciATAC/ATAC-seq
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
# outdirectory=~/Projects/scATAC_Simulator/results/20211126_${filename}_NONINPUT
outdirectory=/home/gayan/Projects/scATAC_Simulator/results/20221130_BoneMarrow_62016_chr1_NONINPUT

INPUT_cells_barcode_file=/home/gayan/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas/cell_metadata_bonemarrow.csv
OUTPUT_cells_barcode_file=${outdirectory}/synthetic_cell_barcode.txt
mkdir ${outdirectory}
cd ${outdirectory}

######################## MACS Call Peaks ######################## 
## MACS3 peak calling results
# Store as referenced-peak.bed
# Take complementary regions and stored as referenced-non-peak.bed 
# MACS3 call peaks
MACS3_peakname_pre=${filename}.MACS3
# Default peak calling method for ground truth peak
time(${macs3_directory}/macs3 callpeak -f BAMPE -t ${INPUT_bamfile} -g mm -n ${outdirectory}/${MACS3_peakname_pre} -B -q 0.01 --outdir ${outdirectory}) #997m39.768s
# Less stringent rule for ground truth non-peak
time(${macs3_directory}/macs3 callpeak -f BAMPE -t ${INPUT_bamfile} -g mm -n ${outdirectory}/${MACS3_peakname_pre}_thr5 -B -q 0.05 --outdir ${outdirectory}) #997m39.768s
time(${macs3_directory}/macs3 callpeak -f BAMPE -t ${INPUT_bamfile} -g mm -n ${outdirectory}/${MACS3_peakname_pre}_thr10 -B -q 0.1 --outdir ${outdirectory}) #997m39.768s
time(${macs3_directory}/macs3 callpeak -f BAMPE -t ${INPUT_bamfile} -g mm -n ${outdirectory}/${MACS3_peakname_pre}_thr50 -B -q 0.5 --outdir ${outdirectory}) #997m39.768s
# Merge intersecting peaks for illustration in IGV browser to identify grey regions
bedtools sort -i ${outdirectory}/${MACS3_peakname_pre}_thr5_peaks.narrowPeak | ${bedtools_directory}/bedtools merge  > ${outdirectory}/${MACS3_peakname_pre}_thr5_peaks.bed
bedtools sort -i ${outdirectory}/${MACS3_peakname_pre}_thr10_peaks.narrowPeak | ${bedtools_directory}/bedtools merge  > ${outdirectory}/${MACS3_peakname_pre}_thr10_peaks.bed
bedtools sort -i ${outdirectory}/${MACS3_peakname_pre}_thr50_peaks.narrowPeak | ${bedtools_directory}/bedtools merge  > ${outdirectory}/${MACS3_peakname_pre}_thr50_peaks.bed

# Generate reference peak set
ref_gt_peak=${MACS3_peakname_pre}_peaks.bed # P1
ref_gt_nonpeak=${MACS3_peakname_pre}_thr10_peaks.COMPLE.bed # P2
ref_gt_greyarea=${MACS3_peakname_pre}_GreyAreas.bed
${bedtools_directory}/bedtools sort -i ${outdirectory}/${MACS3_peakname_pre}_peaks.narrowPeak | ${bedtools_directory}/bedtools merge  > ${ref_peak_directory}/${ref_gt_peak}
# cp 20221121_BoneMarrow_62016_chr1_NONINPUT/${ref_gt_peak} ${outdirectory}
# Grey regoins (P1 U P2C)Cs
# P2C
python3 ${script_directory}/NOINPUT_ATAC_ComplePeakFunction.py ${chr_title} ${outdirectory} ${MACS3_peakname_pre}_thr10_peaks.bed ${outdirectory} ${MACS3_peakname_pre}_thr10_peaks.COMPLE.bed
# P1 U P2C
cat ${MACS3_peakname_pre}_thr10_peaks.COMPLE.bed ${MACS3_peakname_pre}_peaks.bed | bedtools sort | bedtools merge > ${outdirectory}/${MACS3_peakname_pre}_P1UP2C.bed
# Grey areas: (P1 U P2C)C
python3 ${script_directory}/NOINPUT_ATAC_ComplePeakFunction.py ${chr_title} ${outdirectory} ${MACS3_peakname_pre}_P1UP2C.bed ${outdirectory} ${MACS3_peakname_pre}_GreyAreas.bed


######################## Generate Count matrix ######################## 
count_mat_filename=${filename}.countmatrix
count_mat_comple_filename=${filename}.COMPLE.countmatrix
count_mat_format=txt
# Ground truth peak
time(python3 ${script_directory}/NOINPUT_ATAC_BAM2CountMatrix.py ${INPUT_cells_barcode_file} ${outdirectory} ${ref_gt_peak} ${INPUT_bamfile} ${outdirectory} ${count_mat_filename}.${count_mat_format}) 
# cp 20221121_BoneMarrow_62016_chr1_NONINPUT/${count_mat_filename}.txt ${outdirectory}/
# Ground truth non-peak
time(python3 ${script_directory}/NOINPUT_ATAC_BAM2CountMatrix.py ${INPUT_cells_barcode_file} ${outdirectory} ${ref_gt_nonpeak} ${INPUT_bamfile} ${outdirectory} ${count_mat_comple_filename}.${count_mat_format}) # 127m12.784s

## Train scDesign2 on reconstructed count matrix
synthetic_countmat_file=${count_mat_filename}.scDesign2Simulated.${count_mat_format}
synthetic_countmat_file_comple=${count_mat_comple_filename}.scDesign2Simulated.${count_mat_format}
# Ground truth peak
time(Rscript ${script_directory}/NOINPUT_ATAC_SyntheticMat.R ${count_mat_filename} ${count_mat_format} ${outdirectory} ${outdirectory} $INPUT_cells_barcode_file) # 86m8.349s
# cp 20221121_BoneMarrow_62016_chr1_NONINPUT/${synthetic_countmat_file} ${outdirectory}/
# Ground truth non-peak
time(Rscript ${script_directory}/NOINPUT_ATAC_SyntheticMat.R ${count_mat_comple_filename} ${count_mat_format} ${outdirectory} ${outdirectory} $INPUT_cells_barcode_file) # 86m8.349s


######################## Generate BAM file ######################## 
BED_filename_pre=${filename}.syntheticBAM.CBincluded
BED_COMPLE_filename_pre=${filename}.syntheticBAM.COMPLE.CBincluded
BED_GreyArea_filename_pre=${filename}.syntheticBAM.COMPLE.CBincluded.GreyArea
BED_filename_combined_pre=${filename}.syntheticBAM.combined.CBincluded
# BED_filename_pre=${filename}.syntheticBAM
# BED_COMPLE_filename_pre=${filename}.syntheticBAM.COMPLE
# BED_filename_combined_pre=${filename}.syntheticBAM.combined

## Parsing bam files according to referenced features, modify the position according to true features
cd ${outdirectory}

# Peaks, set 'random_noise_mode = False'
# Non-peaks, set 'random_noise_mode = False'
time(python3 ${script_directory}/NOINPUT_ATAC_GenerateBAMCoord_20221130.py )

## Combine peak and comple.peak
# cat ${outdirectory}/${BED_filename_pre}.read1.bed ${outdirectory}/${BED_COMPLE_filename_pre}.read1.bed | sort -k1,1 -k2,2n | cut -f1-5 > ${outdirectory}/${BED_filename_combined_pre}.read1.bed
# cat ${outdirectory}/${BED_filename_pre}.read2.bed ${outdirectory}/${BED_COMPLE_filename_pre}.read2.bed | sort -k1,1 -k2,2n | cut -f1-5 > ${outdirectory}/${BED_filename_combined_pre}.read2.bed
cat ${outdirectory}/${BED_filename_pre}.read1.bed ${outdirectory}/${BED_COMPLE_filename_pre}.read1.bed ${outdirectory}/${BED_GreyArea_filename_pre}.read1.bed > ${outdirectory}/${BED_filename_combined_pre}.read1.unsort.bed
cat ${outdirectory}/${BED_filename_pre}.read2.bed ${outdirectory}/${BED_COMPLE_filename_pre}.read2.bed ${outdirectory}/${BED_GreyArea_filename_pre}.read2.bed > ${outdirectory}/${BED_filename_combined_pre}.read2.unsort.bed

## Convert bed files to FASTQ files
# Double check each fasta read header contains single cell barcode, read 1or2 index, single cell's read pseudo order
# Option: -s Forcing the extracted sequence to reflect the requested strand
# Option: -name Using the BED “name” column as a FASTA header.
# bedtools getfasta -fi ~/Projects/scATAC_Simulator/data/mm9_genome/chr1.fa -bed  ~/Projects/scATAC_Simulator/data/sciATAC_MouseAtlas/BoneMarrow_62016_chr1.bed -fo ~/Projects/scATAC_Simulator/data/bedTobam_output/BoneMarrow_62016_chr1.bed2fa.fa  -s -name
# referenceGenome_name=GRCm38.primary_assembly.genome
# referenceGenome_dir=/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE
referenceGenome_name=NCBIM37.genome
referenceGenome_dir=/home/gayan/Projects/scATAC_Simulator/data/mm9_ref_genome_GECODE
referenceGenome_file=${referenceGenome_dir}/${referenceGenome_name}.fa
output_BAM_pre=${BED_filename_pre}
time(bedtools getfasta  -s -nameOnly -fi ${referenceGenome_file} -bed ${outdirectory}/${BED_filename_combined_pre}.read1.unsort.bed -fo ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.strand.fa) # 4 mins
time(bedtools getfasta  -s -nameOnly -fi ${referenceGenome_file} -bed ${outdirectory}/${BED_filename_combined_pre}.read2.unsort.bed -fo ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.strand.fa) # 4 mins
sed '/^>/s/.\{3\}$//' ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.strand.fa > ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.fa # remove (-)
sed '/^>/s/.\{3\}$//' ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.strand.fa > ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fa # remove (-)

# Convert fasta to fastq
# ~/Tools/seqtk/seqtk/seqtk seq -F '#' ~/Projects/scATAC_Simulator/data/bedTobam_output/BoneMarrow_62016_chr1.bed2fa.fa > ~/Projects/scATAC_Simulator/data/bedTobam_output/BoneMarrow_62016_chr1.bed2fa.fq
time(~/Tools/seqtk/seqtk/seqtk seq -F 'F' ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.fa > ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.fq)
time(~/Tools/seqtk/seqtk/seqtk seq -F 'F' ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fa > ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fq)

time(cat ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq)
time(cat ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq)

# Bowtie align sam
# cd ${referenceGenome_dir}
bowtie2_directory=~/Tools/bowtie2
# bowtie2-build ${referenceGenome_file} ${referenceGenome_name}
# time(bowtie2 -x ${referenceGenome_dir}/${referenceGenome_name} -1 ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq -2 ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq | samtools view -bS - > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.quickversion.bam) # 151m37.660s
# Full version mapping
time(bowtie2 --minins 0 --maxins 1200 -x ${referenceGenome_dir}/${referenceGenome_name} -1 ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq -2 ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq | samtools view -bS - > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.bam) # 400 minutes
# add cell barcode to bam file
time(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.bam -H > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.header.sam)
time(cat <( cat ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.header.sam ) \
<( paste <(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.bam ) <(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.bam | cut -f1 | cut -d':' -f1 | sed -e 's/^/CB:Z:/')) | ${samtools_directory}/samtools view -bS - > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.bam)
rm ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.header.sam 

# Sort and index bam 
${samtools_directory}/samtools sort ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.bam > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam
${samtools_directory}/samtools index ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam

###########################################################################
######################## 20221129 Producing Results  ######################
###########################################################################
### Plot 0: Peak Calling
## MACS3 call peak 
synthetic_MACS3_peakcall_filename=${output_BAM_pre}.synthetic.${chr_title}.sorted.bam.MACS3
${macs3_directory}/macs3 callpeak -f BAMPE -t ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam -g mm -n ${synthetic_MACS3_peakcall_filename} -B -q 0.05
${macs3_directory}/macs3 callpeak -f BAMPE -t ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam -g mm -n ${synthetic_MACS3_peakcall_filename}.stringent01 -B -q 0.01
bedtools sort -i ${outdirectory}/${synthetic_MACS3_peakcall_filename}_peaks.narrowPeak | ${bedtools_directory}/bedtools merge  > ${outdirectory}/${synthetic_MACS3_peakcall_filename}_peaks.narrowPeak.bed 
python3 ${script_directory}/NOINPUT_ATAC_ComplePeakFunction.py ${chr_title} ${outdirectory} ${synthetic_MACS3_peakcall_filename}_peaks.narrowPeak.bed ${outdirectory} ${synthetic_MACS3_peakcall_filename}_peaks.narrowPeak.COMPLE.bed

# More peaks
${macs3_directory}/macs3 callpeak -f BAMPE -t ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam -g mm -n ${synthetic_MACS3_peakcall_filename}.lessStringent -B -q 0.5


## Calculate marginal count for RPKM
ref_gt_peak=${MACS3_peakname_pre}_peaks.bed
ref_gt_nonpeak=${MACS3_peakname_pre}_thr10_peaks.COMPLE.bed
scReadSim_bamfile=${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam
INPUT_bamfile=${directory}/${filename}.bam  
count_mat_filename=${filename}.countmatrix
count_mat_comple_filename=${filename}.COMPLE.countmatrix
# Synthetic BAM read coverage over real peaks and non-peaks
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} ${ref_gt_peak} ${scReadSim_bamfile} ${count_mat_filename}.scDesign2Simulated.MarginalFeatureCount.txt
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} ${ref_gt_nonpeak} ${scReadSim_bamfile} ${count_mat_comple_filename}.scDesign2Simulated.MarginalFeatureCount.txt
# Real BAM read coverage over real peaks and non-peaks
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} ${ref_gt_peak} ${INPUT_bamfile} ${count_mat_filename}.MarginalFeatureCount.txt
python3 ${script_directory}/NOINPUT_ATAC_VerifyRead_20221130.py ${outdirectory} ${ref_gt_nonpeak} ${INPUT_bamfile} ${count_mat_comple_filename}.MarginalFeatureCount.txt


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

# ### fastqc for quality 
# mkdir ${outdirectory}/VerifyQuality/fastqc_scReadSim_manualError
# ${fastqc_dir}/fastqc ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.ErrorIncluded.sorted.bam -o ${outdirectory}/VerifyQuality/fastqc_scReadSim_manualError -f bam






