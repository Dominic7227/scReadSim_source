#!/bin/bash

# Task: mouse 10X scRNA-seq
# Date: 20221011


#####################################################################
############################ Prepare the data #######################
#####################################################################
full_filename=e18_mouse_brain_fresh_5k_gex_possorted_bam
full_directory=/home/gayan/Projects/scATAC_Simulator/data/10X_MOUSE_BRAIN_ATACandRNA
directory=${full_directory}/split.${full_filename}


mkdir ${directory}

# Tools directory
script_directory=/home/gayan/Projects/scATAC_Simulator/scripts/03092022/10X_MOUSE_BRAIN/RNA-seq
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
| ${samtools_directory}/samtools view -bS - > ${full_directory}/${full_filename}.scReadSim.bam) # 241m40.507s
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

# User input parameter
# 1. BAM file (attached barcodes)
chr_title=chr1
filename=${full_filename}_${chr_title}
outdirectory=~/Projects/scATAC_Simulator/results/20221011_${filename}_NONINPUT

INPUT_bamfile=${directory}/${filename}.bam  
INPUT_cells_barcode_file=/home/gayan/Projects/scATAC_Simulator/data/10X_MOUSE_BRAIN_ATACandRNA/filtered_feature_bc_matrix/barcodes.tsv
OUTPUT_cells_barcode_file=${outdirectory}/synthetic_cell_barcode.txt
mkdir ${outdirectory}
cd ${outdirectory}

######################## Extract barcodes from BAM file ######################## 
# time(${samtools_directory}/samtools view ${full_directory}/${full_filename}.bam | cut -f 12- | tr "\t" "\n"  | grep  "^XC:Z:"  | cut -d ':' -f 3 | sort | uniq > ${directory}/${full_filename}.barcodes.txt)

######################## Reference feature set ######################## 
INPUT_ref_peak_directory=/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE/gene_region
feature_set_file=gencode.vM10.annotation.gene_region.${chr_title}.merged.bed
ref_peakfile=${filename}_peaks.bed
ref_comple_peakfile=${filename}_peaks.COMPLE.bed

# ${bedtools_directory}/bedtools sort -i ${INPUT_ref_peak_directory}/${feature_set_file}| ${bedtools_directory}/bedtools merge  > ${outdirectory}/${ref_peakfile}
python3 ${script_directory}/NOINPUT_RNA_ComplePeakFunction_20220603.py $INPUT_bamfile $samtools_directory $bedtools_directory $outdirectory $genome_annotation $genome_size_file $ref_peakfile $ref_comple_peakfile
# Copy from 20220608
# cp /home/gayan/Projects/scATAC_Simulator/results/20220608_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/${ref_peakfile} ${outdirectory}
# cp /home/gayan/Projects/scATAC_Simulator/results/20220608_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/${ref_comple_peakfile} ${outdirectory}

######################## Generate Count matrices ######################## 
## Construct read count matrix
#Extract referenced peaks set and convert to count matrices
UMI_count_mat_filename=${filename}.UMI.countmatrix
UMI_count_mat_comple_filename=${filename}.UMI.COMPLE.countmatrix
# count_mat_filename=${filename}.countmatrix
# count_mat_comple_filename=${filename}.COMPLE.countmatrix
count_mat_format=txt
# Generate UMI count for foreground features
time(python3 ${script_directory}/NOINPUT_RNA_BAM2CountMatrix_20220928.py ${INPUT_cells_barcode_file} ${outdirectory}/${ref_peakfile} ${INPUT_bamfile} ${outdirectory} ${UMI_count_mat_filename})
# Copy from 20220608
# cp /home/gayan/Projects/scATAC_Simulator/results/20220926_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/${UMI_count_mat_filename}.txt ${outdirectory}
# Generate UMI count for background features
time(python3 ${script_directory}/NOINPUT_RNA_BAM2CountMatrix_20220928.py ${INPUT_cells_barcode_file} ${outdirectory}/${ref_comple_peakfile} ${INPUT_bamfile} ${outdirectory} ${UMI_count_mat_comple_filename}) # 30m4.198s
# cp /home/gayan/Projects/scATAC_Simulator/results/20220926_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/${UMI_count_mat_comple_filename}.txt ${outdirectory}

## Train scDesign2 on real count matrix and UMI count matrix
# Gene regions
time(Rscript ${script_directory}/NOINPUT_RNA_SyntheticMat.R ${UMI_count_mat_filename} ${count_mat_format} ${outdirectory} ${outdirectory} 1) 
# InterGenetic regions
time(Rscript ${script_directory}/NOINPUT_RNA_SyntheticMat.R ${UMI_count_mat_comple_filename} ${count_mat_format} ${outdirectory} ${outdirectory} 0) # 59min for 2165 peaks


######################## Generate BED file ######################## 
# cellnumberfile=${outdirectory}/${count_mat_filename}.scDesign2Simulated.nReadRegionmargional.txt
# cellnumberfile_comple=${outdirectory}/${count_mat_comple_filename}.scDesign2Simulated.nReadRegionmargional.txt
# synthetic_countmat_file=${count_mat_filename}.scDesign2Simulated.${count_mat_format}
# synthetic_countmat_file_comple=${count_mat_comple_filename}.scDesign2Simulated.${count_mat_format}
synthetic_UMI_count_mat_file=${UMI_count_mat_filename}.scDesign2Simulated.${count_mat_format}
synthetic_UMI_count_mat_file_comple=${UMI_count_mat_comple_filename}.scDesign2Simulated.${count_mat_format}
coordinate_file=BAMfile_coordinates.txt
coordinate_COMPLE_file=BAMfile_halfsampled_COMPLE_coordinates.txt
# BED_filename_pre=${filename}.syntheticBAM
# BED_COMPLE_filename_pre=${filename}.syntheticBAM.COMPLE
# BED_filename_combined_pre=${filename}.syntheticBAM.combined
BED_filename_pre=${filename}.syntheticBAM.20221011_UMITransLevel2
BED_COMPLE_filename_pre=${filename}.syntheticBAM.COMPLE.20221011_UMITransLevel2
BED_filename_combined_pre=${filename}.syntheticBAM.combined.20221011_UMITransLevel2
# Gene region
python3 ${script_directory}/NOINPUT_RNA_GenerateBAMCoord_20221011_UMITranscriptLevel2.py ${outdirectory}/${ref_peakfile} ${synthetic_UMI_count_mat_file} ${BED_filename_pre}.read.bed ${INPUT_bamfile} ${outdirectory} ${OUTPUT_cells_barcode_file}
# Inter-gene region
python3 ${script_directory}/NOINPUT_RNA_GenerateBAMCoord_20221011_UMITranscriptLevel2.py ${outdirectory}/${ref_comple_peakfile} ${synthetic_UMI_count_mat_file_comple} ${BED_COMPLE_filename_pre}.read.bed ${INPUT_bamfile} ${outdirectory} ${OUTPUT_cells_barcode_file}
## Combine gene and inter-gene region
cat ${outdirectory}/${BED_filename_pre}.read.bed ${outdirectory}/${BED_COMPLE_filename_pre}.read.bed | sort -k1,1 -k2,2n > ${outdirectory}/${BED_filename_combined_pre}.read.bed

######################## Generate FASTQ file ######################## 
referenceGenome_name=GRCm38.primary_assembly.genome
referenceGenome_dir=/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE
referenceGenome_file=${referenceGenome_dir}/${referenceGenome_name}.fa
output_BAM_pre=${BED_filename_pre}
time(bedtools getfasta -s -nameOnly -fi ${referenceGenome_file} -bed  ${outdirectory}/${BED_filename_combined_pre}.read.bed -fo ${outdirectory}/${BED_filename_combined_pre}.read2.strand.bed2fa.fa) # 4 mins
sed '/^>/s/.\{3\}$//' ${outdirectory}/${BED_filename_combined_pre}.read2.strand.bed2fa.fa > ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fa # remove (-)
# Create read 1 in fasta 
time(awk 'NR%2==0 {print substr(p,2,26);} NR%2 {p=$0;print p;}' ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fa > ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.fa)
## Convert fasta to fastq
# ~/Tools/seqtk/seqtk/seqtk seq -F '#' ~/Projects/scATAC_Simulator/data/bedTobam_output/BoneMarrow_62016_chr1.bed2fa.fa > ~/Projects/scATAC_Simulator/data/bedTobam_output/BoneMarrow_62016_chr1.bed2fa.fq
time(~/Tools/seqtk/seqtk/seqtk seq -F 'F' ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.fa > ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.fq)
time(~/Tools/seqtk/seqtk/seqtk seq -F 'F' ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fa > ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fq)
time(cat ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq)
time(cat ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq)

######################## Generate BAM file ######################## 
bowtie2_directory=~/Tools/bowtie2
time(bowtie2 -x ${referenceGenome_dir}/${referenceGenome_name} -U ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq -S ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.sam) # 151m37.660s
${samtools_directory}/samtools view -S -b ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.sam > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.bam
# Create CB and UB tag
# ${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.bam | cut -f1 | cut -d':' -f1 | awk '{s=substr($1,1,16)}{g=substr($1,17,length($1))}{printf "CB:Z:%s\tUB:Z:%s\n",s,g;}'
time(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.bam -H > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.header.sam)
time(cat <( cat ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.header.sam ) \
<( paste <(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.bam ) <(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.bam | cut -f1 | cut -d':' -f1 | awk '{s=substr($1,1,16)}{g=substr($1,17,length($1))}{printf "CB:Z:%s\tUB:Z:%s\n",s,g;}')) | ${samtools_directory}/samtools view -bS - > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.bam)
rm ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.header.sam 
# Convert sam to bam
# ${samtools_directory}/samtools view -S -b ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sam > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.bam
${samtools_directory}/samtools sort ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.bam > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam
${samtools_directory}/samtools index ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam




######################## 20221015 Introduce Error using NONINPUT_ATAC_ReadError.py  ######################## 
mkdir -p ${outdirectory}/VerifyCoverage
scReadSim_bamfile=${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.sorted.bam

### Plot 2: Kmer spectrum
mkdir -p ${outdirectory}/VerifyCoverage/kmer_real
# cp /home/gayan/Projects/scATAC_Simulator/results/20220608_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/VerifyCoverage/kmer_real/* ${outdirectory}/VerifyCoverage/kmer_real/
mkdir -p ${outdirectory}/VerifyCoverage/kmer_scReadSim
for kmer_len in 11 21 31 41 51
do 
  jellyfish count -m ${kmer_len} -o ${outdirectory}/VerifyCoverage/kmer_real/mer${kmer_len}_counts_read.jf -s 100M -t 10 -C ${outdirectory}/real_bam2fasta_read.fa
  jellyfish histo ${outdirectory}/VerifyCoverage/kmer_real/mer${kmer_len}_counts_read.jf > ${outdirectory}/VerifyCoverage/kmer_real/mer${kmer_len}_counts_read_hist.txt
  # Kmer on simulated fa
  jellyfish count -m ${kmer_len} -o ${outdirectory}/VerifyCoverage/kmer_scReadSim/mer${kmer_len}_counts_read2.jf -s 100M -t 10 -C ${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.fa
  jellyfish histo ${outdirectory}/VerifyCoverage/kmer_scReadSim/mer${kmer_len}_counts_read2.jf > ${outdirectory}/VerifyCoverage/kmer_scReadSim/mer${kmer_len}_counts_read2_hist.txt
done


### Plot 3: error rates
## fgbio 
mkdir -p ${outdirectory}/VerifyQuality/fgbio/
# cp  /home/gayan/Projects/scATAC_Simulator/results/20220608_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT/VerifyQuality/fgbio/Real* ${outdirectory}/VerifyQuality/fgbio/
java -jar /home/gayan/Tools/picard/build/libs/picard.jar CreateSequenceDictionary \
      -R $referenceGenome_file \
      -O ${referenceGenome_file}.dict
java -jar /home/gayan/Tools/fgbio/target/scala-2.13/fgbio-2.0.1-e884860-SNAPSHOT.jar ErrorRateByReadPosition -i $INPUT_bamfile -r $referenceGenome_file -o ${outdirectory}/VerifyQuality/fgbio/Real --collapse false

## Introduce error rate bam
# python3 ${script_directory}/NOINPUT_RNA_ReadError_20220519.py ${outdirectory}/VerifyQuality/fgbio/Real.error_rate_by_read_position.txt ${outdirectory} ${BED_filename_combined_pre} 8 
python3 ${script_directory}/NOINPUT_RNA_ReadError_20220421.py ${outdirectory}/VerifyQuality/fgbio/Real.error_rate_by_read_position.txt ${outdirectory} ${BED_filename_combined_pre} ${BED_filename_combined_pre}
time(cat ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > ${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq)
time(cat ${outdirectory}/${BED_filename_combined_pre}.ErrorIncluded.read2.bed2fa.fq | paste - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > ${outdirectory}/${BED_filename_combined_pre}.ErrorIncluded.read2.bed2fa.sorted.fq)
## Bowtie align sam
# cd ${referenceGenome_dir}
bowtie2_directory=~/Tools/bowtie2
time(bowtie2 -x ${referenceGenome_dir}/${referenceGenome_name} -U ${outdirectory}/${BED_filename_combined_pre}.ErrorIncluded.read2.bed2fa.sorted.fq -S ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.sam) # 151m37.660s
${samtools_directory}/samtools view -S -b ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.sam > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.bam
## add cell barcode to bam file
time(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.bam -H > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.header.sam)
# time(cat <( cat ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.header.sam ) \
# <( paste <(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.bam ) <(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.bam | cut -f1 | cut -d':' -f1 | sed -e 's/^/CB:Z:/')) | ${samtools_directory}/samtools view -bS - > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.ErrorIncluded.bam)
# rm ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.header.sam 
time(cat <( cat ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.header.sam ) \
<( paste <(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.bam ) <(${samtools_directory}/samtools view ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.bam | cut -f1 | cut -d':' -f1 | awk '{s=substr($1,1,16)}{g=substr($1,17,length($1))}{printf "CB:Z:%s\tUB:Z:%s\n",s,g;}')) | ${samtools_directory}/samtools view -bS - > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.ErrorIncluded.bam)
rm ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.noCB.ErrorIncluded.header.sam 
## Sort and index bam 
${samtools_directory}/samtools sort ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.ErrorIncluded.bam > ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.ErrorIncluded.sorted.bam
${samtools_directory}/samtools index ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.ErrorIncluded.sorted.bam
# Problematic one sorted bam 1974199293

## fgbio 
# java -jar /home/gayan/Tools/fgbio/target/scala-2.13/fgbio-2.0.1-e884860-SNAPSHOT.jar ErrorRateByReadPosition -i $scReadSim_bamfile -r $referenceGenome_fissle -o ${outdirectory}/VerifyQuality/fgbio/scReadSim --collapse false
java -jar /home/gayan/Tools/fgbio/target/scala-2.13/fgbio-2.0.1-e884860-SNAPSHOT.jar ErrorRateByReadPosition -i ${outdirectory}/${output_BAM_pre}.synthetic.${chr_title}.ErrorIncluded.sorted.bam -r $referenceGenome_file -o ${outdirectory}/VerifyQuality/fgbio/scReadSim_manualError --collapse false

