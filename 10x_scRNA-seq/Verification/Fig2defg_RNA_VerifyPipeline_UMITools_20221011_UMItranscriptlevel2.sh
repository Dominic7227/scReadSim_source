#!/bin/bash
full_filename=e18_mouse_brain_fresh_5k_gex_possorted_bam
full_directory=/home/gayan/Projects/scATAC_Simulator/data/10X_MOUSE_BRAIN_ATACandRNA

referenceGenome_dir=/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE


chr_title=chr1
filename=${full_filename}_${chr_title}
outdirectory=/home/gayan/Projects/scATAC_Simulator/results/20221011_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT
BED_filename_combined_pre=${filename}.syntheticBAM.combined.20221011_UMITransLevel2
synthetic_read1_fq=${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq
synthetic_read2_fq=${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq

## Test
# ${outdirectory}/${filename}.syntheticBAM.20221011_UMITransLevel2.read1.bed2fa.sorted.fq
# umi_tools whitelist --stdin ${outdirectory}/${filename}.syntheticBAM.20221011_UMITransLevel2.read1.bed2fa.sorted.fq \
#                     --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
#                     --set-cell-number=4877 \
#                      -L ${umitools_verification_dir}/extract.log > ${umitools_verification_dir}/whitelist.txt;


umitools_verification_dir=${outdirectory}/verification_umitools
mkdir ${umitools_verification_dir}

export PATH=${PATH}:/home/gayan/.local/bin
# Step 2: Identify correct cell barcodes
umi_tools whitelist --stdin ${synthetic_read1_fq} \
                    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                    --set-cell-number=4877 \
                     -L ${umitools_verification_dir}/extract.log > ${umitools_verification_dir}/whitelist.txt;
                    
# Step 3: Extract barcdoes and UMIs and add to read names
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin ${synthetic_read1_fq} \
                  --stdout ${umitools_verification_dir}/${BED_filename_combined_pre}.read1.bed2fa.sorted.extracted.fq \
                  --read2-in ${synthetic_read2_fq} \
                  --read2-out=${umitools_verification_dir}/${BED_filename_combined_pre}.read2.bed2fa.sorted.extracted.fq \
                  -L ${umitools_verification_dir}/extract.log;

# Step 4: Map reads
samtools_directory=/home/gayan/Tools/samtools/bin
export PATH=$PATH:${samtools_directory}:/home/gayan/Tools/STAR/bin/Linux_x86_64_static
genomeDir_STAR=${referenceGenome_dir}/STAR
cd ${umitools_verification_dir}
# mkdir ${genomeDir_STAR}
# STAR --runThreadN 4 \
# 	 --runMode genomeGenerate \
# 	 --genomeDir ${genomeDir_STAR} \
# 	 --genomeFastaFiles ${referenceGenome_dir}/GRCm38.primary_assembly.genome.fa \
# 	 --sjdbGTFfile ${referenceGenome_dir}/gencode.vM10.primary_assembly.annotation.gtf \
# 	 --sjdbOverhang 79

STAR --runThreadN 4 \
     --genomeDir ${genomeDir_STAR} \
     --readFilesIn ${umitools_verification_dir}/${BED_filename_combined_pre}.read2.bed2fa.sorted.extracted.fq \
     --readFilesCommand cat \
     --outFilterMultimapNmax 1 \
     --outSAMtype BAM SortedByCoordinate;

# bwa_directory=~/Tools/bwa/
# time(${bwa_directory}/bwa mem ${referenceGenome_file} ${umitools_verification_dir}/${BED_filename_combined_pre}.read2.bed2fa.sorted.extracted.fq > ${umitools_verification_dir}/${output_BAM_pre}.synthetic.${chr_title}.umi_tools.bam)

# Step 5: Assign reads to genes
referenceGenome_annotationfile=${referenceGenome_dir}/gencode.vM10.primary_assembly.annotation.gtf
time(featureCounts -a ${referenceGenome_annotationfile} \
              -o ${umitools_verification_dir}/gene_assigned \
              -R BAM ${umitools_verification_dir}/Aligned.sortedByCoord.out.bam \
              -T 4;           ) 
samtools sort Aligned.sortedByCoord.out.bam.featureCounts.bam -o assigned_sorted.bam
samtools index assigned_sorted.bam

# Step 6: Count UMIs per gene per cell
umi_tools count --per-gene --gene-tag=XT --wide-format-cell-counts --assigned-status-tag=XS --per-cell -I assigned_sorted.bam -S counts.tsv
            