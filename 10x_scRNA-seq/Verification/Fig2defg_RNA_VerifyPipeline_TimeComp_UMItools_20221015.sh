#!/bin/bash
full_filename=e18_mouse_brain_fresh_5k_gex_possorted_bam
full_directory=/home/gayan/Projects/scATAC_Simulator/data/10X_MOUSE_BRAIN_ATACandRNA

referenceGenome_dir=/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE
referenceGenome_dir=/home/gayan/Projects/scATAC_Simulator/data/mm10_ref_genome_GECODE
referenceGenome_file=${referenceGenome_dir}/${referenceGenome_name}.fa


chr_title=chr1
filename=${full_filename}_${chr_title}
outdirectory=/home/gayan/Projects/scATAC_Simulator/results/20221013_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_VaryNewCell

umitools_verification_dir=${outdirectory}/verification_umitools
mkdir ${umitools_verification_dir}
export PATH=${PATH}:/home/gayan/.local/bin
samtools_directory=/home/gayan/Tools/samtools/bin
export PATH=$PATH:${samtools_directory}:/home/gayan/Tools/STAR/bin/Linux_x86_64_static
genomeDir_STAR=${referenceGenome_dir}/STAR

###################################################################
######################## VarySeqDepth ################################
###################################################################
timer_umitools_VarySeqDepth (){
local adj_factor=$1
BED_filename_pre=${filename}.syntheticBAM.VarySeqDepth${adj_factor}
BED_COMPLE_filename_pre=${filename}.syntheticBAM.VarySeqDepth${adj_factor}.COMPLE
BED_filename_combined_pre=${filename}.syntheticBAM.VarySeqDepth${adj_factor}.combined
output_BAM_pre=${BED_filename_pre}
verification_dir=${umitools_verification_dir}/Case_VarySeqDepth${adj_factor}
mkdir ${verification_dir}

synthetic_read1_fq=${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq
synthetic_read2_fq=${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq

START=$(date +%s)
# Step 1: Whitelist
umi_tools whitelist --stdin ${synthetic_read1_fq} \
                    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                    --set-cell-number=4877 \
                     -L ${verification_dir}/extract.log > ${verification_dir}/whitelist.txt;
                    
# Step 3: Extract barcdoes and UMIs and add to read names
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin ${synthetic_read1_fq} \
                  --stdout ${verification_dir}/${BED_filename_combined_pre}.read1.bed2fa.sorted.extracted.fq \
                  --read2-in ${synthetic_read2_fq} \
                  --read2-out=${verification_dir}/${BED_filename_combined_pre}.read2.bed2fa.sorted.extracted.fq \
                  -L ${verification_dir}/extract.log;

# Step 4: Map reads
cd ${verification_dir}
STAR --runThreadN 4 \
     --genomeDir ${genomeDir_STAR} \
     --readFilesIn ${verification_dir}/${BED_filename_combined_pre}.read2.bed2fa.sorted.extracted.fq \
     --readFilesCommand cat \
     --outFilterMultimapNmax 1 \
     --outSAMtype BAM SortedByCoordinate;

# Step 5: Assign reads to genes
referenceGenome_annotationfile=${referenceGenome_dir}/gencode.vM10.primary_assembly.annotation.gtf
time(featureCounts -a ${referenceGenome_annotationfile} \
              -o ${verification_dir}/gene_assigned \
              -R BAM ${verification_dir}/Aligned.sortedByCoord.out.bam \
              -T 4;           ) 
samtools sort Aligned.sortedByCoord.out.bam.featureCounts.bam -o assigned_sorted.bam
samtools index assigned_sorted.bam

# Step 6: Count UMIs per gene per cell
umi_tools count --per-gene --gene-tag=XT --wide-format-cell-counts --assigned-status-tag=XS --per-cell -I assigned_sorted.bam -S counts.tsv
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "[Task VarySeq ${adj_factor}] $DIFF seconds"
echo "[Task VarySeq ${adj_factor}] Lapsed $DIFF seconds" >> ${outdirectory}/verification_umitools/TIMER_VarySeqDepth.txt
}          
echo $(date) >> ${outdirectory}/verification_umitools/TIMER_VarySeqDepth.txt
for adj_factor in 0.25 0.5 1 2 4
do 
     timer_umitools_VarySeqDepth "$adj_factor"
done
# [1] 8722
# [2] 8723
# [3] 8724
# [4] 8725
# [5] 8727


###################################################################
######################## VaryCell ################################
###################################################################
timer_umitools_VaryCell (){
local adj_factor=$1
BED_filename_pre=${filename}.syntheticBAM.VaryCellNumber${adj_factor}
BED_COMPLE_filename_pre=${filename}.syntheticBAM.VaryCellNumber${adj_factor}.COMPLE
BED_filename_combined_pre=${filename}.syntheticBAM.VaryCellNumber${adj_factor}.combined
output_BAM_pre=${BED_filename_pre}
verification_dir=${umitools_verification_dir}/Case_VaryCell${adj_factor}
mkdir ${verification_dir}

synthetic_read1_fq=${outdirectory}/${BED_filename_combined_pre}.read1.bed2fa.sorted.fq
synthetic_read2_fq=${outdirectory}/${BED_filename_combined_pre}.read2.bed2fa.sorted.fq

START=$(date +%s)
# Step 1: Whitelist
umi_tools whitelist --stdin ${synthetic_read1_fq} \
                    --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                    --set-cell-number=4877 \
                     -L ${verification_dir}/extract.log > ${verification_dir}/whitelist.txt;
                    
# Step 3: Extract barcdoes and UMIs and add to read names
umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN \
                  --stdin ${synthetic_read1_fq} \
                  --stdout ${verification_dir}/${BED_filename_combined_pre}.read1.bed2fa.sorted.extracted.fq \
                  --read2-in ${synthetic_read2_fq} \
                  --read2-out=${verification_dir}/${BED_filename_combined_pre}.read2.bed2fa.sorted.extracted.fq \
                  -L ${verification_dir}/extract.log;

# Step 4: Map reads
cd ${verification_dir}
STAR --runThreadN 4 \
     --genomeDir ${genomeDir_STAR} \
     --readFilesIn ${verification_dir}/${BED_filename_combined_pre}.read2.bed2fa.sorted.extracted.fq \
     --readFilesCommand cat \
     --outFilterMultimapNmax 1 \
     --outSAMtype BAM SortedByCoordinate;

# Step 5: Assign reads to genes
referenceGenome_annotationfile=${referenceGenome_dir}/gencode.vM10.primary_assembly.annotation.gtf
time(featureCounts -a ${referenceGenome_annotationfile} \
              -o ${verification_dir}/gene_assigned \
              -R BAM ${verification_dir}/Aligned.sortedByCoord.out.bam \
              -T 4;           ) 
samtools sort Aligned.sortedByCoord.out.bam.featureCounts.bam -o assigned_sorted.bam
samtools index assigned_sorted.bam

# Step 6: Count UMIs per gene per cell
umi_tools count --per-gene --gene-tag=XT --wide-format-cell-counts --assigned-status-tag=XS --per-cell -I assigned_sorted.bam -S counts.tsv
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "[Task VaryCell ${adj_factor}] $DIFF seconds"
echo "[Task VaryCell ${adj_factor}] Lapsed $DIFF seconds" >> ${outdirectory}/verification_umitools/TIMER_VaryCell.txt
}          

echo $(date) >> ${outdirectory}/verification_umitools/TIMER_VaryCell.txt
for adj_factor in 0.25 0.5 1 2 4
     do timer_umitools_VaryCell "$adj_factor"
done
# [1] 9379
# [2] 9380
# [3] 9381
# [4] 9382
# [5] 9384











