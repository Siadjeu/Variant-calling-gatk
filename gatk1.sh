#!/bin/bash
#SBATCH --partition mpcb.p           # or mpcp.p
#SBATCH --ntasks 1                          # instead of node
#SBATCH --cpus-per-task 16           # use 40 with mpcp.p
##SBATCH --exclusive                      # not needed if you ask for all the cores in a node
#SBATCH --time 2-0:00:00
##SBATCH --mem-per-cpu 30000  # not needed if you use all the cores, you also get all memory per default
#SBATCH --job-name tunnel
#SBATCH --output gatk4-log-%J.txt


# Load module

ml hpc-env/8.3
ml picard/2.23.6-Java-8
ml GATK/4.1.9.0-GCCcore-8.3.0-Java-8
ml SAMtools/1.9-GCC-8.3.0
ml BWA/0.7.17-GCC-8.3.0
ml R/4.0.2-foss-2019b

## Alignment – Map to Reference

# index reference

bwa index Yam_ref.fasta

## --- adding read group IDs --- #
bwa mem -t 16 -R '@RG\tID:Yam_L001\tLB:Yam_L001\tPL:ILLUMINA\tPM:HISEQ\tSM:Yam_L001' Dioscorea_dumetorum_v1.0.fasta /gss/work/wupf2892/trimodata/Yam_L001_R1paired.fastq.gz /gss/work/wupf2892/trimodata/Yam_L001_R2paired.fastq.gz > Yam_ref_Aligned.reads.sam

samtools faidx Yam_ref.fasta

samtools dict -o Yam_ref.dict Yam_ref.fasta

gatk MarkDuplicatesSpark \
       -I Yam_ref_Aligned.reads.sam \
       -M dedup_metrics.txt \
       -O Yam_sorted_dedup_reads.bam 

java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics -R Yam_ref.fasta -I Yam_sorted_dedup_reads.bam -O alignment_metrics.txt

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics --I Yam_sorted_dedup_reads.bam  --O insert_size_metrics.txt --Histogram_FILE insert_size_histogram.pdf

samtools depth -a Yam_sorted_dedup_reads.bam > depth_out.txt

# Call Variants
gatk HaplotypeCaller \
       -R Yam_ref.fasta \
       -I Yam_sorted_dedup_reads.bam \
       -O Yam_raw_variants.vcf

# Extract SNPs and Indels

gatk SelectVariants \
       -R Yam_ref.fasta \
       -V Yam_raw_variants.vcf \
       --select-type-to-include SNP \
       -O Yam_raw_snps.vcf

gatk SelectVariants \
      -R Yam_ref.fasta \
      -V Yam_raw_variants.vcf \
       --select-type-to-include INDEL \
       -O Yam_raw_indels.vcf

# Filter SNPs

 
gatk VariantFiltration \
        -R Yam_ref.fasta \
        -V Yam_raw_snps.vcf \
        -O Yam_filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

# Filter indels

 gatk VariantFiltration \
  #      -R Yam_ref.fasta \
   #     -V Yam_raw_indels.vcf \
    #    -O Yam_filtered_indels.vcf \
     #   -filter-name "QD_filter" -filter "QD < 2.0" \
      #  -filter-name "FS_filter" -filter "FS > 200.0" \
       # -filter-name "SOR_filter" -filter "SOR > 10.0" 

# Exclude Filtered variants

 gatk SelectVariants \
        --exclude-filtered \
        -V Yam_filtered_snps.vcf \
        -O Yam_bqsr_snps.vcf

 gatk SelectVariants \
       --exclude-filtered \
       -V Yam_filtered_indels.vcf \
       -O Yam_bqsr_indels.vcf


# Base quality score recalibration(BQSR)

 gatk BaseRecalibrator \
        -R Yam_ref.fasta \
        -I Yam_sorted_dedup_reads.bam \
        --known-sites Yam_bqsr_snps.vcf \
        --known-sites Yam_bqsr_indels.vcf \
        -O recal_data.table 

# Apply BQSR

gatk ApplyBQSR \
        -R Yam_ref.fasta \
        -I Yam_sorted_dedup_reads.bam \
        -bqsr recal_data.table \
        -O recal_reads.bam 

# Base Quality Score Recalibration (BQSR) #2

gatk BaseRecalibrator \
        -R Yam_ref.fasta \
        -I recal_reads.bam \
        --known-sites Yam_bqsr_snps.vcf \
        --known-sites Yam_bqsr_indels.vcf \
        -O post_recal_data.table 

# Analyze Covariates

gatk AnalyzeCovariates -before recal_data.table \
        -after post_recal_data.table \
        -plots recalibration_plots.pdf

# Call variants

gatk HaplotypeCaller \
        -R Yam_ref.fasta \
        -I recal_reads.bam \
        -O Yam_raw_variants_recal.vcf

# Extract SNPs & Indels

gatk SelectVariants \
        -R Yam_ref.fasta \
        -V Yam_raw_variants_recal.vcf \
        --select-type-to-include SNP \
        -O Yam_raw_snps_recal.vcf

gatk SelectVariants \
        -R Yam_ref.fasta \
        -V Yam_raw_variants.vcf \
        --select-type-to-include INDEL \
        -O Yam_raw_indels_recal.vcf

# Filter SNPs

gatk VariantFiltration -R Yam_ref.fasta \
       -V Yam_raw_snps_recal.vcf \
        -O Yam_filtered_snps_final.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

gatk VariantFiltration -R Yam_ref.fasta \
        -V Yam_raw_indels_recal.vcf \
        -O Yam_filtered_indels_final.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0" 
  
