#-----------------------------------------------------#
# NEEDS ONE ARGUMENT : PIPELINE_HOME                  #
#-----------------------------------------------------#
if [[ $# -eq 1 ]]; then
  export PIPELINE_HOME=$1
else
  echo "ERROR INITIALIZING PIPELINE: PROVIDE VALID PIPELINE_HOME"
  exit 1
fi
#-----------------------------------------------------#
# System Defaults
umask 002
export PATH=${PIPELINE_HOME}/soft/bin/:$PATH
export JAVA_MEM=4
#export MAX_MEM=46 ; #this should be set in a server environment script
#export THREADS=12 ; #same
export EXOME_DOWNSAMPLE=1000

#-----------------------------------------------------#
# Analysis files and paths
#export GATK_TORQUE_DIR=${PIPELINE_HOME}/soft/packages/gatk07_24_12_torque
#export GATK_DIR=${PIPELINE_HOME}/soft/packages/GATK-git-2.6-4
#export GATK_DIR=${PIPELINE_HOME}/soft/packages/GenomeAnalysisTK-3.5
#export GATK_DIR=${PIPELINE_HOME}/soft/packages/GenomeAnalysisTK-3.6
#export GATK_DIR=/RQexec/dionnela/soft/packages//GenomeAnalysisTK-nightly-2016-11-25-gf5379ae-3.6
#export GATK_DIR=${PIPELINE_HOME}/soft/packages/GenomeAnalysisTK-3.7
#export GATK_DIR=${PIPELINE_HOME}/soft/packages/GenomeAnalysisTK-3.8
export GATK_DIR=${PIPELINE_HOME}/soft/packages/gatk
export PICARD_DIR=${PIPELINE_HOME}/soft/packages/picard
#export BWA_DIR=${PIPELINE_HOME}/soft/src/bwa
export BWA_DIR=${PIPELINE_HOME}/soft/src/bwa
export SAMTOOLS_DIR=${PIPELINE_HOME}/soft/src/samtools
export TRIMMOMATIC_JAR=${PIPELINE_HOME}/soft/packages/trimmomatic/trimmomatic.jar
#-----------------------------------------------------#
export VERIFYBAMID_GENOTYPES=${PIPELINE_HOME}/soft/packages/verifyBamID/data/Omni25_genotypes.1525_samples_v2.b37.PASS.ALL_sites.in_concat_38mb_50mb_V4.vcf.gz
export DBSNP=${PIPELINE_HOME}/data/reference/dbsnp_132.b37.vcf
export INDELS=${PIPELINE_HOME}/data/reference/1000G_indels_for_realignment.b37.vcf
export HAPMAP=${PIPELINE_HOME}/data/reference/hapmap_3.3.b37.sites.vcf
export OMNI=${PIPELINE_HOME}/data/reference/1000G_omni2.5.b37.sites.vcf
export MILLS_KG=${PIPELINE_HOME}/data/reference/Mills_and_1000G_gold_standard.indels.b37.sites.vcf
export REF=${PIPELINE_HOME}/data/reference/human_g1k_v37.fasta
export REF_BWA_SQ_HEADER=${PIPELINE_HOME}/data/reference/human_g1k_v37.fasta.bwa_sq_header
export REF_MT=${PIPELINE_HOME}/data/reference/human_g1k_v37_MT.fasta
#export REF=${PIPELINE_HOME}/data/reference/human_g1k_v37.bwa-0_7_5a.fasta
#export REF_MT=${PIPELINE_HOME}/data/reference/human_g1k_v37_MT.bwa-0_7_5a.fasta
export BED38=${PIPELINE_HOME}/data/targets/SureSelect_All_Exon_38mb.v37.bed
export BED38_FLANK=${PIPELINE_HOME}/data/targets/SureSelect_All_Exon_38mb.v37.500_flanking.bed
export BED50=${PIPELINE_HOME}/data/targets/SureSelect_All_Exon_50mb.v37.bed
export BED50_FLANK=${PIPELINE_HOME}/data/targets/SureSelect_All_Exon_50mb.v37.500_flanking.bed
export BEDV4=${PIPELINE_HOME}/data/targets/SureSelect_All_Exon_V4.v37.bed
export BEDV4_FLANK=${PIPELINE_HOME}/data/targets/SureSelect_All_Exon_V4.v37.500_flanking.bed
export BEDV5=${PIPELINE_HOME}/data/targets/SureSelect_All_Exon_V5.v37.bed
export BEDV5_FLANK=${PIPELINE_HOME}/data/targets/SureSelect_All_Exon_V5.v37.500_flanking.bed
export BEDV5CLINIC=${PIPELINE_HOME}/data/targets/SureSelect_All_Exon_V5_Clinic.v37.bed
export BEDV5CLINIC_FLANK=${PIPELINE_HOME}/data/targets/SureSelect_All_Exon_V5_Clinic.v37.500_flanking.bed
export BEDV7=${PIPELINE_HOME}/data/targets/SureSelect_All_Exon_V7.v37.bed
export BEDV7_FLANK=${PIPELINE_HOME}/data/targets/SureSelect_All_Exon_V7.v37.500_flanking.bed
export BEDTS=${PIPELINE_HOME}/data/targets/TruSeq_Exome_Targeted_Regions.b37.bed
export BEDTS_FLANK=${PIPELINE_HOME}/data/targets/TruSeq_Exome_Targeted_Regions.b37.500_flanks.bed
export BEDNEXTERA_10=${PIPELINE_HOME}/data/targets/Illumina_Nextera_Rapid_Capture_Exome_Targeted_Regions_1.0_37Mbps_2013-03-07.v37.bed
export BEDNEXTERA_10_FLANK=${PIPELINE_HOME}/data/targets/Illumina_Nextera_Rapid_Capture_Exome_Targeted_Regions_1.0_37Mbps_2013-03-07.v37.500_flanking.bed
export BEDNEXTERA_12=${PIPELINE_HOME}/data/targets/Illumina_Nextera_Rapid_Capture_Exome_Targeted_Regions_1.2_54Mbps_2014-03-12.v37.bed
export BEDNEXTERA_12_FLANK=${PIPELINE_HOME}/data/targets/Illumina_Nextera_Rapid_Capture_Exome_Targeted_Regions_1.2_54Mbps_2014-03-12.v37.500_flanking.bed
export BEDNGvCRome=${PIPELINE_HOME}/data/targets/Nimblegen_SeqCap_HGSC_VCRome_Design_Regions.v37.bed
export BEDNGvCRome_FLANK=${PIPELINE_HOME}/data/targets/Nimblegen_SeqCap_HGSC_VCRome_Design_Regions.v37.500_flanking.bed
export BEDNGv3=${PIPELINE_HOME}/data/targets/Nimblegen_SeqCap_EZ_Exome_v3_primary.v37.bed
export BEDNGv3_FLANK=${PIPELINE_HOME}/data/targets/Nimblegen_SeqCap_EZ_Exome_v3_primary.v37.500_flanking.bed
#export BEDALLKITS=${PIPELINE_HOME}/data/targets/CCDS-RefSeq_Genes-SS_38mb_50mb_V4_V5_V5_Clinic-NG_V3_VCRome-TruSeq_Nextera_1_2.v37.bed
#export BEDALLKITS_FLANK=${PIPELINE_HOME}/data/targets/CCDS-RefSeq_Genes-SS_38mb_50mb_V4_V5_V5_Clinic-NG_V3_VCRome-TruSeq_Nextera_1_2.v37.500_flanking.bed
export BEDALLKITS=${PIPELINE_HOME}/data/targets/CCDS-RefSeq_Genes-SS_38mb_50mb_V4_V5_V5_Clinic_V7-NG_V3_VCRome-TruSeq_Nextera_1_2.v37.bed
export BEDALLKITS_FLANK=${PIPELINE_HOME}/data/targets/CCDS-RefSeq_Genes-SS_38mb_50mb_V4_V5_V5_Clinic_V7-NG_V3_VCRome-TruSeq_Nextera_1_2.v37.500_flanking.bed
export BEDCCDS=${PIPELINE_HOME}/data/targets/CCDS.v37.bed
export BEDCCDS_FLANK=${PIPELINE_HOME}/data/targets/CCDS.v37.500_flanking.bed
export CALLING_BED=${BEDALLKITS_FLANK}
export GENE_BED=${PIPELINE_HOME}/data/targets/RefSeq.Genes.v37.coding_exons.no_flanks.sort.bed
export GENE_LIST=${PIPELINE_HOME}/data/targets/RefSeq.Genes.v37.refGene
#-----------------------------------------------------#
export PERL5LIB=${PIPELINE_HOME}/soft/packages/perl/share/perl5
export CMHTOOLS=${PIPELINE_HOME}/soft/packages/cmhtools/
#-----------------------------------------------------#
# GATK Scala files
export EXOME_POST_PROCESSING_SCALA=${PIPELINE_HOME}/data/templates/DataProcessingPipeline.noTorque.exome.fastqs.scala
export CUSTOM_POST_PROCESSING_SCALA_NO_RECAL=${PIPELINE_HOME}/data/templates/DataProcessingPipeline.noTorque.custom.no_recal.fastqs.scala
export WGS_ALIGNED_BAM_POST_PROCESSING_SCALA=${PIPELINE_HOME}/data/templates/DataProcessingPipeline.noTorque.genome.aligned_bam.scala
export GATK_KEY=${PIPELINE_HOME}/data/templates/dan.spiegelman_crchum.qc.ca.key
#-----------------------------------------------------#
# GATK results files
export GATK_WGS_PREFIX=WGS
export GATK_EXOME_PREFIX=Exome
export GATK_CUSTOM_PREFIX=Custom
export GATK_MT_PREFIX=MT
export GATK_WGS_BAMFILE=${GATK_WGS_PREFIX}.${SAMPLE}.clean.dedup.recal.bam
export GATK_EXOME_BAMFILE=${GATK_EXOME_PREFIX}.${SAMPLE}.clean.dedup.recal.bam
export GATK_CUSTOM_BAMFILE=${GATK_CUSTOM_PREFIX}.${SAMPLE}.clean.dedup.recal.bam
export GATK_CUSTOM_BAMFILE_NO_RECAL=${GATK_CUSTOM_PREFIX}.${SAMPLE}.clean.dedup.bam
export GATK_HC_BAMFILE=${GATK_EXOME_PREFIX}.${SAMPLE}.GATK_HC.bam
export GATK_CUSTOM_HC_BAMFILE=${GATK_CUSTOM_PREFIX}.${SAMPLE}.GATK_HC.bam
export GATK_EXOME_HC_VCFFILE=${GATK_EXOME_PREFIX}.${SAMPLE}.GATK_HC.variant_sites.vcf.gz
export GATK_CUSTOM_HC_VCFFILE=${GATK_CUSTOM_PREFIX}.${SAMPLE}.GATK_HC.variant_sites.vcf.gz
export GATK_EXOME_HC_GVCFFILE=${GATK_EXOME_PREFIX}.${SAMPLE}.GATK_HC.g.vcf.gz
export GATK_CUSTOM_HC_GVCFFILE=${GATK_CUSTOM_PREFIX}.${SAMPLE}.GATK_HC.g.vcf.gz
export GATK_EXOME_TOTAL_RECALL_VCFFILE=${GATK_EXOME_PREFIX}.${SAMPLE}.GATK_UG.vcf.gz
export GATK_CUSTOM_TOTAL_RECALL_VCFFILE=${GATK_CUSTOM_PREFIX}.${SAMPLE}.GATK_UG.vcf.gz
export GATK_EXOME_TOTAL_RECALL_VARIANT_SITES_VCFFILE=${GATK_EXOME_PREFIX}.${SAMPLE}.GATK_UG.variant_sites.vcf.gz
export GATK_CUSTOM_TOTAL_RECALL_VARIANT_SITES_VCFFILE=${GATK_CUSTOM_PREFIX}.${SAMPLE}.GATK_UG.variant_sites.vcf.gz
export GATK_EXOME_GVCFFILE=${GATK_EXOME_PREFIX}.${SAMPLE}.GATK_HC.g.vcf.gz
export GATK_CUSTOM_GVCFFILE=${GATK_CUSTOM_PREFIX}.${SAMPLE}.GATK_HC.g.vcf.gz
export GATK_WGS_VCFFILE=${GATK_WGS_PREFIX}.${SAMPLE}.GATK_UG.vcf
export GATK_MT_BAMFILE=${GATK_MT_PREFIX}.${SAMPLE}.clean.dedup.recal.bam
#export GATK_BAMFILE_MT_FILTERED=${GATK_MT_PREFIX}.${SAMPLE}.ref_MT.MT_L30.bam
export GATK_BAMFILE_MT_FILTERED=${GATK_MT_PREFIX}.${SAMPLE}.clean.dedup.recal.filtered.bam
#-----------------------------------------------------#
