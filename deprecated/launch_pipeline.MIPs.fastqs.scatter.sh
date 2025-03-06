#!/bin/bash

# ===============================================
# default variables values
# ===============================================
unset SAMPLE DIR THREADS MAX_MEM PIPELINE_HOME WALLTIME EXTRA_CONF CAPTURE_KIT RUN1F1 RUN1F2 RUN1RG RUN2F1 RUN2F2 RUN2RG RUN3F1 RUN3F2 RUN3RG RUN4F1 RUN4F2 RUN4RG RUN5F1 RUN5F2 RUN5RG RUN6F1 RUN6F2 RUN6RG RUN7F1 RUN7F2 RUN7RG RUN8F1 RUN8F2 RUN8RG RUN9F1 RUN9F2 RUN9RG  QUEUE ACCOUNT VERBOSE CENTER SFS KEEP_INTERMEDIATES 

#Assuming script is in root of PIPELINE_HOME
export PIPELINE_HOME=$(cd $(dirname $0) && pwd -P)
export SFS=0
export KEEP_INTERMEDIATES=0
export DISCARD_UNPAIRED_READS=0


#set queue-specific values, depending on system check
QUEUE=qsub                                                                           # default job scheduler: qsub
if [[ `which msub 2>/dev/null` ]]; then QUEUE=msub; fi                               # assign $QUEUE based on server's scheduler software
if [[ `uname -n` =~ ^gpc ]]; then QUEUE=qsub; fi;                                    # fix for scinet gpc system, which has msub, but diabled it in favor of qsub (i.e. `which msub` test will yield TRUE, and mess everything up)
if [[ $QUEUE == "msub" ]]; then dependMod="x=depend:"; else dependMod="depend="; fi  # assigns scheduler-specific dependency syntax, based on scheduler software 

# create function to handle error messages
# ===============================================
Usage() {
	echo
	echo -e "Usage:\t$0 [arguments]"
	echo -e "\tmandatory arguments:\n" \
          "\t\t-s  (--sample)            = sample name\n" \
          "\t\t-d  (--dir)               = run directory (where all the outputs will be printed)\n" \
          "\t\t--r[1-9]f1 (--run[1-9]f1) = fastq sequence run [1-9] forward\n" \
          "\t\t--r[1-9]f2 (--run[1-9]f2) = fastq sequence run [1-9] reverse\n" \
          "\t\t--r[1-9]rg (--run[1-9]rg) = fastq sequence run [1-9] read group\n" \
          "\t\t-f  (--mip_design_file)   = mip design file (mipgen output)\n" \
          "\t\t-x  (--extra)             = use extra config file to (re)define variables [CURRENT \"$EXTRA_CONF\"]" 
	echo -e "\toptional arguments:\n " \
          "\t\t-h  (--help)              = get the program options and exit\n" \
          "\t\t-c  (--capture_kit)       = capture kit used for exome targeting (38,50,V4 or TruSeq)\n" \
          "\t\t--bed                     = specify the bed file for variant calling and QC metrics\n" \
          "\t\t--bed_flank               = specify the bed file with flanking for variant calling and QC metrics\n" \
          "\t\t--sites                   = specify the sites to use for variant calling\n" \
          "\t\t-a  (--adapters)          = fasta sequence file of adapters\n" \
          "\t\t--molecular_tag_length    = length of molecular tag added to each read pair" \	#PICK UP HERE - ADD THIS VARIABLE THROUGHOUT SCRIPT
          "\t\t--tagged_read_number      = read number (1 or 2) of read with molecular tag" \
          "\t\t-m  (--max_mem)           = max memory on a node (in gigs) to use [CURRENT \"$MAX_MEM\"]\n" \
          "\t\t-t  (--threads)           = number of threads to use [CURRENT \"$THREADS\"]\n" \
          "\t\t-w  (--walltime)          = wall time for the scheduler (format HH:MM:SS) [CURRENT \"$WALLTIME\"]\n" \
          "\t\t-v  (--verbose)           = set verbosity level [CURRENT \"$VERBOSE\"]\n" \
          "\t\t--center                  = specify sequencing center [CURRENT \"$CENTER\"]\n" \
          "\t\t--start_from_scratch      = flag to overwrite existing files [CURRENT \"$SFS\"]\n" \
          "\t\t--discard_unpaired        = flag to discard unpaired reads after pre processing [CURRENT \"$DISCARD_UNPAIRED_READS\"]\n" \
          "\t\t--keep_intermediates      = flag to keep intermediates files [CURRENT \"$KEEP_INTERMEDIATES\"]" 
        echo
}


# ===============================================
# alignment function (bwa mem)
# ===============================================
align() {
  local F1=$1
  local F2=$2
  local MIP=$3
  local RG=$4
  local SAM=$5
  local RUN=$6
  local ID=${RG}-${RUN}
  unset DEPEND_ALN;

#DETERMINE BASE
  export ENCODING=`${PIPELINE_HOME}/soft/bin/determine_base.perl $F1`;
  if [[ "$ENCODING" == "Cannot determine base" ]]; then echo "ERROR: Cannot determine base for $F1"; exit 42;
  else echo "*** ENCODING bases set to $ENCODING"; fi

  echo -e "*** ALIGNMENT \n*** R1:$F1\n*** R2:$F2\n*** RUN ID: $ID"

#PRE PREPROCESSING
  FASTQC_OUTPUT_DIR=fastqc.${SAMPLE}-${RUN}
  export TRIM_PAIRED_FASTQ1=${SAMPLE}-${RUN}.r1.fastq.trimmomaticPE.paired.gz
  export TRIM_PAIRED_FASTQ2=${SAMPLE}-${RUN}.r2.fastq.trimmomaticPE.paired.gz
  export TRIM_UNPAIRED_FASTQ1=${SAMPLE}-${RUN}.r1.fastq.trimmomaticPE.unpaired.gz
  export TRIM_UNPAIRED_FASTQ2=${SAMPLE}-${RUN}.r2.fastq.trimmomaticPE.unpaired.gz
  if [[ -n $MOLECULAR_TAG_LENGTH && -n $TAGGED_READ_NUMBER ]]; then MOLECULAR_TAG_ARGS="MOLECULAR_TAG_LENGTH=$MOLECULAR_TAG_LENGTH,TAGGED_READ_NUMBER=$TAGGED_READ_NUMBER,"; else MOLECULAR_TAG_ARGS=""; fi
  export MOLECULAR_TAG_ARGS

  NAME=pipeline.MIP.part_0_1.preprocessing
  export DONE=.${SAMPLE}-${RUN}.${NAME}.done
  #touch $DONE
  if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
    preprocess=`$QUEUE $ACCOUNT  \
      -l nodes=1:ppn=${THREADS} \
      -l mem=${MAX_MEM}g \
      -l walltime=$WALLTIME \
      -N ${SAMPLE}-${RUN}.${NAME} \
      -V \
      -v FASTQ1=$F1,FASTQ2=$F2,MIP_DESIGN_FILE=$MIP_DESIGN_FILE,${MOLECULAR_TAG_ARGS}TRIM_PAIRED_FASTQ1=$TRIM_PAIRED_FASTQ1,TRIM_PAIRED_FASTQ2=$TRIM_PAIRED_FASTQ2,TRIM_UNPAIRED_FASTQ1=$TRIM_UNPAIRED_FASTQ1,TRIM_UNPAIRED_FASTQ2=$TRIM_UNPAIRED_FASTQ2,FASTQC_OUTPUT_DIR=$FASTQC_OUTPUT_DIR,DONE=$DONE \
      ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
    preprocess=$(echo $preprocess)
    echo "[Q] Preprocessing     : $preprocess"
    DEPEND_PREPROCESS="-W "$dependMod"afterok:"`echo $preprocess`;
  else echo "[D] Preprocessing"
  fi

#ALIGN PAIRED READS (bwa mem)
  NAME=pipeline.MIP.part_1_4.align_bwa_mem_pe.in_fastq
  DONE=.${SAMPLE}-${RUN}.${NAME}.done
  if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
    alignPE=`$QUEUE $ACCOUNT  \
      -l nodes=1:ppn=$THREADS \
      -l mem=${MAX_MEM}g \
      -l walltime=$WALLTIME \
      -N ${SAMPLE}-${RUN}.${NAME} \
      $DEPEND_PREPROCESS \
      -V \
      -v FASTQ1=$TRIM_PAIRED_FASTQ1,FASTQ2=$TRIM_PAIRED_FASTQ2,SAM=$SAM,RUN=$RG,DONE=$DONE \
      ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
    alignPE=$(echo $alignPE)
    echo "[Q] Align PE    : $alignPE [$DEPEND_PREPROCESS]"
    DEPEND_ALN=${DEPEND_ALN}:${alignPE}
  else echo "[D] Align PE"
  fi

}

# ===============================================
# PARSING ARGUMENTS
# ===============================================

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt --name pipeline --alternative --unquoted --options hs:d:t:m:vw:f:c:a:x: --longoptions sample:,dir:,threads:,max_mem:,walltime:,capture_kit:,bed:,bed_flank:,adapters:,mip_design_file:,molecular_tag_length:,tagged_read_number:,r1f1:,run1f1:,r1f2:,run1f2:,r1rg:,run1rg:,r2f1:,run2f1:,r2f2:,run2f2:,r2rg:,run2rg:,r3f1:,run3f1:,r3f2:,run3f2:,r3rg:,run3rg:,r4f1:,run4f1:,r4f2:,run4f2:,r4rg:,run4rg:,r5f1:,run5f1:,r5f2:,run5f2:,r5rg:,run5rg:,r6f1:,run6f1:,r6f2:,run6f2:,r6rg:,run6rg:,r7f1:,run7f1:,r7f2:,run7f2:,r7rg:,run7rg:,r8f1:,run8f1:,r8f2:,run8f2:,r8rg:,run8rg:,r9f1:,run9f1:,r9f2:,run9f2:,r9rg:,run9rg:,extra:,center:,keep_intermediates,discard_unpaired,start_from_scratch,verbose,help -- "$@")
then
    # something went wrong, getopt will put out an error message for us
    echo "Error processing options."
    exit 42
fi

# ===============================================
# LOAD & OVERRIDE EXTRA CONFIG FILE FOR PROJECT
# ===============================================
set -- $options
while [ $# -gt 0 ]
do
    case $1 in
    -x| --extra) 
      EXTRA_CONF="$2" ;
      if [ -f $EXTRA_CONF ]; then
        echo "* LOADING EXTRA CONFIG FILE $EXTRA_CONF";
        . $EXTRA_CONF
      else
        echo "ERROR: invalid EXTRA CONFIG file: $EXTRA_CONF";
        echo "Please check options and try again"; exit 42;
      fi
    esac
    shift
done

# ===============================================
# LOAD ALL OTHER OPTIONS
# ===============================================
set -- $options
while [ $# -gt 0 ]
#echo -e "$1 $2"
do
    case $1 in
    -h| --help) Usage; exit 0;;
    # for options with required arguments, an additional shift is required
    -s| --sample) SAMPLE="$2" ; shift ;;
    -d| --dir) DIR="$2" ; shift ;;
    -t| --threads) THREADS="$2" ; shift ;;
    -m| --max_mem) MAX_MEM="$2" ; shift ;;
    -w| --walltime) WALLTIME="$2" ; shift ;;
    -c| --capture_kit) CAPTURE_KIT="$2" ; shift ;;
    -x| --extra) echo -n "" ; shift ;;
    -a| --adapters) ADAPTERS="$2" ; shift ;;
    -v| --verbose) VERBOSE=1 ;; 
    -f| --mip_design_file) MIP_DESIGN_FILE="$2" ; shift ;;
    --molecular_tag_length) MOLECULAR_TAG_LENGTH="$2" ; shift ;;
    --tagged_read_number) TAGGED_READ_NUMBER="$2" ; shift ;;
    --r1f1| --run1f1) RUN1F1="$2" ; shift ;;
    --r1f2| --run1f2) RUN1F2="$2" ; shift ;;
    --r1rg| --run1rg) RUN1RG="$2" ; shift ;;
    --r2f1| --run2f1) RUN2F1="$2" ; shift ;;
    --r2f2| --run2f2) RUN2F2="$2" ; shift ;;
    --r2rg| --run2rg) RUN2RG="$2" ; shift ;;
    --r3f1| --run3f1) RUN3F1="$2" ; shift ;;
    --r3f2| --run3f2) RUN3F2="$2" ; shift ;;
    --r3rg| --run3rg) RUN3RG="$2" ; shift ;;
    --r4f1| --run4f1) RUN4F1="$2" ; shift ;;
    --r4f2| --run4f2) RUN4F2="$2" ; shift ;;
    --r4rg| --run4rg) RUN4RG="$2" ; shift ;;
    --r5f1| --run5f1) RUN5F1="$2" ; shift ;;
    --r5f2| --run5f2) RUN5F2="$2" ; shift ;;
    --r5rg| --run5rg) RUN5RG="$2" ; shift ;;
    --r6f1| --run6f1) RUN6F1="$2" ; shift ;;
    --r6f2| --run6f2) RUN6F2="$2" ; shift ;;
    --r6rg| --run6rg) RUN6RG="$2" ; shift ;;
    --r7f1| --run7f1) RUN7F1="$2" ; shift ;;
    --r7f2| --run7f2) RUN7F2="$2" ; shift ;;
    --r7rg| --run7rg) RUN7RG="$2" ; shift ;;
    --r8f1| --run8f1) RUN8F1="$2" ; shift ;;
    --r8f2| --run8f2) RUN8F2="$2" ; shift ;;
    --r8rg| --run8rg) RUN8RG="$2" ; shift ;;
    --r9f1| --run9f1) RUN9F1="$2" ; shift ;;
    --r9f2| --run9f2) RUN9F2="$2" ; shift ;;
    --r9rg| --run9rg) RUN9RG="$2" ; shift ;;
    --bed) BED="$2" ; shift ;;
    --bed_flank) BED_FLANK="$2" ; shift ;;
    --center) CENTER="$2" ; shift ;;
    --start_from_scratch) SFS=1 ;;
    --discard_unpaired) DISCARD_UNPAIRED_READS=1 ;; 
    --keep_intermediates) KEEP_INTERMEDIATES=1 ;; 
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 42;;
    (*) break;;
    esac
    shift
done

FOUND_ERROR=0
# ===============================================
# CHECKING VARIABLES
# ===============================================
#check to ensure all mandatory arguments have been entered
if [ -z $SAMPLE ]; then echo "ERROR: missing mandatory option: -s (sample name) must be specified"; FOUND_ERROR=1; fi
if [ -z $DIR ]; then echo "ERROR: missing mandatory option: -d (--dir) must be specified"; FOUND_ERROR=1; fi
if [ ! -d $DIR ]; then echo "ERROR: mandatory option: -d (--dir) is invalid - directory must exist"; FOUND_ERROR=1; fi
if [ -z $CAPTURE_KIT ]; then echo "ERROR: missing mandatory option: -c (capture kit) must be specified"; FOUND_ERROR=1; fi
if [ -z $ADAPTERS ]; then echo "ERROR: missing mandatory option: -a (adapters sequence) must be specified"; FOUND_ERROR=1; 
else
ADAPTERS=$(cd $(dirname $ADAPTERS) && pwd -P)/$(basename $ADAPTERS);
if [ ! -f $ADAPTERS ]; then echo "ERROR: invalid adapters file option: '$ADAPTERS' could not be found"; FOUND_ERROR=1; fi
fi
#f [ -z $FASTQ1 ]; then echo "ERROR: missing mandatory option: --fastq1 must be specified"; FOUND_ERROR=1; fi
#f [ -z $FASTQ2 ]; then echo "ERROR: missing mandatory option: --fastq2 must be specified"; FOUND_ERROR=1; fi
if [ -z $EXTRA_CONF ]; then echo "ERROR: missing mandatory option: --extra must be specified"; FOUND_ERROR=1; fi
if [ ! -s $BED -a ! -s $BED.targets-only.bed -a ! -s $BED.targets-only.non-overlapping.bed ]; then echo "ERROR: specified bedfile or sub-bedfiles [ $BED.bed.targets-only.bed and $BED.bed-targets-only.non-overlapping.bed ] do not exist"; FOUND_ERROR=1; fi
#if [ -z $SITES ]; then echo "ERROR: missing mandatory option: --sites must be specified"; FOUND_ERROR=1; fi
#if (( $FOUND_ERROR )); then echo "Please check options and try again"; exit 42; fi
# ===============================================
#setting path and execution variables
# ===============================================
if [[ -n $ACCOUNT ]]; then ACCOUNT="-A $ACCOUNT"; fi
OLD_PATH=`pwd -P`
[ `echo $DIR | grep -P "^\."` ] && DIR=$OLD_PATH/$DIR; 
DIR=$(cd $DIR 2>/dev/null && pwd -P);
if [ ! -d $DIR ]; then echo "ERROR: wrong mandatory option: -d (run directory $DIR) is invalid"; FOUND_ERROR=1; fi

checkRawFastq() {
  read forward reverse rg <<< $@
  if [[ -n $forward || -n $reverse || -n $rg ]]; then
    if [[ -n $forward && -n $reverse && -n $rg ]]; then
      [ `echo $forward | grep -P "^\."` ] && forward=$OLD_PATH/$forward;
      forward=$(cd $(dirname $forward) && pwd -P)/$(basename $forward);
      if [ ! -f $forward ]; then echo "ERROR: invalid forward file option: '$forward' could not be found"; FOUND_ERROR=1; fi
      [ `echo $reverse | grep -P "^\."` ] && reverse=$OLD_PATH/$reverse;
      reverse=$(cd $(dirname $reverse) && pwd -P)/$(basename $reverse);
      if [ ! -f $reverse ]; then echo "ERROR: invalid reverse file option: '$reverse' could not be found"; FOUND_ERROR=1; fi
    else echo "ERROR: Please define all variables for RUN $rg"; FOUND_ERROR=1; fi
  fi
}

checkRawFastq $RUN1F1 $RUN1F2 $RUN1RG
checkRawFastq $RUN2F1 $RUN2F2 $RUN2RG
checkRawFastq $RUN3F1 $RUN3F2 $RUN3RG
checkRawFastq $RUN4F1 $RUN4F2 $RUN4RG
checkRawFastq $RUN5F1 $RUN5F2 $RUN5RG
checkRawFastq $RUN6F1 $RUN6F2 $RUN6RG
checkRawFastq $RUN7F1 $RUN7F2 $RUN7RG
checkRawFastq $RUN8F1 $RUN8F2 $RUN8RG
checkRawFastq $RUN8F1 $RUN8F2 $RUN8RG

if (( $FOUND_ERROR )); then echo "Please check options and try again"; exit 42; fi

export SAMPLE DIR THREADS MAX_MEM WALLTIME CAPTURE_KIT BED BED_FLANKS RUN1F1 RUN1F2 RUN1RG RUN2F1 RUN2F2 RUN2RG RUN3F1 RUN3F2 RUN3RG RUN4F1 RUN4F2 RUN4RG RUN5F1 RUN5F2 RUN5RG RUN6F1 RUN6F2 RUN6RG RUN7F1 RUN7F2 RUN7RG RUN8F1 RUN8F2 RUN8RG RUN9F1 RUN9F2 RUN9RG VERBOSE SFS KEEP_INTERMEDIATES ADAPTERS CENTER

# ===============================================
# QUEUE THE JOB IN CLUSTER
# ===============================================
cd $DIR; # for job outputs

# ALIGNMENTS
unset MERGE_SORT_SAMS
unset ALIGN_IDS

unset DEPEND_ALN;

alignFastq() {
  read fastq1 fastq2 rg mip_file id <<< $@
  if [[ -n $fastq1 && -n $fastq2 && -n $rg && -n $mip_file ]]; then
    SAM=${SAMPLE}.${rg}.paired.aligned.sam
    MERGE_SORT_SAMS="${MERGE_SORT_SAMS}@INPUT=${SAM}"
    align $fastq1 $fastq2 $mip_file $rg $SAM $id
    if [[ ${DEPEND_ALN} ]]; then ALIGN_IDS=${ALIGN_IDS}${DEPEND_ALN}; fi
  fi
}

alignFastq $RUN1F1 $RUN1F2 $RUN1RG $MIP_DESIGN_FILE 1
alignFastq $RUN2F1 $RUN2F2 $RUN2RG $MIP_DESIGN_FILE 2
alignFastq $RUN3F1 $RUN3F2 $RUN3RG $MIP_DESIGN_FILE 3
alignFastq $RUN4F1 $RUN4F2 $RUN4RG $MIP_DESIGN_FILE 4
alignFastq $RUN5F1 $RUN5F2 $RUN5RG $MIP_DESIGN_FILE 5
alignFastq $RUN6F1 $RUN6F2 $RUN6RG $MIP_DESIGN_FILE 6
alignFastq $RUN7F1 $RUN7F2 $RUN7RG $MIP_DESIGN_FILE 7
alignFastq $RUN8F1 $RUN8F2 $RUN8RG $MIP_DESIGN_FILE 8
alignFastq $RUN9F1 $RUN9F2 $RUN9RG $MIP_DESIGN_FILE 9

echo "*** Post Alignment"

# MERGE SAM FILES
NAME=pipeline.generic.part_2_0.combine_sample_libs.primary_aln
DONE=.${SAMPLE}.${NAME}.done
#touch $DONE
MERGE_SORT_BAM_PRIMARY_ALN=${SAMPLE}.aligned.merged.sorted.bam
MERGE_SORT_BAM_NON_PRIMARY_ALN=${SAMPLE}.aligned.merged.sorted.non_pri_aln.bam
if [[ ! -z $ALIGN_IDS ]]; then DEPEND="-W "$dependMod"afterok"`echo ${ALIGN_IDS}`; else unset DEPEND; fi
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  merge_sample_libs=`$QUEUE $ACCOUNT  \
    -l nodes=1:ppn=$THREADS \
    -l mem=${MAX_MEM}g \
    -l walltime=$WALLTIME \
    -N ${SAMPLE}.${NAME} \
    $DEPEND \
    -V \
    -v MERGE_SORT_SAMS=${MERGE_SORT_SAMS},MERGE_SORT_BAM=${MERGE_SORT_BAM_PRIMARY_ALN},MERGE_SORT_BAM_NON_PRIMARY_ALN=${MERGE_SORT_BAM_NON_PRIMARY_ALN},DONE=${DONE} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  merge_sample_libs=$(echo $merge_sample_libs)
  echo "[Q] Merge Sort  : $merge_sample_libs [$DEPEND]"
  MERGE_DEPEND="-W "$dependMod"afterok:"`echo $merge_sample_libs`;
else 
  echo "[D] Merge Sort"
fi

# GATK EXOME PRE-PROCESSING
NAME=pipeline.MIP.part_2_1.GATK_data_processing.in_bam
DONE=.${SAMPLE}.${NAME}.done
#touch $DONE
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  gatk=`$QUEUE $ACCOUNT  \
    -l nodes=1:ppn=$THREADS \
    -l mem=${MAX_MEM}g \
    -l walltime=$WALLTIME \
    -N ${SAMPLE}.${NAME} \
    $MERGE_DEPEND \
    -V \
    -v ALIGNED_BAM=${MERGE_SORT_BAM_PRIMARY_ALN},USE_INDELS=1,DONE=${DONE} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  gatk=$(echo $gatk)
  echo "[Q] GATK Exome      : $gatk [$MERGE_DEPEND]"
  GATK_EXOME_DEPEND="-W "$dependMod"afterok:"`echo $gatk`;
else echo "[D] GATK Exome"
fi

## Exome STATS AND QC
NAME=pipeline.MIP.part_3_1.Exome_stats_QC
DONE=.${SAMPLE}.${NAME}.done
#touch $DONE
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  Exome_sample_stats_QC=`$QUEUE $ACCOUNT  \
    -l nodes=1:ppn=$THREADS \
    -l mem=${MAX_MEM}g \
    -l walltime=$WALLTIME \
    -N ${SAMPLE}.${NAME}\
    $GATK_EXOME_DEPEND \
    -V \
    -v MIP_DESIGN_FILE=$MIP_DESIGN_FILE,CAPTURE_KIT=${CAPTURE_KIT},BED=$BED,DONE=${DONE} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  Exome_sample_stats_QC=$(echo $Exome_sample_stats_QC)
  echo "[Q] Ex Stats and QC : $Exome_sample_stats_QC [$GATK_EXOME_DEPEND]"
  EXOME_STATS_QC_DEPEND="-W "$dependMod"afterok:"`echo $Exome_sample_stats_QC`;
else echo "[D] Exome Stats and QC"
fi

## GATK UG total recall 
NAME=pipeline.MIP.part_3_5.GATK_UG_total_recall
DONE=.${SAMPLE}.${NAME}.done
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  gatk_UG_total_recall=`$QUEUE $ACCOUNT  \
    -l nodes=1:ppn=$THREADS \
    -l mem=${MAX_MEM}g \
    -l walltime=$WALLTIME \
    -N ${SAMPLE}.${NAME}\
    $GATK_EXOME_DEPEND \
    -V \
    -v DONE=${DONE} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  gatk_UG_total_recall=$(echo $gatk_UG_total_recall)
  echo "[Q] GATK UG tot rec : $gatk_UG_total_recall [$GATK_EXOME_DEPEND]"
else echo "[D] GATK UG total recall"
fi

cd $OLD_PATH
exit 0
