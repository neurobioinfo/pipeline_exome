#!/bin/bash

#hard-coded version number. Don't forget to update!!
VERSION=0.2.0; echo "custom capture pipeline version $VERSION"
#last updated version 2017-09-25

# ===============================================
# default variables values
# ===============================================
unset SAMPLE DIR THREADS MAX_MEM PIPELINE_HOME WALLTIME EXTRA_CONF CAPTURE_KIT RUNF1 RUNF2 RUNLIB RUNRUN QUEUE ACCOUNT VERBOSE CENTER SFS KEEP_INTERMEDIATES REF FILTER CROP_COMMAND CROP_SIZE MATCH_SIZE

#Assuming script is in root of PIPELINE_HOME
export PIPELINE_HOME=$(cd $(dirname $0) && pwd -P)
export SFS=0
export KEEP_INTERMEDIATES=0
export DISCARD_UNPAIRED_READS=0
export FILTER=0
export EXPECTED_DONE_FILES='expected.done.files.txt'
export NO_RECAL=0
export CROP_SIZE=1000
export CROP_COMMAND="CROP:1000"
export MATCH_SIZE=1

rm -f $EXPECTED_DONE_FILES

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
          "\t\t-s  (--sample)           = sample name\n" \
          "\t\t-d  (--dir)              = run directory (where all the outputs will be printed)\n" \
          "\t\t-x  (--extra)            = use extra config file to (re)define variables [CURRENT \"$EXTRA_CONF\"]\n" \
          "\t\t-rf1                     = raw fastq sequence forward\n" \
          "\t\t-rf2                     = raw fastq sequence reverse\n" \
          "\t\t-rl                      = raw fastq sequence library name\n" \
          "\t\t-rr                      = raw fastq sequence run name\n"
	echo -e "\toptional arguments:\n" \
          "\t\t-h  (--help)             = get the program options and exit\n" \
          "\t\t--bed                    = specify the bed file for QC metrics\n" \
          "\t\t--ref                    = specify the reference file to use\n" \
          "\t\t--filter                 = use filtering to keep only highly specific reads\n" \
          "\t\t--match_size             = when filtering, keep only reads with at least that match size [CURRENT \"$MATCH_SIZE\"]\n" \
          "\t\t--norecal                = use filtering to keep only highly specific reads\n" \
          "\t\t--crop_size              = max size of the reads to use, by cropping raw reads [CURRENT \"$CROP_SIZE\"]\n" \
          "\t\t-a  (--adapters)         = fasta sequence file of adapters\n" \
          "\t\t-m  (--max_mem)          = max memory on a node (in gigs) to use [CURRENT \"$MAX_MEM\"]\n" \
          "\t\t-t  (--threads)          = number of threads to use [CURRENT \"$THREADS\"]\n" \
          "\t\t-w  (--walltime)         = wall time for the scheduler (format HH:MM:SS) [CURRENT \"$WALLTIME\"]\n" \
          "\t\t-v  (--verbose)          = set verbosity level [CURRENT \"$VERBOSE\"]\n" \
          "\t\t--account                = set account on cluster (to be used with qsub with -A [CURRENT \"$ACCOUNT\"]\n" \
          "\t\t--center                 = specify sequencing center [CURRENT \"$CENTER\"]\n" \
          "\t\t--start_from_scratch     = flag to overwrite existing files [CURRENT \"$SFS\"]\n" \
          "\t\t--discard_unpaired       = flag to discard unpaired reads after pre processing [CURRENT \"$DISCARD_UNPAIRED_READS\"]\n" \
          "\t\t--keep_intermediates     = flag to keep intermediates files [CURRENT \"$KEEP_INTERMEDIATES\"]" 
        echo
}

# ===============================================
# alignment function (bwa mem)
# ===============================================
align() {
  local F1=$1
  local F2=$2
  local LIB=$3
  local RUN=$4
  local CROP_SIZE=$5
  local MATCH_SIZE=$6
  local SAM_V37=$7
  local SAM_V37_UNPAIRED=$8

  unset DEPEND_V37_ALN;

  export CROP_COMMAND="CROP:$CROP_SIZE" ;

#DETERMINE BASE
  export ENCODING=$(${PIPELINE_HOME}/soft/bin/determine_base.perl $F1);
  if [[ "$ENCODING" == "Cannot determine base" ]]; then echo "ERROR: Cannot determine base for $F1"; exit 42;
  else echo "*** ENCODING bases set to $ENCODING"; fi

  echo -e "*** ALIGNMENT v37\n*** R1:$F1\n*** R2:$F2\n*** LIB:$LIB\n*** RUN:$RUN"
  if [[ -n $BED ]]; then echo -e "*** BED: $BED"; fi
  if [[ -n $REF ]]; then echo -e "*** REF: $REF"; fi
  if [[ -n $CROP_SIZE ]]; then echo -e "*** CROP_SIZE: $CROP_SIZE"; fi
  if [[ -n $MATCH_SIZE ]]; then echo -e "*** MATCH_SIZE: $MATCH_SIZE"; fi

#PRE PREPROCESSING
  FASTQC_OUTPUT_DIR=fastqc.${LIB}-${RUN}
  export TRIM_PAIRED_FASTQ1=${LIB}-${RUN}.r1.fastq.trimmomaticPE.paired.gz
  export TRIM_PAIRED_FASTQ2=${LIB}-${RUN}.r2.fastq.trimmomaticPE.paired.gz
  export TRIM_UNPAIRED_FASTQ1=${LIB}-${RUN}.r1.fastq.trimmomaticPE.unpaired.gz
  export TRIM_UNPAIRED_FASTQ2=${LIB}-${RUN}.r2.fastq.trimmomaticPE.unpaired.gz
  export CROP_COMMAND="CROP:$CROP_SIZE" ;

  NAME=pipeline.generic.part_0_1.trimmomatic
  export DONE=.${LIB}-${RUN}.${NAME}.done
  echo $DONE >> $EXPECTED_DONE_FILES
  #touch $DONE
  if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
    trimm=`$QUEUE $ACCOUNT  \
      -l nodes=1:ppn=${THREADS} \
      -l mem=${MAX_MEM}g \
      -l walltime=$WALLTIME \
      $DEPEND_SAMTOFASTQ \
      -N ${LIB}-${RUN}.${NAME} \
      -V \
      -v FASTQ1=$F1,FASTQ2=$F2,TRIM_PAIRED_FASTQ1=$TRIM_PAIRED_FASTQ1,TRIM_PAIRED_FASTQ2=$TRIM_PAIRED_FASTQ2,TRIM_UNPAIRED_FASTQ1=$TRIM_UNPAIRED_FASTQ1,TRIM_UNPAIRED_FASTQ2=$TRIM_UNPAIRED_FASTQ2,FASTQC_OUTPUT_DIR=$FASTQC_OUTPUT_DIR,DONE=$DONE,CROP_COMMAND=$CROP_COMMAND \
      ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
    trimm=$(echo $trimm)
    echo "[Q] Trimmomatic     : $trimm [$DEPEND_SAMTOFASTQ]"
    DEPEND_TRIM="-W "$dependMod"afterok:"`echo $trimm`;
  else echo "[D] Trimmomatic"
  fi

#ALIGN PAIRED READS (bwa mem)
  if [[ -n $REF ]]; then GENOME=$REF; else GENOME=v37; fi
  NAME=pipeline.generic.part_1_4.align_bwa_mem_pe.in_fastq
  DONE=.${LIB}-${RUN}.${NAME}.v37.done
  echo $DONE >> $EXPECTED_DONE_FILES
  if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
    v37_alignPE=`$QUEUE $ACCOUNT  \
      -l nodes=1:ppn=$THREADS \
      -l mem=${MAX_MEM}g \
      -l walltime=$WALLTIME \
      -N ${LIB}-${RUN}.${NAME}.v37 \
      $DEPEND_TRIM \
      -V \
      -v FASTQ1=$TRIM_PAIRED_FASTQ1,FASTQ2=$TRIM_PAIRED_FASTQ2,SAM=$SAM_V37,LIB=$LIB,RUN=$RUN,GENOME=$GENOME,DONE=$DONE \
      ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
    v37_alignPE=$(echo $v37_alignPE)
    echo "[Q] Align v37 PE    : $v37_alignPE [$DEPEND_TRIM]"
    DEPEND_V37_ALN=${DEPEND_V37_ALN}:${v37_alignPE}
  else echo "[D] Align v37 PE"
  fi

if [[ $DISCARD_UNPAIRED_READS -eq 0 ]]; then
  #ALIGN UNPAIRED (bwa mem)
  NAME=pipeline.generic.part_1_5.align_bwa_mem_se.in_fastq
  DONE=.${LIB}-${RUN}.${NAME}.v37.done
  echo $DONE >> $EXPECTED_DONE_FILES
  #touch $DONE
  if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
    v37_alignSE=`$QUEUE $ACCOUNT  \
      -l nodes=1:ppn=$THREADS \
      -l mem=${MAX_MEM}g \
      -l walltime=$WALLTIME \
      -N ${LIB}-${RUN}.${NAME}.v37 \
      $DEPEND_TRIM \
      -V \
      -v INPUT="custom",FASTQ1=$TRIM_UNPAIRED_FASTQ1,FASTQ2=$TRIM_UNPAIRED_FASTQ2,SAM=$SAM_V37_UNPAIRED,LIB=$LIB,RUN=$RUN,GENOME=v37,DONE=$DONE \
      ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
    v37_alignSE=$(echo $v37_alignSE)
    echo "[Q] Align v37 SE    : $v37_alignSE [$DEPEND_TRIM]"
    DEPEND_V37_ALN=${DEPEND_V37_ALN}:${v37_alignSE}
  else echo "[D] Align v37 SE"
  fi
fi


if [[ $FILTER -eq 1 ]]; then
  if [[ ! -z $DEPEND_V37_ALN ]]; then DEPEND="-W "$dependMod"afterok"`echo ${DEPEND_V37_ALN}`; else unset DEPEND; fi
  NAME=pipeline.generic.part_1_8.filter_align
  DONE=.Filter-${LIB}-${RUN}.${NAME}.v37.done
  echo $DONE >> $EXPECTED_DONE_FILES
  if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  filter_reads=`$QUEUE $ACCOUNT  \
      -l nodes=1:ppn=$THREADS \
      -l mem=${MAX_MEM}g \
      -l walltime=$WALLTIME \
      -N Filter-${LIB}-${RUN}.${NAME}.v37 \
      $DEPEND \
      -V \
      -v SAM=$SAM_V37,MIN_MATCH_SIZE=$MATCH_SIZE \
      ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  filter_reads=$(echo $filter_reads)
  echo "[Q] Filter PE       : $filter_reads [$DEPEND_V37_ALN]"
  DEPEND_V37_ALN=${DEPEND_V37_ALN}:${filter_reads}
  else echo "[D] Filter PE"
  fi
fi
}

# ===============================================
# PARSING ARGUMENTS
# ===============================================

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt --name pipeline --alternative --unquoted --options hs:d:t:m:vw:c:a:x: --longoptions sample:,dir:,threads:,max_mem:,walltime:,capture_kit:,mode:,bed:,ref:,filter,norecal,sites:,crop_size:,match_size:,adapters:,rf1:,rf2:,rr:,rl:,extra:,center:,account:,keep_intermediates,discard_unpaired,start_from_scratch,verbose,help -- "$@")
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
    -x| --extra) echo -n "" ; shift ;;
    -a| --adapters) ADAPTERS="$2" ; shift ;;
    -v| --verbose) VERBOSE=1 ;; 
    --account) ACCOUNT="$2" ; shift ;;
    --rf1) RUNF1="$2" ; shift ;;
    --rf2) RUNF2="$2" ; shift ;;
    --rr) RUNRUN="$2" ; shift ;;
    --rl) RUNLIB="$2" ; shift ;;
    --bed) BED="$2" ; shift ;;
    --ref) REF="$2" ; shift ;;
    --filter) FILTER=1;;
    --crop_size) CROP_SIZE="$2" ; shift ;;
    --match_size) MATCH_SIZE="$2" ; shift ;;
    --norecal) NO_RECAL=1;;
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
if [[ -z $RUNF1 || -z $RUNF2 ]]; then echo "ERROR: must supply at least 1 set of raw forward and reverse fastq files if using \"--mode fastq\""; FOUND_ERROR=1; fi
if [ -z $ADAPTERS ]; then 
  echo "ERROR: missing mandatory option: -a (adapters sequence) must be specified"; FOUND_ERROR=1; 
else
  ADAPTERS=$(cd $(dirname $ADAPTERS) && pwd -P)/$(basename $ADAPTERS);
  if [ ! -f $ADAPTERS ]; then echo "ERROR: invalid adapters file option: '$ADAPTERS' could not be found"; FOUND_ERROR=1; fi
fi
if [ -z $EXTRA_CONF ]; then echo "ERROR: missing mandatory option: --extra must be specified"; FOUND_ERROR=1; fi
# ===============================================
#setting path and execution variables
# ===============================================
if [[ -n $ACCOUNT ]]; then ACCOUNT="-A $ACCOUNT"; fi
OLD_PATH=`pwd -P`
[ `echo $DIR | grep -P "^\."` ] && DIR=$OLD_PATH/$DIR; 
DIR=$(cd $DIR 2>/dev/null && pwd -P);
if [ ! -d $DIR ]; then echo "ERROR: wrong mandatory option: -d (run directory $DIR) is invalid"; FOUND_ERROR=1; fi

if [[ -n $BED ]]; then 
  BED=$(cd $(dirname $BED) && pwd -P)/$(basename $BED);
  if [ ! -f $BED ]; then echo "ERROR: invalid bed file option: '$BED' could not be found"; FOUND_ERROR=1; fi
fi

if [[ -n $REF ]]; then
  REF=$(cd $(dirname $REF) && pwd -P)/$(basename $REF);
  if [ ! -f $REF ]; then echo "ERROR: invalid ref file option: '$REF' could not be found"; FOUND_ERROR=1; fi
fi

checkRawFastq() {
  read forward reverse run lib <<< $@
  if [[ -n $forward || -n $reverse || -n $run || -n $lib ]]; then
    if [[ -n $forward && -n $reverse && -n $run && -n $lib ]]; then
      [ `echo $forward | grep -P "^\."` ] && forward=$OLD_PATH/$forward;
      forward=$(cd $(dirname $forward) && pwd -P)/$(basename $forward);
      if [ ! -f $forward ]; then echo "ERROR: invalid forward file option: '$forward' could not be found"; FOUND_ERROR=1; fi
      [ `echo $reverse | grep -P "^\."` ] && reverse=$OLD_PATH/$reverse;
      reverse=$(cd $(dirname $reverse) && pwd -P)/$(basename $reverse);
      if [ ! -f $reverse ]; then echo "ERROR: invalid reverse file option: '$reverse' could not be found"; FOUND_ERROR=1; fi
    else echo "ERROR: Please define all variables for RUN $run"; FOUND_ERROR=1; fi
  fi
}

checkRawFastq $RUNF1 $RUNF2 $RUNRUN $RUNLIB

if (( $FOUND_ERROR )); then echo "Please check options and try again"; exit 42; fi

export SAMPLE DIR THREADS MAX_MEM WALLTIME CAPTURE_KIT BED BED_FLANKS RUNF1 RUNF2 RUNRUN RUNLIB VERBOSE SFS KEEP_INTERMEDIATES ADAPTERS CENTER

# ===============================================
# QUEUE THE JOB IN CLUSTER
# ===============================================
cd $DIR; # for job outputs

# ALIGNMENTS
unset MERGE_V37_SORT_SAMS
unset ALIGN_V37_IDS

alignFastq() {
  read forward reverse run lib crop_size match_size<<< $@
  if [[ -n $forward && -n $reverse && -n $run && -n $lib ]]; then
    V37SAM=${lib}-${run}.paired.aligned.v37.sam
    MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${V37SAM}"
    if [[ $DISCARD_UNPAIRED_READS == 0 ]]; then 
      V37SAMUNPAIRED=${lib}-${run}.unpaired.aligned.v37.sam;
      MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${V37SAMUNPAIRED}"
    fi
    align $forward $reverse $lib $run $crop_size $match_size $V37SAM $V37SAMUNPAIRED;
    if [[ ${DEPEND_V37_ALN} ]]; then ALIGN_V37_IDS=${ALIGN_V37_IDS}${DEPEND_V37_ALN}; fi
  fi
}

unset DEPEND;
alignFastq $RUNF1 $RUNF2 $RUNRUN $RUNLIB $CROP_SIZE $MATCH_SIZE

echo "*** Post Alignment"

# MERGE SAM FILES
NAME=pipeline.generic.part_2_0.combine_sample_libs.primary_aln
DONE=.${SAMPLE}.${NAME}.v37.done
echo $DONE >> $EXPECTED_DONE_FILES
MERGE_V37_SORT_BAM_PRIMARY_ALN=${SAMPLE}.aligned.v37.merged.sorted.bam
MERGE_V37_SORT_BAM_NON_PRIMARY_ALN=${SAMPLE}.aligned.v37.merged.sorted.non_pri_aln.bam
if [[ ! -z $ALIGN_V37_IDS ]]; then DEPEND="-W "$dependMod"afterok"`echo ${ALIGN_V37_IDS}`; else unset DEPEND; fi
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  v37_merge_sample_libs=`$QUEUE $ACCOUNT \
    -l nodes=1:ppn=$THREADS \
    -l mem=${MAX_MEM}g \
    -l walltime=$WALLTIME \
    -N ${SAMPLE}.${NAME}.v37 \
    $DEPEND \
    -V \
    -v MERGE_SORT_SAMS=${MERGE_V37_SORT_SAMS},MERGE_SORT_BAM=${MERGE_V37_SORT_BAM_PRIMARY_ALN},MERGE_SORT_BAM_NON_PRIMARY_ALN=${MERGE_V37_SORT_BAM_NON_PRIMARY_ALN},DONE=${DONE} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  v37_merge_sample_libs=$(echo $v37_merge_sample_libs)
  echo "[Q] Merge Sort v37  : $v37_merge_sample_libs [$DEPEND]"
  V37_MERGE_DEPEND="-W "$dependMod"afterok:"`echo $v37_merge_sample_libs`;
else 
  echo "[D] Merge Sort v37"
fi

#Set the Scala no recal thing here
if [[ $NO_RECAL -eq 1 ]]; then 
  SCALA_FILE="NO_RECAL";#$CUSTOM_POST_PROCESSING_SCALA_NO_RECAL;
else
  SCALA_FILE="EXOME";#$EXOME_POST_PROCESSING_SCALA;
fi

# GATK EXOME PRE-PROCESSING
NAME=pipeline.generic.part_2_1.GATK_data_processing.in_bam
DONE=.${SAMPLE}.${NAME}.v37.done
echo $DONE >> $EXPECTED_DONE_FILES
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  gatk_v37=`$QUEUE $ACCOUNT  \
    -l nodes=1:ppn=$THREADS \
    -l mem=${MAX_MEM}g \
    -l walltime=$WALLTIME \
    -N ${SAMPLE}.${NAME}.v37 \
    $V37_MERGE_DEPEND \
    -V \
    -v ALIGNED_BAM=${MERGE_V37_SORT_BAM_PRIMARY_ALN},USE_INDELS=1,PREFIX="Custom",DONE=${DONE},SCALA=${SCALA_FILE} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  gatk_v37=$(echo $gatk_v37)
  echo "[Q] GATK custom     : $gatk_v37 [$V37_MERGE_DEPEND]"
  GATK_EXOME_DEPEND="-W "$dependMod"afterok:"`echo $gatk_v37`;
else echo "[D] GATK custom"
fi

## STATS AND QC
NAME=pipeline.generic.part_3_0.Capture_stats_QC
DONE=.${SAMPLE}.${NAME}.done
echo $DONE >> $EXPECTED_DONE_FILES
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  Exome_sample_stats_QC=`$QUEUE $ACCOUNT  \
    -l nodes=1:ppn=$THREADS \
    -l mem=${MAX_MEM}g \
    -l walltime=$WALLTIME \
    -N ${SAMPLE}.${NAME}\
    $GATK_EXOME_DEPEND \
    -V \
    -v BED=${BED},DONE=${DONE} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  Exome_sample_stats_QC=$(echo $Exome_sample_stats_QC)
  echo "[Q] Ex Stats and QC : $Exome_sample_stats_QC [$GATK_EXOME_DEPEND]"
  EXOME_STATS_QC_DEPEND="-W "$dependMod"afterok:"`echo $Exome_sample_stats_QC`;
else echo "[D] Exome Stats and QC"
fi

## GATK UG total recall 
NAME=pipeline.generic.part_3_5.GATK_UG_total_recall
DONE=.${SAMPLE}.${NAME}.done
echo $DONE >> $EXPECTED_DONE_FILES
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  gatk_UG_total_recall=`$QUEUE $ACCOUNT  \
    -l nodes=1:ppn=$THREADS \
    -l mem=${MAX_MEM}g \
    -l walltime=24:0:0 \
    -N ${SAMPLE}.${NAME}\
    $GATK_EXOME_DEPEND \
    -V \
    -v BEDFILE=${BED},DONE=${DONE} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  gatk_UG_total_recall=$(echo $gatk_UG_total_recall)
  echo "[Q] GATK UG         : $gatk_UG_total_recall [$GATK_EXOME_DEPEND]"
else echo "[D] GATK UG "
fi

## GATK HC total recall 
NAME=pipeline.generic.part_3_6.GATK_HC
DONE=.${SAMPLE}.${NAME}.done
echo $DONE >> $EXPECTED_DONE_FILES
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  gatk_HC=`$QUEUE $ACCOUNT  \
    -l nodes=1:ppn=$THREADS \
    -l mem=${MAX_MEM}g \
    -l walltime=48:0:0 \
    -N ${SAMPLE}.${NAME}\
    $GATK_EXOME_DEPEND \
    -V \
    -v BEDFILE=${BED},DONE=${DONE} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  gatk_HC=$(echo $gatk_HC)
  echo "[Q] GATK HC         : $gatk_HC [$GATK_EXOME_DEPEND]"
else echo "[D] GATK HC"
fi

cd $OLD_PATH
exit 0
