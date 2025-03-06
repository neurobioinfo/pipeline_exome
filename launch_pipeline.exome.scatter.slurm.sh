#!/bin/bash

#hard-coded version number. Don't forget to update!!
VERSION=1.2.0; echo "exome analysis pipeline version $VERSION"
#last updated version 2017-09-13

# ===============================================
# default variables values
# ===============================================
unset SAMPLE DIR THREADS PIPELINE_HOME WALLTIME MAX_MEM EXTRA_CONF CAPTURE_KIT RUN1F1 RUN1F2 RUN1LIB RUN1RUN RUN2F1 RUN2F2 RUN2LIB RUN2RUN RUN3F1 RUN3F2 RUN3LIB RUN3RUN RUN4F1 RUN4F2 RUN4LIB RUN4RUN RUN5F1 RUN5F2 RUN5LIB RUN5RUN RUN6F1 RUN6F2 RUN6LIB RUN6RUN RUN7F1 RUN7F2 RUN7LIB RUN7RUN RUN8F1 RUN8F2 RUN8LIB RUN8RUN RUN9F1 RUN9F2 RUN9LIB RUN9RUN QUEUE ACCOUNT VERBOSE CENTER SFS KEEP_INTERMEDIATES 

#Assuming script is in root of PIPELINE_HOME
export PIPELINE_HOME=$(cd $(dirname $0) && pwd -P)
export PATH=$PIPELINE_HOME/soft/bin:$PATH 
export SFS=0
export KEEP_INTERMEDIATES=0
export DISCARD_UNPAIRED_READS=0
export NOMT=0
export EXPECTED_DONE_FILES='expected.done.files.txt'

rm -f $EXPECTED_DONE_FILES

#set queue-specific values, depending on system check
QUEUE="sbatch"                                                                         # default job scheduler: slurm

# create function to handle error messages
# ===============================================
Usage() {
        echo
        echo -e "Usage:\t$0 [arguments]"
        echo -e "\tmandatory arguments:\n" \
          "\t\t-s  (--sample)           = sample name\n" \
          "\t\t-d  (--dir)              = run directory (where all the outputs will be printed)\n" \
          "\t\t-x  (--extra)            = use extra config file to (re)define variables [CURRENT \"$EXTRA_CONF\"]\n" \
          "\t\t--mode                   = raw input mode (bam or fastq)\n" \
          "\t     NOTE: must specify either bamfile(s) or fastq files for each read group with the following arguments\n" \
          "\t\t-r[1-9]b  (--run[1-9]b)  = raw bam   sequence run [1-9]\n" \
          "\t\t-r[1-9]f1 (--run[1-9]f1) = raw fastq sequence run [1-9] forward\n" \
          "\t\t-r[1-9]f2 (--run[1-9]f2) = raw fastq sequence run [1-9] reverse\n" \
          "\t\t-r[1-9]l (--run[1-9]lib) = raw fastq sequence run [1-9] library name\n" \
          "\t\t-r[1-9]r (--run[1-9]run) = raw fastq sequence run [1-9] name\n"
	echo -e "\toptional arguments:\n" \
          "\t\t-h  (--help)             = get the program options and exit\n" \
          "\t\t-c  (--capture_kit)      = capture kit used for exome targeting (SS38,SS50,SSV4,SSV5,SSV5C,NGv3,NGvCRome,TruSeq,Nextera10,Nextera12 or CCDS)\n" \
          "\t\t--bed                    = specify the bed file for QC metrics\n" \
          "\t\t--bed_flank              = specify the bed file with flanking for QC metrics\n" \
          "\t\t--sites                  = specify the sites to use for variant calling\n" \
          "\t\t-a  (--adapters)         = fasta sequence file of adapters\n" \
          "\t\t-v  (--verbose)          = set verbosity level [CURRENT \"$VERBOSE\"]\n" \
          "\t\t--account                = set account on cluster (to be used with qsub with -A [CURRENT \"$ACCOUNT\"]\n" \
          "\t\t--center                 = specify sequencing center [CURRENT \"$CENTER\"]\n" \
          "\t\t--nomt                   = do not run the MT alignment pipeline [CURRENT \"$NOMT\"]\n" \
          "\t\t--start_from_scratch     = flag to overwrite existing files [CURRENT \"$SFS\"]\n" \
          "\t\t--discard_unpaired       = flag to discard unpaired reads after pre processing [CURRENT \"$DISCARD_UNPAIRED_READS\"]\n" \
          "\t\t--keep_intermediates     = flag to keep intermediates files [CURRENT \"$KEEP_INTERMEDIATES\"]" 
        echo
}

# ===============================================
# alignment function (bwa mem)
# ===============================================
align() {
  local input=$1
  local F1=$2
  local F2=$3
  local LIB=$4
  local RUN=$5
  local SAM_V37=$6
  local SAM_MT=$7
  local SAM_V37_UNPAIRED=$8
  local NOMT=$9

  unset DEPEND_V37_ALN;

#DETERMINE BASE
  if [[ $MODE == "fastq" ]]; then 
    export ENCODING=$(${PIPELINE_HOME}/soft/bin/determine_base.perl $F1);
  else 
    export ENCODING=$(samtools view $RUN1BAM|${PIPELINE_HOME}/soft/bin/determine_base.perl -);
  fi
  if [[ "$ENCODING" == "Cannot determine base" ]]; then echo "ERROR: Cannot determine base for $F1"; exit 42;
  else echo "*** ENCODING bases set to $ENCODING"; fi

  echo -e "*** ALIGNMENT v37\n*** R1:$F1\n*** R2:$F2\n*** LIB:$LIB\n*** RUN:$RUN"
#GENERATE FASTQ (picard samToFastq)
  NAME=pipeline.generic.part_0_0.samToFastq
  if [[ $MODE == "bam" ]]; then
    export DONE=.${input}-${LIB}-${RUN}.${NAME}.done
    export BAM=$(echo $F1|sed 's/_read1.fastq.gz//g')
    echo $DONE >> $EXPECTED_DONE_FILES
    if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
      samToFastq="$QUEUE $ACCOUNT  \
        --ntasks-per-node=${THREADS} \
	    --mem-per-cpu=4g \
        --time=${WALLTIME} \
        --job-name ${input}-${LIB}-${RUN}.${NAME} \
        --output %x.o%j \
        ${PIPELINE_HOME}/data/qsub/${NAME}.qsub"
      samToFastq=$($samToFastq | grep -oP "\d+")
      echo "[Q] SamToFastq      : $samToFastq"
      DEPEND_SAMTOFASTQ="--dependency=afterok:$samToFastq"
    else echo "[D] SamToFastq"
    fi
  fi

#PRE PREPROCESSING
  export FASTQC_OUTPUT_DIR=fastqc.${input}-${LIB}-${RUN}
  export TRIM_PAIRED_FASTQ1=${input}-${LIB}-${RUN}.r1.fastq.trimmomaticPE.paired.gz
  export TRIM_PAIRED_FASTQ2=${input}-${LIB}-${RUN}.r2.fastq.trimmomaticPE.paired.gz
  export TRIM_UNPAIRED_FASTQ1=${input}-${LIB}-${RUN}.r1.fastq.trimmomaticPE.unpaired.gz
  export TRIM_UNPAIRED_FASTQ2=${input}-${LIB}-${RUN}.r2.fastq.trimmomaticPE.unpaired.gz

  NAME=pipeline.generic.part_0_1.trimmomatic
  STEP=0_1
  export THREADS=${THREADS_ARRAY[$STEP]}
  export MAX_MEM=$[$THREADS*4]
  export WALLTIME=${WALLTIME_ARRAY[$STEP]}
  export DONE=.${input}-${LIB}-${RUN}.${NAME}.done
  echo $DONE >> $EXPECTED_DONE_FILES
  if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
    export FASTQ1=$F1
    export FASTQ2=$F2
    trimm="$QUEUE $ACCOUNT  \
      --ntasks-per-node=${THREADS} \
      --mem-per-cpu=4g \
      --time=${WALLTIME} \
      $DEPEND_SAMTOFASTQ \
      --job-name ${input}-${LIB}-${RUN}.${NAME} \
      --output %x.o%j \
      ${PIPELINE_HOME}/data/qsub/${NAME}.qsub"
    trimm=$($trimm | grep -oP "\d+")
    echo "[Q] Trimmomatic     : $trimm [$DEPEND_SAMTOFASTQ]"
    DEPEND_TRIM="--dependency=afterok:$trimm";
  else echo "[D] Trimmomatic"
  fi

#ALIGN PAIRED READS (bwa mem)
  NAME=pipeline.generic.part_1_4.align_bwa_mem_pe.in_fastq
  STEP=1_4g
  export THREADS=${THREADS_ARRAY[$STEP]}
  export MAX_MEM=$[$THREADS*4]
  export WALLTIME=${WALLTIME_ARRAY[$STEP]}
  export DONE=.${input}-${LIB}-${RUN}.${NAME}.v37.done
  echo $DONE >> $EXPECTED_DONE_FILES
  if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
    export FASTQ1=$TRIM_PAIRED_FASTQ1
    export FASTQ2=$TRIM_PAIRED_FASTQ2
    export SAM=$SAM_V37
    export LIB=$LIB
    export RUN=$RUN
    export GENOME=v37
    v37_alignPE="$QUEUE $ACCOUNT  \
      --ntasks-per-node=${THREADS} \
      --mem-per-cpu=4g \
      --time=${WALLTIME} \
      --job-name ${input}-${LIB}-${RUN}.${NAME}.v37 \
      --output %x.o%j \
      $DEPEND_TRIM \
      ${PIPELINE_HOME}/data/qsub/${NAME}.qsub"
    v37_alignPE=$($v37_alignPE | grep -oP "\d+")
    echo "[Q] Align v37 PE    : $v37_alignPE [$DEPEND_TRIM]"
    DEPEND_V37_ALN=${DEPEND_V37_ALN}:${v37_alignPE}
  else echo "[D] Align v37 PE"
  fi

if [[ $DISCARD_UNPAIRED_READS -eq 0 ]]; then
  #ALIGN UNPAIRED (bwa mem)
  NAME=pipeline.generic.part_1_5.align_bwa_mem_se.in_fastq
  STEP=1_5g
  export THREADS=${THREADS_ARRAY[$STEP]}
  export MAX_MEM=$[$THREADS*4]
  export WALLTIME=${WALLTIME_ARRAY[$STEP]}
  export DONE=.${input}-${LIB}-${RUN}.${NAME}.v37.done
  echo $DONE >> $EXPECTED_DONE_FILES
  #touch $DONE
  if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
    export FASTQ1=$TRIM_UNPAIRED_FASTQ1
    export FASTQ2=$TRIM_UNPAIRED_FASTQ2
    export SAM=$SAM_V37_UNPAIRED
    export LIB=$LIB
    export RUN=$RUN
    export INPUT=$input
    export GENOME=v37
    v37_alignSE="$QUEUE $ACCOUNT \
      --ntasks-per-node=${THREADS} \
      --mem-per-cpu=4g \
      --time=${WALLTIME} \
      --job-name ${input}-${LIB}-${RUN}.${NAME}.v37 \
      --output %x.o%j \
      $DEPEND_TRIM \
      ${PIPELINE_HOME}/data/qsub/${NAME}.qsub"
    v37_alignSE=$($v37_alignSE | grep -oP "\d+")
    echo "[Q] Align v37 SE    : $v37_alignSE [$DEPEND_TRIM]"
    DEPEND_V37_ALN=${DEPEND_V37_ALN}:${v37_alignSE}
  else echo "[D] Align v37 SE"
  fi
fi

# MTCHONDRIA
echo -e "*** ALIGNMENT MT\n*** R1:$F1\n*** R2:$F2\n*** LIB:$LIB\n*** RUN:$RUN"
#PARSING ALIGNMENT (bwa sampe)
  NAME=pipeline.generic.part_1_4.align_bwa_mem_pe.in_fastq
  STEP=1_4m
  export THREADS=${THREADS_ARRAY[$STEP]}
  export MAX_MEM=$[$THREADS*4]
  export WALLTIME=${WALLTIME_ARRAY[$STEP]}
  DONE=.${input}-${LIB}-${RUN}.${NAME}.MT.done
  echo $DONE >> $EXPECTED_DONE_FILES
  if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
    export FASTQ1=$TRIM_PAIRED_FASTQ1
    export FASTQ2=$TRIM_PAIRED_FASTQ2
    export SAM=$SAM_MT
    export LIB=$LIB
    export RUN=$RUN
    export GENOME=MT
    MT_alignPE="$QUEUE $ACCOUNT  \
      --ntasks-per-node=${THREADS} \
      --mem-per-cpu=4g \
      --time=${WALLTIME} \
      --job-name ${input}-${LIB}-${RUN}.${NAME}.MT \
      --output %x.o%j \
      $DEPEND_TRIM \
      ${PIPELINE_HOME}/data/qsub/${NAME}.qsub"
    MT_alignPE=$($MT_alignPE | grep -oP "\d+")
    echo "[Q] Align MT PE     : $MT_alignPE [$DEPEND_TRIM]"
    DEPEND_MT_ALN=:${MT_alignPE}
  else echo "[D] Align MT PE"
  fi
}

# ===============================================
# PARSING ARGUMENTS
# ===============================================

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt --name pipeline --alternative --unquoted --options hs:d:vc:a:x: --longoptions sample:,dir:,capture_kit:,mode:,bed:,bed_flank:,sites:,adapters:,r1b:,run1bam:,r2b:,run2bam:,r3b:,run3bam:,r4b:,run4bam:,r5b:,run5bam:,r6b:,run6bam:,r7b:,run7bam:,r8b:,run8bam:,r9b:,run9bam:,r1f1:,run1f1:,r1f2:,run1f2:,r1r:,run1run:,r1l:,run1lib:,r2f1:,run2f1:,r2f2:,run2f2:,r2r:,run2run:,r2l:,run2lib:,r3f1:,run3f1:,r3f2:,run3f2:,r3r:,run3run:,r3l:,run3lib:,r4f1:,run4f1:,r4f2:,run4f2:,r4r:,run4run:,r4l:,run4lib:,r5f1:,run5f1:,r5f2:,run5f2:,r5r:,run5run:,r5l:,run5lib:,r6f1:,run6f1:,r6f2:,run6f2:,r6r:,run6run:,r6l:,run6lib:,r7f1:,run7f1:,r7f2:,run7f2:,r7r:,run7run:,r7l:,run7lib:,r8f1:,run8f1:,r8f2:,run8f2:,r8r:,run8run:,r8l:,run8lib:,r9f1:,run9f1:,r9f2:,run9f2:,r9r:,run9run:,r9l:,run9lib:,extra:,center:,account:,keep_intermediates,discard_unpaired,start_from_scratch,verbose,nomt,help -- "$@")
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
    -c| --capture_kit) CAPTURE_KIT="$2" ; shift ;;
    --sites) CALLING_BED="$2" ; shift ;; 
    -x| --extra) echo -n "" ; shift ;;
    -a| --adapters) ADAPTERS="$2" ; shift ;;
    -v| --verbose) VERBOSE=1 ;; 
    --account) ACCOUNT="$2" ; shift ;;
    --mode) MODE="$2"; shift ;;
    --r1b| --run1bam) RUN1BAM="$2"; shift ;;
    --r1f1| --run1f1) RUN1F1="$2" ; shift ;;
    --r1f2| --run1f2) RUN1F2="$2" ; shift ;;
    --r1r| --run1run) RUN1RUN="$2" ; shift ;;
    --r1l| --run1lib) RUN1LIB="$2" ; shift ;;
    --r2b| --run2bam) RUN2BAM="$2"; shift ;;
    --r2f1| --run2f1) RUN2F1="$2" ; shift ;;
    --r2f2| --run2f2) RUN2F2="$2" ; shift ;;
    --r2r| --run2run) RUN2RUN="$2" ; shift ;;
    --r2l| --run2lib) RUN2LIB="$2" ; shift ;;
    --r3b| --run3bam) RUN3BAM="$2"; shift ;;
    --r3f1| --run3f1) RUN3F1="$2" ; shift ;;
    --r3f2| --run3f2) RUN3F2="$2" ; shift ;;
    --r3r| --run3run) RUN3RUN="$2" ; shift ;;
    --r3l| --run3lib) RUN3LIB="$2" ; shift ;;
    --r4b| --run4bam) RUN4BAM="$2"; shift ;;
    --r4f1| --run4f1) RUN4F1="$2" ; shift ;;
    --r4f2| --run4f2) RUN4F2="$2" ; shift ;;
    --r4r| --run4run) RUN4RUN="$2" ; shift ;;
    --r4l| --run4lib) RUN4LIB="$2" ; shift ;;
    --r5b| --run5bam) RUN5BAM="$2"; shift ;;
    --r5f1| --run5f1) RUN5F1="$2" ; shift ;;
    --r5f2| --run5f2) RUN5F2="$2" ; shift ;;
    --r5r| --run5run) RUN5RUN="$2" ; shift ;;
    --r5l| --run5lib) RUN5LIB="$2" ; shift ;;
    --r6b| --run6bam) RUN6BAM="$2"; shift ;;
    --r6f1| --run6f1) RUN6F1="$2" ; shift ;;
    --r6f2| --run6f2) RUN6F2="$2" ; shift ;;
    --r6r| --run6run) RUN6RUN="$2" ; shift ;;
    --r6l| --run6lib) RUN6LIB="$2" ; shift ;;
    --r7b| --run7bam) RUN7BAM="$2"; shift ;;
    --r7f1| --run7f1) RUN7F1="$2" ; shift ;;
    --r7f2| --run7f2) RUN7F2="$2" ; shift ;;
    --r7r| --run7run) RUN7RUN="$2" ; shift ;;
    --r7l| --run7lib) RUN7LIB="$2" ; shift ;;
    --r8b| --run8bam) RUN8BAM="$2"; shift ;;
    --r8f1| --run8f1) RUN8F1="$2" ; shift ;;
    --r8f2| --run8f2) RUN8F2="$2" ; shift ;;
    --r8r| --run8run) RUN8RUN="$2" ; shift ;;
    --r8l| --run8lib) RUN8LIB="$2" ; shift ;;
    --r9b| --run9bam) RUN9BAM="$2"; shift ;;
    --r9f1| --run9f1) RUN9F1="$2" ; shift ;;
    --r9f2| --run9f2) RUN9F2="$2" ; shift ;;
    --r9r| --run9run) RUN9RUN="$2" ; shift ;;
    --r9l| --run9lib) RUN9LIB="$2" ; shift ;;
    --bed) BED="$2" ; shift ;;
    --bed_flank) BED_FLANK="$2" ; shift ;;
    --center) CENTER="$2" ; shift ;;
    --nomt) NOMT=1 ;;
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
if [ -z $MODE ]; then echo "ERROR: missing mandatory option: --mode (bam or fastq raw inputs) must be specified"; FOUND_ERROR=1; fi
if [[ $MODE != "bam" && $MODE != "fastq" ]]; then echo "ERROR: --mode must be one of \"bam\" or \"fastq\""; FOUND_ERROR=1; fi
if [[ $MODE == "bam" && -z $RUN1BAM ]]; then echo "ERROR: must supply at least 1 raw bam file if using \"--mode bam\""; FOUND_ERROR=1; fi
if [[ $MODE == "fastq" && (-z $RUN1F1 || -z $RUN1F2) ]]; then echo "ERROR: must supply at least 1 set of raw forward and reverse fastq files if using \"--mode fastq\""; FOUND_ERROR=1; fi
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
if [[ -n $ACCOUNT ]]; then ACCOUNT="--account=$ACCOUNT"; fi
OLD_PATH=`pwd -P`
[ `echo $DIR | grep -P "^\."` ] && DIR=$OLD_PATH/$DIR; 
DIR=$(cd $DIR 2>/dev/null && pwd -P);
if [ ! -d $DIR ]; then echo "ERROR: wrong mandatory option: -d (run directory $DIR) is invalid"; FOUND_ERROR=1; fi

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

checkRawBam() {
  read bam run lib <<< $@
  if [[ -n $bam || -n $run || -n $lib ]]; then
    if [[ -n $bam && -n $run && -n $lib ]]; then
      [ `echo $bam | grep -P "^\."` ] && bam=$OLD_PATH/$bam;
      bam=$(cd $(dirname $bam) && pwd -P)/$(basename $bam);
      if [ ! -f $bam ]; then echo "ERROR: invalid bam file option: '$bam' could not be found"; FOUND_ERROR=1; fi
      [ `echo $reverse | grep -P "^\."` ] && reverse=$OLD_PATH/$reverse;
    else echo "ERROR: Please define all variables for RUN $run"; FOUND_ERROR=1; fi
  fi
}

if [[ $MODE == "fastq" ]]; then
  checkRawFastq $RUN1F1 $RUN1F2 $RUN1RUN $RUN1LIB
  checkRawFastq $RUN2F1 $RUN2F2 $RUN2RUN $RUN2LIB
  checkRawFastq $RUN3F1 $RUN3F2 $RUN3RUN $RUN3LIB
  checkRawFastq $RUN4F1 $RUN4F2 $RUN4RUN $RUN4LIB
  checkRawFastq $RUN5F1 $RUN5F2 $RUN5RUN $RUN5LIB
  checkRawFastq $RUN6F1 $RUN6F2 $RUN6RUN $RUN6LIB
  checkRawFastq $RUN7F1 $RUN7F2 $RUN7RUN $RUN7LIB
  checkRawFastq $RUN8F1 $RUN8F2 $RUN8RUN $RUN8LIB
  checkRawFastq $RUN9F1 $RUN9F2 $RUN9RUN $RUN9LIB
elif [[ $MODE == "bam" ]]; then
  checkRawBam $RUN1BAM $RUN1RUN $RUN1LIB
  checkRawBam $RUN2BAM $RUN2RUN $RUN2LIB
  checkRawBam $RUN3BAM $RUN3RUN $RUN3LIB
  checkRawBam $RUN4BAM $RUN4RUN $RUN4LIB
  checkRawBam $RUN5BAM $RUN5RUN $RUN5LIB
  checkRawBam $RUN6BAM $RUN6RUN $RUN6LIB
  checkRawBam $RUN7BAM $RUN7RUN $RUN7LIB
  checkRawBam $RUN8BAM $RUN8RUN $RUN8LIB
  checkRawBam $RUN9BAM $RUN9RUN $RUN9LIB
fi

if (( $FOUND_ERROR )); then echo "Please check options and try again"; exit 42; fi

export SAMPLE DIR THREADS WALLTIME CAPTURE_KIT BED BED_FLANKS RUN1BAM1 RUN1F1 RUN1F2 RUN1RUN RUN1LIB RUN2BAM RUN2F1 RUN2F2 RUN2RUN RUN2LIB RUN3BAM RUN3F1 RUN3F2 RUN3RUN RUN3LIB RUN4BAM RUN4F1 RUN4F2 RUN4RUN RUN4LIB RUN5BAM RUN5F1 RUN5F2 RUN5RUN RUN5LIB RUN6BAM RUN6F1 RUN6F2 RUN6RUN RUN6LIB RUN7BAM RUN7F1 RUN7F2 RUN7RUN RUN7LIB RUN8BAM RUN8F1 RUN8F2 RUN8RUN RUN8LIB RUN9BAM RUN9F1 RUN9F2 RUN9RUN RUN9LIB VERBOSE SFS KEEP_INTERMEDIATES ADAPTERS CENTER

# ===============================================
# QUEUE THE JOB IN CLUSTER
# ===============================================
cd $DIR; # for job outputs

# ALIGNMENTS
unset MERGE_V37_SORT_SAMS
unset MERGE_MT_SORT_SAMS
unset ALIGN_V37_IDS
unset ALIGN_MT_IDS

alignFastq(){
  read input forward reverse run lib <<< $@
  if [[ -n $forward && -n $reverse && -n $run && -n $lib ]]; then
    V37SAM=${input}-${lib}-${run}.paired.aligned.v37.sam
    MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${V37SAM}"
    MTSAM=${input}-${lib}-${run}.paired.aligned.MT.sam
    MERGE_MT_SORT_SAMS="${MERGE_MT_SORT_SAMS}@INPUT=${MTSAM}"
    if [[ $DISCARD_UNPAIRED_READS == 0 ]]; then 
      V37SAMUNPAIRED=${input}-${lib}-${run}.unpaired.aligned.v37.sam;
      MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${V37SAMUNPAIRED}"
    fi
    align $input $forward $reverse $lib $run $V37SAM $MTSAM $V37SAMUNPAIRED;
    if [[ ${DEPEND_V37_ALN} ]]; then ALIGN_V37_IDS=${ALIGN_V37_IDS}${DEPEND_V37_ALN}; fi
    if [[ ${DEPEND_MT_ALN} ]]; then ALIGN_MT_IDS=${ALIGN_MT_IDS}${DEPEND_MT_ALN}; fi
  fi
}

alignBam() {
  read input bam run lib <<< $@
  if [[ -n $bam && -n $run && -n $lib ]]; then
    V37SAM=${input}-${lib}-${run}.paired.aligned.v37.sam
    MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${V37SAM}"
    MTSAM=${input}-${lib}-${run}.paired.aligned.MT.sam
    MERGE_MT_SORT_SAMS="${MERGE_MT_SORT_SAMS}@INPUT=${MTSAM}"
    if [[ $DISCARD_UNPAIRED_READS == 0 ]]; then 
      V37SAMUNPAIRED=${input}-${lib}-${run}.unpaired.aligned.v37.sam;
      MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${V37SAMUNPAIRED}"
    fi
    align $input ${bam}_read1.fastq.gz ${bam}_read2.fastq.gz $lib $run $V37SAM $MTSAM $V37SAMUNPAIRED;
    if [[ ${DEPEND_V37_ALN} ]]; then ALIGN_V37_IDS=${ALIGN_V37_IDS}${DEPEND_V37_ALN}; fi
    if [[ ${DEPEND_MT_ALN} ]]; then ALIGN_MT_IDS=${ALIGN_MT_IDS}${DEPEND_MT_ALN}; fi
  fi
}

unset DEPEND;
if [[ $MODE == "fastq" ]]; then
  alignFastq input1 $RUN1F1 $RUN1F2 $RUN1RUN $RUN1LIB
  alignFastq input2 $RUN2F1 $RUN2F2 $RUN2RUN $RUN2LIB
  alignFastq input3 $RUN3F1 $RUN3F2 $RUN3RUN $RUN3LIB
  alignFastq input4 $RUN4F1 $RUN4F2 $RUN4RUN $RUN4LIB
  alignFastq input5 $RUN5F1 $RUN5F2 $RUN5RUN $RUN5LIB
  alignFastq input6 $RUN6F1 $RUN6F2 $RUN6RUN $RUN6LIB
  alignFastq input7 $RUN7F1 $RUN7F2 $RUN7RUN $RUN7LIB
  alignFastq input8 $RUN8F1 $RUN8F2 $RUN8RUN $RUN8LIB
  alignFastq input9 $RUN9F1 $RUN9F2 $RUN9RUN $RUN9LIB
elif [[ $MODE == "bam" ]]; then
  alignBam input1 $RUN1BAM $RUN1RUN $RUN1LIB
  alignBam input2 $RUN2BAM $RUN2RUN $RUN2LIB
  alignBam input3 $RUN3BAM $RUN3RUN $RUN3LIB
  alignBam input4 $RUN4BAM $RUN4RUN $RUN4LIB
  alignBam input5 $RUN5BAM $RUN5RUN $RUN5LIB
  alignBam input6 $RUN6BAM $RUN6RUN $RUN6LIB
  alignBam input7 $RUN7BAM $RUN7RUN $RUN7LIB
  alignBam input8 $RUN8BAM $RUN8RUN $RUN8LIB
  alignBam input9 $RUN9BAM $RUN9RUN $RUN9LIB
fi

echo "*** Post Alignment"

# MERGE SAM FILES
NAME=pipeline.generic.part_2_0.combine_sample_libs.primary_aln
STEP=2_0g
export THREADS=${THREADS_ARRAY[$STEP]}
export MAX_MEM=$[$THREADS*4]
export WALLTIME=${WALLTIME_ARRAY[$STEP]}
export DONE=.${SAMPLE}.${NAME}.v37.done
echo $DONE >> $EXPECTED_DONE_FILES
MERGE_V37_SORT_BAM_PRIMARY_ALN=${SAMPLE}.aligned.v37.merged.sorted.bam
MERGE_V37_SORT_BAM_NON_PRIMARY_ALN=${SAMPLE}.aligned.v37.merged.sorted.non_pri_aln.bam
if [[ ! -z $ALIGN_V37_IDS ]]; then DEPEND="--dependency=afterok${ALIGN_V37_IDS}"; else unset DEPEND; fi
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  export MERGE_SORT_SAMS=${MERGE_V37_SORT_SAMS}
  export MERGE_SORT_BAM=${MERGE_V37_SORT_BAM_PRIMARY_ALN}
  export MERGE_SORT_BAM_NON_PRIMARY_ALN=${MERGE_V37_SORT_BAM_NON_PRIMARY_ALN}
  v37_merge_sample_libs="$QUEUE $ACCOUNT \
    --ntasks-per-node=${THREADS} \
    --mem-per-cpu=4g \
    --time=${WALLTIME} \
    --job-name ${SAMPLE}.${NAME}.v37 \
    --output %x.o%j \
    $DEPEND \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub"
  v37_merge_sample_libs=$($v37_merge_sample_libs | grep -oP "\d+")
  echo "[Q] Merge Sort v37  : $v37_merge_sample_libs [$DEPEND]"
  V37_MERGE_DEPEND="--dependency=afterok:$v37_merge_sample_libs";
else 
  echo "[D] Merge Sort v37"
fi

# GATK EXOME PRE-PROCESSING
NAME=pipeline.exome.part_2_1.GATK_data_processing.in_bam
STEP=2_1g
export THREADS=${THREADS_ARRAY[$STEP]}
export MAX_MEM=$[$THREADS*4]
export WALLTIME=${WALLTIME_ARRAY[$STEP]}
export DONE=.${SAMPLE}.${NAME}.v37.done
echo $DONE >> $EXPECTED_DONE_FILES
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  export ALIGNED_BAM=${MERGE_V37_SORT_BAM_PRIMARY_ALN}
  export USE_INDELS=1
  export GENOME=v37
  gatk_v37="$QUEUE $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem-per-cpu=4g \
    --time=${WALLTIME} \
    --job-name ${SAMPLE}.${NAME}.v37 \
    --output %x.o%j \
    $V37_MERGE_DEPEND \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub"
  gatk_v37=$($gatk_v37 | grep -oP "\d+")
  echo "[Q] GATK Exome      : $gatk_v37 [$V37_MERGE_DEPEND]"
  GATK_EXOME_DEPEND="--dependency=afterok:$gatk_v37";
else echo "[D] GATK Exome"
fi

# MTCHONDRIA PROCESSING
# MERGE SAM FILES
NAME=pipeline.generic.part_2_0.combine_sample_libs
STEP=2_0m
export THREADS=${THREADS_ARRAY[$STEP]}
export MAX_MEM=$[$THREADS*4]
export WALLTIME=${WALLTIME_ARRAY[$STEP]}
export DONE=.${SAMPLE}.${NAME}.MT.done
echo $DONE >> $EXPECTED_DONE_FILES
export MERGE_MT_SORT_BAM=${SAMPLE}.aligned.MT.merged.sorted.bam
if [[ ! -z $ALIGN_MT_IDS ]]; then DEPEND="--dependency=afterok${ALIGN_MT_IDS}"; else unset DEPEND; fi
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  export MERGE_SORT_SAMS=${MERGE_MT_SORT_SAMS}
  export MERGE_SORT_BAM=${MERGE_MT_SORT_BAM}
  MT_merge_sample_libs="$QUEUE $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem-per-cpu=4g \
    --time=${WALLTIME} \
    --job-name ${SAMPLE}.${NAME}.MT\
    --output %x.o%j \
    $DEPEND \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub"
  MT_merge_sample_libs=$($MT_merge_sample_libs | grep -oP "\d+")
  echo "[Q] Merge Sort MT   : $MT_merge_sample_libs [$DEPEND]"
  MT_MERGE_DEPEND="--dependency=afterok:$MT_merge_sample_libs";
else
  echo "[D] Merge Sort MT"
fi

# GATK MT PRE-PROCESSING
NAME=pipeline.exome.part_2_1.GATK_data_processing.in_bam
STEP=2_1m
export THREADS=${THREADS_ARRAY[$STEP]}
export MAX_MEM=$[$THREADS*4]
export WALLTIME=${WALLTIME_ARRAY[$STEP]}
export DONE=.${SAMPLE}.${NAME}.MT.done
echo $DONE >> $EXPECTED_DONE_FILES
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  export ALIGNED_BAM=${MERGE_MT_SORT_BAM}
  export USE_INDELS=0
  export GENOME=MT
  gatk_MT="$QUEUE $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem-per-cpu=4g \
    --time=${WALLTIME} \
    --job-name ${SAMPLE}.${NAME}.MT\
    --output %x.o%j \
    $MT_MERGE_DEPEND \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub"
  gatk_MT=$($gatk_MT | grep -oP "\d+")
  echo "[Q] GATK MT         : $gatk_MT [$MT_MERGE_DEPEND]"
  GATK_MT_DEPEND="--dependency=afterok:$gatk_MT";
else echo "[D] GATK MT"
fi

## FILTER MT
NAME=pipeline.exome.part_2_2.MT_Filtering
STEP=2_2
export THREADS=${THREADS_ARRAY[$STEP]}
export MAX_MEM=$[$THREADS*4]
export WALLTIME=${WALLTIME_ARRAY[$STEP]}
export DONE=.${SAMPLE}.${NAME}.done
echo $DONE >> $EXPECTED_DONE_FILES
if [[ ! `find ./ -maxdepth 1 -name "$DONE"` || $SFS -eq 1  ]]; then
  sample_filter_MT="$QUEUE $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem-per-cpu=4g \
    --time=${WALLTIME} \
    --job-name ${SAMPLE}.${NAME}\
    --output %x.o%j \
    $GATK_MT_DEPEND \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub"
  sample_filter_MT=$($sample_filter_MT | grep -oP "\d+")
  echo "[Q] Filter MT       : $sample_filter_MT [$GATK_MT_DEPEND]"
  FILTER_MT_DEPEND="--dependency=afterok:$sample_filter_MT";
else echo "[D] Filter MT"
fi

## MT STATS AND QC
NAME=pipeline.exome.part_3_0.MT_stats_QC
STEP=3_0
export THREADS=${THREADS_ARRAY[$STEP]}
export MAX_MEM=$[$THREADS*4]
export WALLTIME=${WALLTIME_ARRAY[$STEP]}
export DONE=.${SAMPLE}.${NAME}.done
echo $DONE >> $EXPECTED_DONE_FILES
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  export CAPTURE_KIT=${CAPTURE_KIT}
  MT_sample_stats_QC="$QUEUE $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem-per-cpu=4g \
    --time=${WALLTIME} \
    --job-name ${SAMPLE}.${NAME}\
    --output %x.o%j \
    ${FILTER_MT_DEPEND} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub"
  MT_sample_stats_QC=$($MT_sample_stats_QC | grep -oP "\d+")
  echo "[Q] MT Stats and QC : $MT_sample_stats_QC [$FILTER_MT_DEPEND]"
  MT_STATS_QC_DEPEND="--dependency=afterok:$MT_sample_stats_QC";
else echo "[D] MT Stats and QC"
fi

## Exome STATS AND QC
NAME=pipeline.exome.part_3_1.Exome_stats_QC
STEP=3_1
export THREADS=${THREADS_ARRAY[$STEP]}
export MAX_MEM=$[$THREADS*4]
export WALLTIME=${WALLTIME_ARRAY[$STEP]}
export DONE=.${SAMPLE}.${NAME}.done
echo $DONE >> $EXPECTED_DONE_FILES
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  export CAPTURE_KIT=${CAPTURE_KIT}
  Exome_sample_stats_QC="$QUEUE $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem-per-cpu=4g \
    --time=${WALLTIME} \
    --job-name ${SAMPLE}.${NAME}\
    --output %x.o%j \
    $GATK_EXOME_DEPEND \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub"
  Exome_sample_stats_QC=$($Exome_sample_stats_QC | grep -oP "\d+")
  echo "[Q] Ex Stats and QC : $Exome_sample_stats_QC [$GATK_EXOME_DEPEND]"
  EXOME_STATS_QC_DEPEND="--dependency=afterok:$Exome_sample_stats_QC";
else echo "[D] Exome Stats and QC"
fi

### GATK UG total recall 
#NAME=pipeline.exome.part_3_5.GATK_UG_total_recall
#STEP=3_5
#export THREADS=${THREADS_ARRAY[$STEP]}
#export MAX_MEM=$[$THREADS*4]
#export WALLTIME=${WALLTIME_ARRAY[$STEP]}
#export DONE=.${SAMPLE}.${NAME}.done
#echo $DONE >> $EXPECTED_DONE_FILES
#if [[ ! -z $CALLING_BED ]]; then export BEDFILE=$CALLING_BED; else unset CALLING_BED_ARG; fi
#if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
#  gatk_UG_total_recall="$QUEUE $ACCOUNT  \
#    --ntasks-per-node=${THREADS} \
#    --mem-per-cpu=4g \
#    --time=${WALLTIME} \
#    --job-name ${SAMPLE}.${NAME}\
#    --output %x.o%j \
#    $GATK_EXOME_DEPEND \
#    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub"
#  gatk_UG_total_recall=$($gatk_UG_total_recall | grep -oP "\d+")
#  echo "[Q] GATK UG         : $gatk_UG_total_recall [$GATK_EXOME_DEPEND]"
#else echo "[D] GATK UG "
#fi

## GATK HC total recall 
NAME=pipeline.exome.part_3_6.GATK_HC
STEP=3_6
export THREADS=${THREADS_ARRAY[$STEP]}
export MAX_MEM=$[$THREADS*4]
export WALLTIME=${WALLTIME_ARRAY[$STEP]}
export DONE=.${SAMPLE}.${NAME}.done
echo $DONE >> $EXPECTED_DONE_FILES
if [[ ! -z $CALLING_BED ]]; then export BEDFILE=$CALLING_BED; else unset CALLING_BED_ARG; fi
if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  gatk_HC="$QUEUE $ACCOUNT  \
    --ntasks-per-node=${THREADS} \
    --mem-per-cpu=4g \
    --time=${WALLTIME} \
    --job-name ${SAMPLE}.${NAME}\
    --output %x.o%j \
    $GATK_EXOME_DEPEND \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub"
  gatk_HC=$($gatk_HC | grep -oP "\d+")
  echo "[Q] GATK HC         : $gatk_HC [$GATK_EXOME_DEPEND]"
else echo "[D] GATK HC"
fi

cd $OLD_PATH
exit 0
