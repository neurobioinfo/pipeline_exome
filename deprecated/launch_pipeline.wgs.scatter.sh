#!/bin/bash
##DONE: update QUEUE assignation
##NOPE: add functionality for bam inputs (Illumina) -> uses qsub 0_1 to generate fastqs (includes update of usage and getopts --> do fastq version first, make sure it works, then add bam step to a new script "wgs.bams.scatter"
##DONE: add -A | --account option to getops (and usage)
##DONE: fix depend syntax in job executions
##DONE: fix 2_1 launch + qsub script
##DONE: fix 3_1 launch + qsub script
##DONE: create 3_2 launch + qsub script (variant calling)
##DONE: fix 4_1 launch + qsub script
##DONE: add ability to handle up to 10 sets of input gz files
##DONE: update usage function
##TODO: cleanup function.library - remove old functions

# ===============================================
# default variables values
# ===============================================
unset SAMPLE DIR THREADS MAX_MEM PIPELINE_HOME WALLTIME EXTRA_CONF RUN1F1 RUN1F2 RUN1RG RUN2F1 RUN2F2 RUN2RG RUN3F1 RUN3F2 RUN3RG RUN4F1 RUN4F2 RUN4RG RUN5F1 RUN5F2 RUN5RG RUN6F1 RUN6F2 RUN6RG RUN7F1 RUN7F2 RUN7RG RUN8F1 RUN8F2 RUN8RG RUN9F1 RUN9F2 RUN9RG RUN10F1 RUN10F2 RUN10RG QUEUE ACCOUNT VERBOSE CENTER SFS KEEP_INTERMEDIATES

#Assuming script is in root of PIPELINE_HOME
export PIPELINE_HOME=$(cd $(dirname $0) && pwd -P)
export SFS=0
export KEEP_INTERMEDIATES=0
export DISCARD_UNPAIRED_READS=0

#set scheduler-specific variables
QUEUE=qsub                                                                           # default job scheduler: qsub
if [[ `which msub 2>/dev/null` ]]; then QUEUE=msub; fi                               # assign $QUEUE based on server's scheduler software
if [[ `uname -n` =~ ^gpc ]]; then QUEUE=qsub; fi;                                    # fix for scinet gpc system, which has msub, but diabled it in favor of qsub (i.e. `which msub` test will yield TRUE, and mess everything up)
if [[ $QUEUE == "msub" ]]; then dependMod="x=depend:"; else dependMod="depend="; fi  # assigns scheduler-specific dependency syntax, based on scheduler software

#create function to handle error messages
Usage() {
        echo
        echo -e "Usage:\t$0 [arguments]"
        echo -e "\tmandatory arguments:\n" \
          "\t\t-s  (--sample)           = sample name\n" \
          "\t\t-d  (--dir)              = run directory (where all the outputs will be printed)\n" \
          "\t\t-r[1-10]f1 (--run[1-10]f1) = fastq sequence run [1-10] forward\n" \
          "\t\t-r[1-10]f2 (--run[1-10]f2) = fastq sequence run [1-10] reverse\n" \
          "\t\t-r[1-10]rg (--run[1-10]rg) = fastq sequence run [1-10] read group\n" \
          "\t\t-x  (--extra)            = use extra config file to (re)define variables [CURRENT \"$EXTRA_CONF\"]" 
        echo -e "\toptional arguments:\n" \
          "\t\t-h  (--help)             = get the program options and exit\n" \
          "\t\t--bed                    = specify the bed file for QC metrics\n" \
          "\t\t--bed_flank              = specify the bed file with flanking for QC metrics\n" \
          "\t\t--sites                  = specify the sites to use for variant calling\n" \
          "\t\t-a  (--adapters)         = fasta sequence file of adapters\n" \
          "\t\t-A  (--account)          = account name\n" \
          "\t\t-m  (--max_mem)          = max memory on a node (in gigs) to use [CURRENT \"$MAX_MEM\"]\n" \
          "\t\t-t  (--threads)          = number of threads to use [CURRENT \"$THREADS\"]\n" \
          "\t\t-w  (--walltime)         = wall time for the scheduler (format HH:MM:SS) [CURRENT \"$WALLTIME\"]\n" \
          "\t\t-v  (--verbose)          = set verbosity level [CURRENT \"$VERBOSE\"]\n" \
          "\t\t--center                 = specify sequencing center [CURRENT \"$CENTER\"]\n" \
          "\t\t--start_from_scratch     = flag to overwrite existing files [CURRENT \"$SFS\"]\n" \
          "\t\t--discard_unpaired       = flag to discard unpaired reads after pre processing [CURRENT \"$DISCARD_UNPAIRED_READS\"]\n" \
          "\t\t--keep_intermediates     = flag to keep intermediates files [CURRENT \"$KEEP_INTERMEDIATES\"]" 
        echo
}

# ===============================================
# alignment function
# ===============================================
align() {
  local F1=$1
  local F2=$2
  local RUN=$3
  local SAM_V37=$4
  local SAM_V37_UNPAIRED=$5

#DETERMINE BASE
  export ENCODING=`${PIPELINE_HOME}/soft/bin/determine_base.perl $F1`;
  if [[ "$ENCODING" == "Cannot determine base" ]]; then echo "ERROR: Cannot determine base for $F1"; exit 42;
  else echo "*** ENCODING bases set to $ENCODING"; fi

  echo -e "*** ALIGNMENT v37\n*** R1:$F1\n*** R2:$F2\n*** RUN ID: $RUN"
#PRE PREPROCESSING
  FASTQC_OUTPUT_DIR=${DIR}/fastqc
  export TRIM_PAIRED_FASTQ1=${DIR}/${SAMPLE}.${RUN}.r1.fastq.trimmomaticPE.paired.gz
  export TRIM_PAIRED_FASTQ2=${DIR}/${SAMPLE}.${RUN}.r2.fastq.trimmomaticPE.paired.gz
  export TRIM_UNPAIRED_FASTQ1=${DIR}/${SAMPLE}.${RUN}.r1.fastq.trimmomaticPE.unpaired.gz
  export TRIM_UNPAIRED_FASTQ2=${DIR}/${SAMPLE}.${RUN}.r2.fastq.trimmomaticPE.unpaired.gz

  NAME=pipeline.generic.part_0_1.trimmomatic
  DONE=.${SAMPLE}.${NAME}.${RUN}.done
  if [[ ! `find ./ -name ${DONE}` || $SFS -eq 1 ]]; then
    trimm=`${QUEUE} \
      -l nodes=1:ppn=$THREADS \
      -l mem=${MAX_MEM}g \
      -l walltime=${WALLTIME} \
      -N ${SAMPLE}.${NAME}.${RUN} \
      -V \
      -v FASTQ1=$F1,FASTQ2=$F2,TRIM_PAIRED_FASTQ1=$TRIM_PAIRED_FASTQ1,TRIM_PAIRED_FASTQ2=$TRIM_PAIRED_FASTQ2,TRIM_UNPAIRED_FASTQ1=$TRIM_UNPAIRED_FASTQ1,TRIM_UNPAIRED_FASTQ2=$TRIM_UNPAIRED_FASTQ2,FASTQC_OUTPUT_DIR=$FASTQC_OUTPUT_DIR,DONE=$DONE,RUN=$RUN \
      ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
    trimm=$(echo $trimm)
    echo "[Q] Trimmomatic     : $trimm"
    DEPEND_TRIM="-W "$dependMod"afterok:$trimm"
  else echo "[D] Trimmomatic"
  fi

# ALIGN FORWARD READ
  SAI1=${TRIM_PAIRED_FASTQ1}.v37.sai
  NAME=pipeline.generic.part_1_1.align_aln.in_fastq
  DONE=.${SAMPLE}.${NAME}.${RUN}.v37.r1.done
  if [[ ! `find ./ -name ${DONE}` || $SFS -eq 1 ]]; then
    v37_align_r1=`$QUEUE $ACCOUNT \
      -l nodes=1:ppn=$THREADS \
      -l mem=${MAX_MEM}g \
      -l walltime=$WALLTIME \
      -N ${SAMPLE}.${NAME}.${RUN}.v37.r1 \
      $DEPEND_TRIM \
      -V \
      -v FASTQ=$TRIM_PAIRED_FASTQ1,SAI=$SAI1,GENOME=v37,DONE=$DONE \
      ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
    v37_align_r1=$(echo $v37_align_r1)
    echo "[Q] Align v37 r1    : $v37_align_r1 [$DEPEND_TRIM]"
  else echo "[D] Align v37 r1"
  fi

# ALIGN REVERSE READ
  SAI2=${TRIM_PAIRED_FASTQ2}.v37.sai
  NAME=pipeline.generic.part_1_1.align_aln.in_fastq
  DONE=.${SAMPLE}.${NAME}.${RUN}.v37.r2.done
  if [[ ! `find ./ -name ${DONE}` || $SFS -eq 1 ]]; then
    v37_align_r2=`$QUEUE $ACCOUNT  \
      -l nodes=1:ppn=$THREADS \
      -l mem=${MAX_MEM}g \
      -l walltime=$WALLTIME \
      -N ${SAMPLE}.${NAME}.${RUN}.v37.r2 \
      $DEPEND_TRIM \
      -V \
      -v FASTQ=$TRIM_PAIRED_FASTQ2,SAI=$SAI2,GENOME=v37,DONE=$DONE \
      ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
    v37_align_r2=$(echo $v37_align_r2)
    echo "[Q] Align v37 r2    : $v37_align_r2 [$DEPEND_TRIM]"
  else echo "[D] Align v37 r2"
  fi

#PARING ALIGNMENT (bwa sampe)
  NAME=pipeline.generic.part_1_2.align_sampe.in_fastq
  DONE=.${SAMPLE}.${NAME}.${RUN}.v37.done
  if [[ ! `find ./ -name ${DONE}` || $SFS -eq 1 ]]; then
    DEPEND=
    if [[ $v37_align_r1 && $v37_align_r2 ]]; then DEPEND="-W "$dependMod"afterok:$v37_align_r1:$v37_align_r2";
    elif [[ $v37_align_r1 ]]; then DEPEND="-W "$dependMod"afterok:$v37_align_r1";
    elif [[ $v37_align_r2 ]]; then DEPEND="-W "$dependMod"afterok:$v37_align_r2";
    fi
    v37_alignPE=`$QUEUE $ACCOUNT  \
      -l nodes=1:ppn=$THREADS \
      -l mem=${MAX_MEM}g \
      -l walltime=$WALLTIME \
      -N ${SAMPLE}.${NAME}.${RUN}.v37 \
      $DEPEND \
      -V \
      -v FASTQ1=$TRIM_PAIRED_FASTQ1,FASTQ2=$TRIM_PAIRED_FASTQ2,SAI1=$SAI1,SAI2=$SAI2,SAM=$SAM_V37,RUN=$RUN,GENOME=v37,DONE=$DONE \
      ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
    v37_alignPE=$(echo $v37_alignPE)
    echo "[Q] Align v37 PE    : $v37_alignPE [$DEPEND]"
    DEPEND_V37_ALN=${DEPEND_V37_ALN}:${v37_alignPE}
  else echo "[D] Align v37 PE"
  fi

if [[ $DISCARD_UNPAIRED_READS -eq 0 ]]; then
  #ALIGN UNPAIRED (bwa aln + samse)
  NAME=pipeline.generic.part_1_3.aln_samse.in_fastq
  DONE=.${SAMPLE}.${NAME}.${RUN}.v37.done
  if [[ ! `find ./ -name ${DONE}` || $SFS -eq 1 ]]; then
    v37_aln_samse=`$QUEUE $ACCOUNT  \
      -l nodes=1:ppn=$THREADS \
      -l mem=${MAX_MEM}g \
      -l walltime=$WALLTIME \
      -N ${SAMPLE}.${NAME}.${RUN}.v37 \
      $DEPEND_TRIM \
      -V \
      -v FASTQ1=$TRIM_UNPAIRED_FASTQ1,FASTQ2=$TRIM_UNPAIRED_FASTQ2,SAM=$SAM_V37_UNPAIRED,RUN=$RUN,GENOME=v37,DONE=$DONE \
      ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
    v37_aln_samse=$(echo $v37_aln_samse)
    echo "[Q] Align v37 samse: $v37_aln_samse [$DEPEND_TRIM]"
    DEPEND_V37_ALN=${DEPEND_V37_ALN}:${v37_aln_samse}
  else echo "[D] Align v37 samse"
  fi
fi
}


# ===============================================
# PARSING ARGUMENTS
# ===============================================

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt --name pipeline --alternative --unquoted --options hs:d:t:m:vw:c:a:A:x: --longoptions sample:,dir:,threads:,max_mem:,walltime:,adapters:,r1f1:,r1f2:,r1rg:,r2f1:,r2f2:,r2rg:,r3f1:,r3f2:,r3rg:,r4f1:,r4f2:,r4rg:,r5f1:,r5f2:,r5rg:,r6f1:,r6f2:,r6rg:,r7f1:,r7f2:,r7rg:,r8f1:,r8f2:,r8rg:,r9f1:,r9f2:,r9rg:,r10f1:,r10f2:,r10rg:,account:,extra:,center:,keep_intermediates,discard_unpaired,start_from_scratch,verbose,help -- "$@")
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
    -A| --account) ACCOUNT="$2"; shift ;;
    -v| --verbose) VERBOSE=1 ;;
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
    --r10f1| --run10f1) RUN10F1="$2" ; shift ;;
    --r10f2| --run10f2) RUN10F2="$2" ; shift ;;
    --r10rg| --run10rg) RUN10RG="$2" ; shift ;;
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

# ===============================================
# CHECKING VARIABLES
# ===============================================
#check to ensure all mandatory arguments have been entered
FOUND_ERROR=0
if [ -z $SAMPLE ]; then echo "ERROR: missing mandatory option: -s (sample name) must be specified"; FOUND_ERROR=1; fi
if [ -z $DIR ]; then echo "ERROR: missing mandatory option: -d (--dir) must be specified"; FOUND_ERROR=1; fi
if [ ! -d $DIR ]; then echo "ERROR: mandatory option: -d (--dir) is invalid - directory must exist"; FOUND_ERROR=1; fi
#if [ -z $ADAPTERS ]; then echo "ERROR: missing mandatory option: -a (adapters sequence) must be specified"; FOUND_ERROR=1;
#else
#  ADAPTERS=$(cd $(dirname $ADAPTERS) && pwd -P)/$(basename $ADAPTERS);
#  if [ ! -f $ADAPTERS ]; then echo "ERROR: invalid adapters file option: '$ADAPTERS' could not be found"; FOUND_ERROR=1; fi
#fi
if [ -z $EXTRA_CONF ]; then echo "ERROR: missing mandatory option: --extra must be specified"; FOUND_ERROR=1; fi

# ===============================================
#setting path and execution variables
# ===============================================
if [[ -n $ACCOUNT ]]; then ACCOUNT="-A $ACCOUNT"; fi
OLD_PATH=`pwd -P`
[ `echo $DIR | grep -P "^\."` ] && DIR=$OLD_PATH/$DIR;
DIR=$(cd $DIR 2>/dev/null && pwd -P);
if [ ! -d $DIR ]; then echo "ERROR: wrong mandatory option: -d (run directory $DIR) is invalid"; FOUND_ERROR=1; fi

if [[ -n $RUN1F1 || -n $RUN1F2 || -n $RUN1RG ]]; then
  if [[ -n $RUN1F1 && -n $RUN1F2 && -n $RUN1RG ]]; then
    [ `echo $RUN1F1 | grep -P "^\."` ] && RUN1F1=$OLD_PATH/$RUN1F1;
    RUN1F1=$(cd $(dirname $RUN1F1) && pwd -P)/$(basename $RUN1F1);
    if [ ! -f $RUN1F1 ]; then echo "ERROR: invalid RUN1F1 file option: '$RUN1F1' could not be found"; FOUND_ERROR=1; fi
    [ `echo $RUN1F2 | grep -P "^\."` ] && RUN1F2=$OLD_PATH/$RUN1F2;
    RUN1F2=$(cd $(dirname $RUN1F2) && pwd -P)/$(basename $RUN1F2);
    if [ ! -f $RUN1F2 ]; then echo "ERROR: invalid RUN1F2 file option: '$RUN1F2' could not be found"; FOUND_ERROR=1; fi
  else echo "ERROR: Please define all variables for RUN1"; FOUND_ERROR=1; fi
fi
if [[ -n $RUN2F1 || -n $RUN2F2 || -n $RUN2RG ]]; then
  if [[ -n $RUN2F1 && -n $RUN2F2 && -n $RUN2RG ]]; then
    [ `echo $RUN2F1 | grep -P "^\."` ] && RUN2F1=$OLD_PATH/$RUN2F1;
    RUN2F1=$(cd $(dirname $RUN2F1) && pwd -P)/$(basename $RUN2F1);
    if [ ! -f $RUN2F1 ]; then echo "ERROR: invalid RUN2F1 file option: '$RUN2F1' could not be found"; FOUND_ERROR=1; fi
    [ `echo $RUN2F2 | grep -P "^\."` ] && RUN2F2=$OLD_PATH/$RUN2F2;
    RUN2F2=$(cd $(dirname $RUN2F2) && pwd -P)/$(basename $RUN2F2);
    if [ ! -f $RUN2F2 ]; then echo "ERROR: invalid RUN2F2 file option: '$RUN2F2' could not be found"; FOUND_ERROR=1; fi
  else echo "ERROR: Please define all variables for RUN2"; FOUND_ERROR=1; fi
fi
if [[ -n $RUN3F1 || -n $RUN3F2 || -n $RUN3RG ]]; then
  if [[ -n $RUN3F1 && -n $RUN3F2 && -n $RUN3RG ]]; then
    [ `echo $RUN3F1 | grep -P "^\."` ] && RUN3F1=$OLD_PATH/$RUN3F1;
    RUN3F1=$(cd $(dirname $RUN3F1) && pwd -P)/$(basename $RUN3F1);
    if [ ! -f $RUN3F1 ]; then echo "ERROR: invalid RUN3F1 file option: '$RUN3F1' could not be found"; FOUND_ERROR=1; fi
    [ `echo $RUN3F2 | grep -P "^\."` ] && RUN3F2=$OLD_PATH/$RUN3F2;
    RUN3F2=$(cd $(dirname $RUN3F2) && pwd -P)/$(basename $RUN3F2);
    if [ ! -f $RUN3F2 ]; then echo "ERROR: invalid RUN3F2 file option: '$RUN3F2' could not be found"; FOUND_ERROR=1; fi
  else echo "ERROR: Please define all variables for RUN3"; FOUND_ERROR=1; fi
fi
if [[ -n $RUN4F1 || -n $RUN4F2 || -n $RUN4RG ]]; then
  if [[ -n $RUN4F1 && -n $RUN4F2 && -n $RUN4RG ]]; then
    [ `echo $RUN4F1 | grep -P "^\."` ] && RUN4F1=$OLD_PATH/$RUN4F1;
    RUN4F1=$(cd $(dirname $RUN4F1) && pwd -P)/$(basename $RUN4F1);
    if [ ! -f $RUN4F1 ]; then echo "ERROR: invalid RUN4F1 file option: '$RUN4F1' could not be found"; FOUND_ERROR=1; fi
    [ `echo $RUN4F2 | grep -P "^\."` ] && RUN4F2=$OLD_PATH/$RUN4F2;
    RUN4F2=$(cd $(dirname $RUN4F2) && pwd -P)/$(basename $RUN4F2);
    if [ ! -f $RUN4F2 ]; then echo "ERROR: invalid RUN4F2 file option: '$RUN4F2' could not be found"; FOUND_ERROR=1; fi
  else echo "ERROR: Please define all variables for RUN4"; FOUND_ERROR=1; fi
fi
if [[ -n $RUN5F1 || -n $RUN5F2 || -n $RUN5RG ]]; then
  if [[ -n $RUN5F1 && -n $RUN5F2 && -n $RUN5RG ]]; then
    [ `echo $RUN5F1 | grep -P "^\."` ] && RUN5F1=$OLD_PATH/$RUN5F1;
    RUN5F1=$(cd $(dirname $RUN5F1) && pwd -P)/$(basename $RUN5F1);
    if [ ! -f $RUN5F1 ]; then echo "ERROR: invalid RUN5F1 file option: '$RUN5F1' could not be found"; FOUND_ERROR=1; fi
    [ `echo $RUN5F2 | grep -P "^\."` ] && RUN5F2=$OLD_PATH/$RUN5F2;
    RUN5F2=$(cd $(dirname $RUN5F2) && pwd -P)/$(basename $RUN5F2);
    if [ ! -f $RUN5F2 ]; then echo "ERROR: invalid RUN5F2 file option: '$RUN5F2' could not be found"; FOUND_ERROR=1; fi
  else echo "ERROR: Please define all variables for RUN5"; FOUND_ERROR=1; fi
fi
if [[ -n $RUN6F1 && -n $RUN6F2 && -n $RUN6RG ]]; then
  if [[ -n $RUN6F1 && -n $RUN6F2 && -n $RUN6RG ]]; then
    [ `echo $RUN6F1 | grep -P "^\."` ] && RUN6F1=$OLD_PATH/$RUN6F1;
    RUN6F1=$(cd $(dirname $RUN6F1) && pwd -P)/$(basename $RUN6F1);
    if [ ! -f $RUN6F1 ]; then echo "ERROR: invalid RUN6F1 file option: '$RUN6F1' could not be found"; FOUND_ERROR=1; fi
    [ `echo $RUN6F2 | grep -P "^\."` ] && RUN6F2=$OLD_PATH/$RUN6F2;
    RUN6F2=$(cd $(dirname $RUN6F2) && pwd -P)/$(basename $RUN6F2);
    if [ ! -f $RUN6F2 ]; then echo "ERROR: invalid RUN6F2 file option: '$RUN6F2' could not be found"; FOUND_ERROR=1; fi
  else echo "ERROR: Please define all variables for RUN6"; FOUND_ERROR=1; fi
fi
if [[ -n $RUN7F1 || -n $RUN7F2 || -n $RUN7RG ]]; then
  if [[ -n $RUN7F1 && -n $RUN7F2 && -n $RUN7RG ]]; then
    [ `echo $RUN7F1 | grep -P "^\."` ] && RUN7F1=$OLD_PATH/$RUN7F1;
    RUN7F1=$(cd $(dirname $RUN7F1) && pwd -P)/$(basename $RUN7F1);
    if [ ! -f $RUN7F1 ]; then echo "ERROR: invalid RUN7F1 file option: '$RUN7F1' could not be found"; FOUND_ERROR=1; fi
    [ `echo $RUN7F2 | grep -P "^\."` ] && RUN7F2=$OLD_PATH/$RUN7F2;
    RUN7F2=$(cd $(dirname $RUN7F2) && pwd -P)/$(basename $RUN7F2);
    if [ ! -f $RUN7F2 ]; then echo "ERROR: invalid RUN7F2 file option: '$RUN7F2' could not be found"; FOUND_ERROR=1; fi
  else echo "ERROR: Please define all variables for RUN7"; FOUND_ERROR=1; fi
fi
if [[ -n $RUN8F1 || -n $RUN8F2 || -n $RUN8RG ]]; then
  if [[ -n $RUN8F1 && -n $RUN8F2 && -n $RUN8RG ]]; then
    [ `echo $RUN8F1 | grep -P "^\."` ] && RUN8F1=$OLD_PATH/$RUN8F1;
    RUN8F1=$(cd $(dirname $RUN8F1) && pwd -P)/$(basename $RUN8F1);
    if [ ! -f $RUN8F1 ]; then echo "ERROR: invalid RUN8F1 file option: '$RUN8F1' could not be found"; FOUND_ERROR=1; fi
    [ `echo $RUN8F2 | grep -P "^\."` ] && RUN8F2=$OLD_PATH/$RUN8F2;
    RUN8F2=$(cd $(dirname $RUN8F2) && pwd -P)/$(basename $RUN8F2);
    if [ ! -f $RUN8F2 ]; then echo "ERROR: invalid RUN8F2 file option: '$RUN8F2' could not be found"; FOUND_ERROR=1; fi
  else echo "ERROR: Please define all variables for RUN8"; FOUND_ERROR=1; fi
fi
if [[ -n $RUN9F1 || -n $RUN9F2 || -n $RUN9RG ]]; then
  if [[ -n $RUN9F1 && -n $RUN9F2 && -n $RUN9RG ]]; then
    [ `echo $RUN9F1 | grep -P "^\."` ] && RUN9F1=$OLD_PATH/$RUN9F1;
    RUN9F1=$(cd $(dirname $RUN9F1) && pwd -P)/$(basename $RUN9F1);
    if [ ! -f $RUN9F1 ]; then echo "ERROR: invalid RUN9F1 file option: '$RUN9F1' could not be found"; FOUND_ERROR=1; fi
    [ `echo $RUN9F2 | grep -P "^\."` ] && RUN9F2=$OLD_PATH/$RUN9F2;
    RUN9F2=$(cd $(dirname $RUN9F2) && pwd -P)/$(basename $RUN9F2);
    if [ ! -f $RUN9F2 ]; then echo "ERROR: invalid RUN9F2 file option: '$RUN9F2' could not be found"; FOUND_ERROR=1; fi
  else echo "ERROR: Please define all variables for RUN9"; FOUND_ERROR=1; fi
fi
if [[ -n $RUN10F1 || -n $RUN10F2 || -n $RUN10RG ]]; then
  if [[ -n $RUN10F1 && -n $RUN10F2 && -n $RUN10RG ]]; then
    [ `echo $RUN10F1 | grep -P "^\."` ] && RUN10F1=$OLD_PATH/$RUN10F1;
    RUN10F1=$(cd $(dirname $RUN10F1) && pwd -P)/$(basename $RUN10F1);
    if [ ! -f $RUN10F1 ]; then echo "ERROR: invalid RUN10F1 file option: '$RUN10F1' could not be found"; FOUND_ERROR=1; fi
    [ `echo $RUN10F2 | grep -P "^\."` ] && RUN10F2=$OLD_PATH/$RUN10F2;
    RUN10F2=$(cd $(dirname $RUN10F2) && pwd -P)/$(basename $RUN10F2);
    if [ ! -f $RUN10F2 ]; then echo "ERROR: invalid RUN10F2 file option: '$RUN10F2' could not be found"; FOUND_ERROR=1; fi
  else echo "ERROR: Please define all variables for RUN10"; FOUND_ERROR=1; fi
fi

if (( $FOUND_ERROR )); then echo "Please check options and try again"; exit 42; fi

export SAMPLE DIR THREADS MAX_MEM WALLTIME RUN1F1 RUN1F2 RUN1RG RUN2F1 RUN2F2 RUN2RG RUN3F1 RUN3F2 RUN3RG RUN4F1 RUN4F2 RUN4RG RUN5F1 RUN5F2 RUN5RG RUN6F1 RUN6F2 RUN6RG RUN7F1 RUN7F2 RUN7RG RUN8F1 RUN8F2 RUN8RG RUN9F1 RUN9F2 RUN9RG RUN10F1 RUN10F2 RUN10RG VERBOSE SFS KEEP_INTERMEDIATES ADAPTERS CENTER

# ===============================================
# QUEUE THE JOB IN CLUSTER
# ===============================================
cd $DIR; # for job outputs

# ALIGNMENTS
unset MERGE_V37_SORT_SAMS
unset ALIGN_V37_IDS
unset DEPEND;

if [[ -n $RUN1F1 && -n $RUN1F2 && -n $RUN1RG ]]; then
  RUN1V37SAM=${SAMPLE}.${RUN1RG}.paired.aligned.v37.sam
  MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN1V37SAM}"
  if [[ $DISCARD_UNPAIRED_READS == 0 ]]; then
    RUN1V37SAMUNPAIRED=${SAMPLE}.${RUN1RG}.unpaired.aligned.v37.sam;
    MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN1V37SAMUNPAIRED}"
  fi
  align $RUN1F1 $RUN1F2 $RUN1RG $RUN1V37SAM $RUN1V37SAMUNPAIRED
  if [[ ${DEPEND_V37_ALN} ]]; then ALIGN_V37_IDS=${ALIGN_V37_IDS}${DEPEND_V37_ALN}; fi
fi
if [[ -n $RUN2F1 && -n $RUN2F2 && -n $RUN2RG ]]; then
  RUN2V37SAM=${SAMPLE}.${RUN2RG}.paired.aligned.v37.sam
  MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN2V37SAM}"
  if [[ $DISCARD_UNPAIRED_READS == 0 ]]; then
    RUN2V37SAMUNPAIRED=${SAMPLE}.${RUN2RG}.unpaired.aligned.v37.sam;
    MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN2V37SAMUNPAIRED}"
  fi
  align $RUN2F1 $RUN2F2 $RUN2RG $RUN2V37SAM $RUN2V37SAMUNPAIRED
  if [[ ${DEPEND_V37_ALN} ]]; then ALIGN_V37_IDS=${ALIGN_V37_IDS}${DEPEND_V37_ALN}; fi
fi
if [[ -n $RUN3F1 && -n $RUN3F2 && -n $RUN3RG ]]; then
  RUN3V37SAM=${SAMPLE}.${RUN3RG}.paired.aligned.v37.sam
  MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN3V37SAM}"
  if [[ $DISCARD_UNPAIRED_READS == 0 ]]; then
    RUN3V37SAMUNPAIRED=${SAMPLE}.${RUN3RG}.unpaired.aligned.v37.sam;
    MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN3V37SAMUNPAIRED}"
  fi
  align $RUN3F1 $RUN3F2 $RUN3RG $RUN3V37SAM $RUN3V37SAMUNPAIRED
  if [[ ${DEPEND_V37_ALN} ]]; then ALIGN_V37_IDS=${ALIGN_V37_IDS}${DEPEND_V37_ALN}; fi
fi
if [[ -n $RUN4F1 && -n $RUN4F2 && -n $RUN4RG ]]; then
  RUN4V37SAM=${SAMPLE}.${RUN4RG}.paired.aligned.v37.sam
  MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN4V37SAM}"
  if [[ $DISCARD_UNPAIRED_READS == 0 ]]; then
    RUN4V37SAMUNPAIRED=${SAMPLE}.${RUN4RG}.unpaired.aligned.v37.sam;
    MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN4V37SAMUNPAIRED}"
  fi
  align $RUN4F1 $RUN4F2 $RUN4RG $RUN4V37SAM $RUN4V37SAMUNPAIRED
  if [[ ${DEPEND_V37_ALN} ]]; then ALIGN_V37_IDS=${ALIGN_V37_IDS}${DEPEND_V37_ALN}; fi
fi
if [[ -n $RUN5F1 && -n $RUN5F2 && -n $RUN5RG ]]; then
  RUN5V37SAM=${SAMPLE}.${RUN5RG}.paired.aligned.v37.sam
  MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN5V37SAM}"
  if [[ $DISCARD_UNPAIRED_READS == 0 ]]; then
    RUN5V37SAMUNPAIRED=${SAMPLE}.${RUN5RG}.unpaired.aligned.v37.sam;
    MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN5V37SAMUNPAIRED}"
  fi
  align $RUN5F1 $RUN5F2 $RUN5RG $RUN5V37SAM $RUN5V37SAMUNPAIRED
  if [[ ${DEPEND_V37_ALN} ]]; then ALIGN_V37_IDS=${ALIGN_V37_IDS}${DEPEND_V37_ALN}; fi
fi
if [[ -n $RUN6F1 && -n $RUN6F2 && -n $RUN6RG ]]; then
  RUN6V37SAM=${SAMPLE}.${RUN6RG}.paired.aligned.v37.sam
  MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN6V37SAM}"
  if [[ $DISCARD_UNPAIRED_READS == 0 ]]; then
    RUN6V37SAMUNPAIRED=${SAMPLE}.${RUN6RG}.unpaired.aligned.v37.sam;
    MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN6V37SAMUNPAIRED}"
  fi
  align $RUN6F1 $RUN6F2 $RUN6RG $RUN6V37SAM $RUN6V37SAMUNPAIRED
  if [[ ${DEPEND_V37_ALN} ]]; then ALIGN_V37_IDS=${ALIGN_V37_IDS}${DEPEND_V37_ALN}; fi
fi
if [[ -n $RUN7F1 && -n $RUN7F2 && -n $RUN7RG ]]; then
  RUN7V37SAM=${SAMPLE}.${RUN7RG}.paired.aligned.v37.sam
  MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN7V37SAM}"
  if [[ $DISCARD_UNPAIRED_READS == 0 ]]; then
    RUN7V37SAMUNPAIRED=${SAMPLE}.${RUN7RG}.unpaired.aligned.v37.sam;
    MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN7V37SAMUNPAIRED}"
  fi
  align $RUN7F1 $RUN7F2 $RUN7RG $RUN7V37SAM $RUN7V37SAMUNPAIRED
  if [[ ${DEPEND_V37_ALN} ]]; then ALIGN_V37_IDS=${ALIGN_V37_IDS}${DEPEND_V37_ALN}; fi
fi
if [[ -n $RUN8F1 && -n $RUN8F2 && -n $RUN8RG ]]; then
  RUN8V37SAM=${SAMPLE}.${RUN8RG}.paired.aligned.v37.sam
  MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN8V37SAM}"
  if [[ $DISCARD_UNPAIRED_READS == 0 ]]; then
    RUN8V37SAMUNPAIRED=${SAMPLE}.${RUN8RG}.unpaired.aligned.v37.sam;
    MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN8V37SAMUNPAIRED}"
  fi
  align $RUN8F1 $RUN8F2 $RUN8RG $RUN8V37SAM $RUN8V37SAMUNPAIRED
  if [[ ${DEPEND_V37_ALN} ]]; then ALIGN_V37_IDS=${ALIGN_V37_IDS}${DEPEND_V37_ALN}; fi
fi
if [[ -n $RUN9F1 && -n $RUN9F2 && -n $RUN9RG ]]; then
  RUN9V37SAM=${SAMPLE}.${RUN9RG}.paired.aligned.v37.sam
  MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN9V37SAM}"
  if [[ $DISCARD_UNPAIRED_READS == 0 ]]; then
    RUN9V37SAMUNPAIRED=${SAMPLE}.${RUN9RG}.unpaired.aligned.v37.sam;
    MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN9V37SAMUNPAIRED}"
  fi
  align $RUN9F1 $RUN9F2 $RUN9RG $RUN9V37SAM $RUN9V37SAMUNPAIRED
  if [[ ${DEPEND_V37_ALN} ]]; then ALIGN_V37_IDS=${ALIGN_V37_IDS}${DEPEND_V37_ALN}; fi
fi
if [[ -n $RUN10F1 && -n $RUN10F2 && -n $RUN10RG ]]; then
  RUN10V37SAM=${SAMPLE}.${RUN10RG}.paired.aligned.v37.sam
  MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN10V37SAM}"
  if [[ $DISCARD_UNPAIRED_READS == 0 ]]; then
    RUN10V37SAMUNPAIRED=${SAMPLE}.${RUN10RG}.unpaired.aligned.v37.sam;
    MERGE_V37_SORT_SAMS="${MERGE_V37_SORT_SAMS}@INPUT=${RUN10V37SAMUNPAIRED}"
  fi
  align $RUN10F1 $RUN10F2 $RUN10RG $RUN10V37SAM $RUN10V37SAMUNPAIRED
  if [[ ${DEPEND_V37_ALN} ]]; then ALIGN_V37_IDS=${ALIGN_V37_IDS}${DEPEND_V37_ALN}; fi
fi

echo "*** Post Alignment"

# MERGE SAM FILES
NAME=pipeline.generic.part_2_0.combine_sample_libs
DONE=.${SAMPLE}.${NAME}.v37.done
MERGE_V37_SORT_BAM=${SAMPLE}.aligned.v37.merged.sorted.bam
if [[ ! -z $ALIGN_V37_IDS ]]; then DEPEND="-W "$dependMod"afterok${ALIGN_V37_IDS}"; else unset DEPEND; fi
if [[ ! `find ./ -name ${DONE}` || $SFS -eq 1 ]]; then
  v37_merge_sample_libs=`$QUEUE $ACCOUNT \
    -l nodes=1:ppn=$THREADS \
    -l mem=${MAX_MEM}g \
    -l walltime=$WALLTIME \
    -N ${SAMPLE}.${NAME}.v37 \
    $DEPEND \
    -V \
    -v MERGE_SORT_SAMS=${MERGE_V37_SORT_SAMS},MERGE_SORT_BAM=${MERGE_V37_SORT_BAM},DONE=${DONE} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  v37_merge_sample_libs=$(echo $v37_merge_sample_libs)
  echo "[Q] Merge Sort v37  : $v37_merge_sample_libs [$DEPEND]"
  V37_MERGE_DEPEND="-W "$dependMod"afterok:"$v37_merge_sample_libs;
else
  echo "[D] Merge Sort v37"
fi

# GATK WGS PRE-PROCESSING
# 2013-10-11: dan changed walltime for this function to 1 week (hard coded); this is because the job too often exceeds the default 2-day limit
NAME=pipeline.wgs.part_2_1.GATK_data_processing.in_bam
DONE=.${SAMPLE}.${NAME}.v37.done
if [[ ! `find ./ -name ${DONE}` || $SFS -eq 1 ]]; then
  gatk_v37=`$QUEUE $ACCOUNT  \
    -l nodes=1:ppn=$THREADS \
    -l mem=${MAX_MEM}g \
    -l walltime=168:0:0 \
    -N ${SAMPLE}.${NAME}.v37 \
    $V37_MERGE_DEPEND \
    -V \
    -v ALIGNED_BAM=${MERGE_V37_SORT_BAM},USE_INDELS=1,GENOME=v37,DONE=${DONE} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  gatk_v37=$(echo $gatk_v37)
  echo "[Q] GATK WGS        : $gatk_v37 [$V37_MERGE_DEPEND]"
  GATK_WGS_DEPEND="-W "$dependMod"afterok:"$gatk_v37;
else echo "[D] GATK WGS"
fi

## STATS AND QC
NAME=pipeline.wgs.part_3_1.WGS_stats_QC
DONE=.${SAMPLE}.${NAME}.done
if [[ ! `find ./ -name ${DONE}` || $SFS -eq 1 ]]; then
  WGS_sample_stats_QC=`$QUEUE $ACCOUNT  \
    -l nodes=1:ppn=$THREADS \
    -l mem=${MAX_MEM}g \
    -l walltime=$WALLTIME \
    -N ${SAMPLE}.${NAME}\
    $GATK_WGS_DEPEND \
    -V \
    -v DONE=${DONE} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  WGS_sample_stats_QC=$(echo $WGS_sample_stats_QC)
  echo "[Q] Stats and QC    : $WGS_sample_stats_QC [$GATK_WGS_DEPEND]"
  STATS_DEPEND="-W "$dependMod"afterok:"$WGS_sample_stats_QC;
else echo "[D] Stats and QC"
fi

## VARIANT CALLING
NAME=pipeline.wgs.part_3_2.GATK_UG
DONE=.${SAMPLE}.${NAME}.done
if [[ ! `find ./ -name ${DONE}` || $SFS -eq 1 ]]; then
  WGS_sample_variant_calling=`$QUEUE $ACCOUNT  \
    -l nodes=1:ppn=$THREADS \
    -l mem=${MAX_MEM}g \
    -l walltime=$WALLTIME \
    -N ${SAMPLE}.${NAME}\
    $GATK_WGS_DEPEND \
    -V \
    -v DONE=${DONE} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  WGS_sample_variant_calling=$(echo $WGS_sample_variant_calling)
  echo "[Q] Variant calling    : $WGS_sample_variant_calling [$GATK_WGS_DEPEND]"
  CALLING_DEPEND="-W "$dependMod"afterok:"$WGS_sample_variant_calling;
else echo "[D] Variant calling"
fi

# CLEANUP STAGE
NAME=pipeline.wgs.part_4_1.cleanup
DONE=.${SAMPLE}.${NAME}.done
if [[ ! `find ./ -name ${DONE}` || $SFS -eq 1 ]]; then
  run_cleanup=`$QUEUE $ACCOUNT  \
    -l nodes=1:ppn=$THREADS \
    -l mem=${MAX_MEM}g \
    -l walltime=1:0:0 \
    -N ${SAMPLE}.${NAME}\
    $STATS_DEPEND \
    -V \
    -v DONE=${DONE} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  run_cleanup=$(echo $run_cleanup)
  echo "[Q] Run Cleanup     : $run_cleanup [$STATS_DEPEND]"
else echo "[D] Run Cleanup"
fi

cd $OLD_PATH
exit 0


