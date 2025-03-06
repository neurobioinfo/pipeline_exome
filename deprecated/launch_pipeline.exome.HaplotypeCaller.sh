#!/bin/bash

#hard-coded version number. Don't forget to update!!
VERSION=1.1.1; echo "exome analysis pipeline version $VERSION"
#last updated version 2014-06-12

# ===============================================
# default variables values
# ===============================================
unset SAMPLE DIR THREADS MAX_MEM PIPELINE_HOME WALLTIME EXTRA_CONF CAPTURE_KIT RUN1F1 RUN1F2 RUN1LIB RUN1RUN RUN2F1 RUN2F2 RUN2LIB RUN2RUN RUN3F1 RUN3F2 RUN3LIB RUN3RUN RUN4F1 RUN4F2 RUN4LIB RUN4RUN RUN5F1 RUN5F2 RUN5LIB RUN5RUN RUN6F1 RUN6F2 RUN6LIB RUN6RUN RUN7F1 RUN7F2 RUN7LIB RUN7RUN RUN8F1 RUN8F2 RUN8LIB RUN8RUN RUN9F1 RUN9F2 RUN9LIB RUN9RUN QUEUE ACCOUNT VERBOSE CENTER SFS KEEP_INTERMEDIATES 

#Assuming script is in root of PIPELINE_HOME
export PIPELINE_HOME=$(cd $(dirname $0) && pwd -P)
export SFS=0
export KEEP_INTERMEDIATES=0
export DISCARD_UNPAIRED_READS=0
export EXPECTED_DONE_FILES='expected.done.files.txt'

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
# PARSING ARGUMENTS
# ===============================================

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt --name pipeline --alternative --unquoted --options hs:d:t:m:vw:c:a:x: --longoptions sample:,dir:,threads:,max_mem:,walltime:,capture_kit:,mode:,bed:,bed_flank:,sites:,adapters:,r1b:,run1bam:,r2b:,run2bam:,r3b:,run3bam:,r4b:,run4bam:,r5b:,run5bam:,r6b:,run6bam:,r7b:,run7bam:,r8b:,run8bam:,r9b:,run9bam:,r1f1:,run1f1:,r1f2:,run1f2:,r1r:,run1run:,r1l:,run1lib:,r2f1:,run2f1:,r2f2:,run2f2:,r2r:,run2run:,r2l:,run2lib:,r3f1:,run3f1:,r3f2:,run3f2:,r3r:,run3run:,r3l:,run3lib:,r4f1:,run4f1:,r4f2:,run4f2:,r4r:,run4run:,r4l:,run4lib:,r5f1:,run5f1:,r5f2:,run5f2:,r5r:,run5run:,r5l:,run5lib:,r6f1:,run6f1:,r6f2:,run6f2:,r6r:,run6run:,r6l:,run6lib:,r7f1:,run7f1:,r7f2:,run7f2:,r7r:,run7run:,r7l:,run7lib:,r8f1:,run8f1:,r8f2:,run8f2:,r8r:,run8run:,r8l:,run8lib:,r9f1:,run9f1:,r9f2:,run9f2:,r9r:,run9run:,r9l:,run9lib:,extra:,center:,account:,keep_intermediates,discard_unpaired,start_from_scratch,verbose,help -- "$@")
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
    --sites) CALLING_BED="$2" ; shift ;; 
    -x| --extra) echo -n "" ; shift ;;
    -a| --adapters) ADAPTERS="$2" ; shift ;;
    --account) ACCOUNT="$2" ; shift ;;
    -v| --verbose) VERBOSE=1 ;; 
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
if [[ -n $ACCOUNT ]]; then ACCOUNT="-A $ACCOUNT"; fi
OLD_PATH=`pwd -P`
[ `echo $DIR | grep -P "^\."` ] && DIR=$OLD_PATH/$DIR; 
DIR=$(cd $DIR 2>/dev/null && pwd -P);
if [ ! -d $DIR ]; then echo "ERROR: wrong mandatory option: -d (run directory $DIR) is invalid"; FOUND_ERROR=1; fi

export SAMPLE DIR THREADS MAX_MEM WALLTIME CAPTURE_KIT BED BED_FLANKS RUN1BAM1 RUN1F1 RUN1F2 RUN1RUN RUN1LIB RUN2BAM RUN2F1 RUN2F2 RUN2RUN RUN2LIB RUN3BAM RUN3F1 RUN3F2 RUN3RUN RUN3LIB RUN4BAM RUN4F1 RUN4F2 RUN4RUN RUN4LIB RUN5BAM RUN5F1 RUN5F2 RUN5RUN RUN5LIB RUN6BAM RUN6F1 RUN6F2 RUN6RUN RUN6LIB RUN7BAM RUN7F1 RUN7F2 RUN7RUN RUN7LIB RUN8BAM RUN8F1 RUN8F2 RUN8RUN RUN8LIB RUN9BAM RUN9F1 RUN9F2 RUN9RUN RUN9LIB VERBOSE SFS KEEP_INTERMEDIATES ADAPTERS CENTER

# ===============================================
# QUEUE THE JOB IN CLUSTER
# ===============================================
cd $DIR; # for job outputs

# ALIGNMENTS
unset MERGE_V37_SORT_SAMS
unset MERGE_MT_SORT_SAMS
unset ALIGN_V37_IDS
unset ALIGN_MT_IDS

unset DEPEND;

## GATK UG total recall 
NAME=pipeline.exome.part_3_6.GATK_HC
DONE=.${SAMPLE}.${NAME}.done

if [[ ! -z $CALLING_BED ]]; then CALLING_BED_ARG="BEDFILE=$CALLING_BED,"; else unset CALLING_BED_ARG; fi

if [[ ! `find ./ -maxdepth 1 -name ${DONE}` || $SFS -eq 1 ]]; then
  gatk_HC=`$QUEUE $ACCOUNT  \
    -l nodes=1:ppn=$THREADS \
    -l mem=${MAX_MEM}g \
    -l walltime=48:0:0 \
    -N ${SAMPLE}.${NAME}\
    $GATK_EXOME_DEPEND \
    -V \
    -v ${CALLING_BED_ARG}DONE=${DONE} \
    ${PIPELINE_HOME}/data/qsub/${NAME}.qsub`
  gatk_HC=$(echo $gatk_HC)
  echo "[Q] GATK HC : $gatk_HC"
else echo "[D] GATK HC"
fi

cd $OLD_PATH
exit 0
