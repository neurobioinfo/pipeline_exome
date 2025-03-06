#!/bin/bash

#define function to handle error messages
Usage() {
        echo; echo -e "usage: sh launch.txt [optional arguments (-cbd)] OR -h\n"
	echo -e "\toptional arguments:\n\t\t-c = capture kit (38Mbp, 50Mbp, V4, V5, V5_Clinic, TruSeq or NGv3)\n\t\t-b = batch identifier\n\t\t-d = date (in format yyyy-mm-dd)\n\t\t-w = walltime (in format h:m:s)\n\t\t-g = GAIIx sequencing (triggers GAIIx pipeline instead of default HiSeq pipeline)\n\t\t-P sample from Pharmacogenomics centre (triggers special variables for this unique case)\n\t\t-A = pbs account for job submissions (default = ish-284-ad)\n\t\t-F = Freeze name (use to override default behaviour)\n\t\t-D = run directory (use to override default behaviour)\n\t\t-L use LD_PRELOAD during java executions\n\t\t-M modified PBS arguments (put in double quotes)"
}

unset CAPTURE BATCH DATE GAIIx
WALLTIME=12:00:00
ACCOUNT=ish-284-ad

#process options and and pass arguments to script variables (where applicable)
while getopts ":hc:b:d:w:gPA:LF:D:M:" opt; do
case $opt in
h ) Usage; exit 1;;
c ) CAPTURE=$OPTARG;;
b ) BATCH=$OPTARG;;
d ) DATE=$OPTARG;;
g ) GAIIx=1;;
w ) WALLTIME=$OPTARG;;
A ) ACCOUNT=$OPTARG;;
P ) PGX=1;;
F ) FREEZE=$OPTARG;;
D ) DIR_RUN=$OPTARG;;
M ) MODIFIED_PBS_ARGS=$OPTARG;;
L ) USELD=1;;
\?) echo "ERROR: Invalid option: -$OPTARG"; exit 1;;
: ) echo "ERROR: Option -$OPTARG requires an argument"; exit 1;;
esac
done
#check for missing variables and nonsensical combinations; set other script variables
export PIPELINE_HOME=$(cd $(dirname $0) && pwd -P)
#init the pipeline vars (GATK_DIR)
. $PIPELINE_HOME/data/templates/init_pipeline.sh $PIPELINE_HOME

echo $DIR_RUN

DIR_SCRIPT=$PIPELINE_HOME/data/qsub
#if [[ -z $CAPTURE ]]; then echo "ERROR: please provide a capture kit (with -c)"; exit 1; fi
#if [[ $CAPTURE != 38Mbp && $CAPTURE != 50Mbp && $CAPTURE != V4 && $CAPTURE != V5 && $CAPTURE != V5_Clinic && $CAPTURE != TruSeq && $CAPTURE != NGv3 ]]; then echo "ERROR: capture must be one of \"38Mbp\" \"50Mbp\" \"V4\" \"V5\" \"TruSeq\" or \"NGv3\""; exit 42; fi
if [[ -z $DATE    ]]; then echo "date will be set by default to today's date"; DATE=`date +%F`; fi
if [[ $PGX == 1 ]]; then CAPTURE=PGX.V4-XT2; fi
if [[ -z $FREEZE ]]; then if [[ -n $BATCH   ]]; then FREEZE=HiSeq.$CAPTURE.$BATCH; else FREEZE=HiSeq.$CAPTURE; fi; fi
if [[ -z $DIR_RUN ]]; then DIR_RUN=/alloc/grouleau/COMMUN/pipeline_multisample_analyses/variants.FREEZE/GATK/Illumina/HiSeq.$CAPTURE/next.run
  if [[ $GAIIx == 1 ]]; then FREEZE=GAIIx.38Mbp; DIR_RUN=/alloc/grouleau/COMMUN/pipeline_multisample_analyses/variants.FREEZE/GATK/Illumina/$FREEZE/next.run; fi
  if [[ $GAIIx == 1 && (-n $BATCH || $CAPTURE != "38Mbp") ]]; then echo "incorrect parameters; GAIIx (flagged by -g) can only have capture kit \"38Mbp\""; exit 42; fi
fi
if [ ! -d $DIR_RUN ]; then echo "run directory $DIR_RUN does not exist"; exit 42; fi
if [[ $GAIIx && $PGX ]]; then echo "GAIIx and Pharmacogenomics triggers cannot both be selected"; exit 42; fi
if [[ $USELD == 1 ]]; then echo "LD_PRELOAD will be used during java executions"; else USELD=0; fi
#if [[ -z $FREEZE ]]; then if [[ $PGX == 1 ]]; then FREEZE=HiSeq.PGX.V4-XT2.$BATCH; DIR_RUN=/alloc/grouleau/COMMUN/pipeline_multisample_analyses/variants.FREEZE/GATK/Illumina/HiSeq.PGX.V4-XT2/next.run; fi; fi

#quick check to make sure bam list exists in the right place
BAMS=$FREEZE.bam.list
if [ ! -s $DIR_RUN/$BAMS ]; then echo "ERROR: make sure bam file list \"$BAMS\" exists in run directory $DIR_RUN (and is not empty)"; exit 42; fi

export Q="qsub -A $ACCOUNT"
export W="-l walltime=$WALLTIME"

launch() {
  for i in {1..22} X N; do
    if [[ ! -e $DATE.$FREEZE.chr.${i}.vcf.done ]]; then
        CURRENTCHR=`$Q $W $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.command.qsub -N freeze.$FREEZE.part1.chr.${i} -V -v CHROMOSOME=$i,DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN,USELD=$USELD -a $(date --date '3 minute' +%H%M) $MODIFIED_PBS_ARGS`
        AFTEROK=$AFTEROK:$CURRENTCHR
        echo -e "$CURRENTCHR\tchr${i}"
    fi
  done
  if [[ $AFTEROK ]]; then DEP="-W depend=afterok${AFTEROK}"; fi
  PART2=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part2.command.qsub -N freeze.$FREEZE.part2 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $DEP -a $(date --date '3 minute' +%H%M) $MODIFIED_PBS_ARGS`
  echo $Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part2.command.qsub -N freeze.$FREEZE.part2 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $DEP -a $(date --date '3 minute' +%H%M) $MODIFIED_PBS_ARGS
  echo -e "$PART2\tcombination/calibration/annotation"
}

#launch the qsub jobs
cd $DIR_RUN
if [[ $GAIIx == 1 ]]; then 
  QSUB=`$Q $DIR_SCRIPT/GAIIx.GATK.command.qsub -N freeze.$FREEZE -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -a $(date --date '3 minute' +%H%M)`
  echo $QSUB
else 
  launch
fi

echo "check using \"for suffix in \$\(cat \$PIPELINE_HOME/launch_pipeline.batch_exome_variant_calling_GATK_UG.check_file); do find \*\$\{FREEZE\}.\$\{suffix\}; done\""

