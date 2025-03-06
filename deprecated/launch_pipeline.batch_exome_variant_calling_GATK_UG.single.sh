#!/bin/bash

#define function to handle error messages
Usage() {
        echo; echo -e "usage: sh launch.txt [optional arguments (-cbd)] OR -h\n"
	echo -e "\toptional arguments:\n\t\t-c = capture kit (38Mbp, 50Mbp, V4, V5, V5_Clinic, TruSeq or NGv3)\n\t\t-b = batch identifier\n\t\t-d = date (in format yyyy-mm-dd)\n\t\t-w = walltime (in format h:m:s)\n\t\t-g = GAIIx sequencing (triggers GAIIx pipeline instead of default HiSeq pipeline)\n\t\t-P sample from Pharmacogenomics centre (triggers special variables for this unique case)\n\t\t-A = accelerated run; will launch all chromosomes simultaneously. DO NOT USE IF RUNNING >2 HiSeq FREEZE RUNS\n\t\t-F = Freeze name (use to override default behaviour)\n\t\t-D = run directory (use to override default behaviour)\n\t\t-L use LD_PRELOAD during java executions\n\t\t-M modified PBS arguments (put in double quotes)"
}

unset CAPTURE BATCH DATE GAIIx
WALLTIME=12:00:00

#process options and and pass arguments to script variables (where applicable)
while getopts ":hc:b:d:w:gPALF:D:M:" opt; do
case $opt in
h ) Usage; exit 1;;
c ) CAPTURE=$OPTARG;;
b ) BATCH=$OPTARG;;
d ) DATE=$OPTARG;;
g ) GAIIx=1;;
w ) WALLTIME=$OPTARG;;
P ) PGX=1;;
A ) ACCELERATE=1;;
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

export Q=$PIPELINE_HOME/soft/src/queueInterpreter/QueueInterpreter
export W="-l walltime=$WALLTIME"

normal() {
  CHR1=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.1.command.qsub -N freeze.$FREEZE.part1.chr.1 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
  CHR2=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.2.command.qsub -N freeze.$FREEZE.part1.chr.2 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
  CHR3=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.3.command.qsub -N freeze.$FREEZE.part1.chr.3 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
  CHR4=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.4.command.qsub -N freeze.$FREEZE.part1.chr.4 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR1 $MODIFIED_PBS_ARGS`
  CHR5=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.5.command.qsub -N freeze.$FREEZE.part1.chr.5 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR2 $MODIFIED_PBS_ARGS`
  CHR6=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.6.command.qsub -N freeze.$FREEZE.part1.chr.6 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR3 $MODIFIED_PBS_ARGS`
  CHR7=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.7.command.qsub -N freeze.$FREEZE.part1.chr.7 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR4 $MODIFIED_PBS_ARGS`
  CHR8=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.8.command.qsub -N freeze.$FREEZE.part1.chr.8 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR5 $MODIFIED_PBS_ARGS`
  CHR9=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.9.command.qsub -N freeze.$FREEZE.part1.chr.9 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR6 $MODIFIED_PBS_ARGS`
  CHR10=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.10.command.qsub -N freeze.$FREEZE.part1.chr.10 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR7 $MODIFIED_PBS_ARGS`
  CHR11=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.11.command.qsub -N freeze.$FREEZE.part1.chr.11 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR8 $MODIFIED_PBS_ARGS`
  CHR12=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.12.command.qsub -N freeze.$FREEZE.part1.chr.12 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR9 $MODIFIED_PBS_ARGS`
  CHR13=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.13.command.qsub -N freeze.$FREEZE.part1.chr.13 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR10 $MODIFIED_PBS_ARGS`
  CHR14=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.14.command.qsub -N freeze.$FREEZE.part1.chr.14 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR11 $MODIFIED_PBS_ARGS`
  CHR15=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.15.command.qsub -N freeze.$FREEZE.part1.chr.15 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR12 $MODIFIED_PBS_ARGS`
  CHR16=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.16.command.qsub -N freeze.$FREEZE.part1.chr.16 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR13 $MODIFIED_PBS_ARGS`
  CHR17=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.17.command.qsub -N freeze.$FREEZE.part1.chr.17 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR14 $MODIFIED_PBS_ARGS`
  CHR18=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.18.command.qsub -N freeze.$FREEZE.part1.chr.18 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR15 $MODIFIED_PBS_ARGS`
  CHR19=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.19.command.qsub -N freeze.$FREEZE.part1.chr.19 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR16 $MODIFIED_PBS_ARGS`
  CHR20=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.20.command.qsub -N freeze.$FREEZE.part1.chr.20 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR17 $MODIFIED_PBS_ARGS`
  CHR21=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.21.command.qsub -N freeze.$FREEZE.part1.chr.21 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR18 $MODIFIED_PBS_ARGS`
  CHR22=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.22.command.qsub -N freeze.$FREEZE.part1.chr.22 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR19 $MODIFIED_PBS_ARGS`
  CHRX=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.X.command.qsub -N freeze.$FREEZE.part1.chr.X -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR20 $MODIFIED_PBS_ARGS`
  CHRN=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.the_rest.command.qsub -N freeze.$FREEZE.part1.chr.the_rest -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR21 $MODIFIED_PBS_ARGS`
  PART2=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part2.command.qsub -N freeze.$FREEZE.part2 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok:$CHR22:$CHRX:$CHRN $MODIFIED_PBS_ARGS`
  for i in $CHR{1..22} $CHRX $CHRN $PART2; do echo $i; done
}

accelerated() {
#  CHR1=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.1.command.qsub -N freeze.$FREEZE.part1.chr.1 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR2=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.2.command.qsub -N freeze.$FREEZE.part1.chr.2 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR3=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.3.command.qsub -N freeze.$FREEZE.part1.chr.3 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR4=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.4.command.qsub -N freeze.$FREEZE.part1.chr.4 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR5=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.5.command.qsub -N freeze.$FREEZE.part1.chr.5 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR6=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.6.command.qsub -N freeze.$FREEZE.part1.chr.6 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR7=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.7.command.qsub -N freeze.$FREEZE.part1.chr.7 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR8=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.8.command.qsub -N freeze.$FREEZE.part1.chr.8 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR9=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.9.command.qsub -N freeze.$FREEZE.part1.chr.9 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR10=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.10.command.qsub -N freeze.$FREEZE.part1.chr.10 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR11=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.11.command.qsub -N freeze.$FREEZE.part1.chr.11 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR12=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.12.command.qsub -N freeze.$FREEZE.part1.chr.12 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR13=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.13.command.qsub -N freeze.$FREEZE.part1.chr.13 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR14=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.14.command.qsub -N freeze.$FREEZE.part1.chr.14 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR15=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.15.command.qsub -N freeze.$FREEZE.part1.chr.15 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR16=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.16.command.qsub -N freeze.$FREEZE.part1.chr.16 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR17=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.17.command.qsub -N freeze.$FREEZE.part1.chr.17 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR18=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.18.command.qsub -N freeze.$FREEZE.part1.chr.18 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR19=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.19.command.qsub -N freeze.$FREEZE.part1.chr.19 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR20=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.20.command.qsub -N freeze.$FREEZE.part1.chr.20 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR21=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.21.command.qsub -N freeze.$FREEZE.part1.chr.21 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHR22=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.22.command.qsub -N freeze.$FREEZE.part1.chr.22 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHRX=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.X.command.qsub -N freeze.$FREEZE.part1.chr.X -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  CHRN=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.chr.the_rest.command.qsub -N freeze.$FREEZE.part1.chr.the_rest -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN $MODIFIED_PBS_ARGS`
#  PART2=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part2.command.qsub -N freeze.$FREEZE.part2 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok\`for i in $CHR{1..22} $CHRX $CHRN; do printf ":%s" $i; done\` $MODIFIED_PBS_ARGS`
#  for i in $CHR{1..22} $CHRX $CHRN $PART2; do echo $i; done  
  for i in {1..22} X N; do
    if [[ ! -e $DATE.$FREEZE.chr.${i}.vcf.done ]]; then
        CURRENTCHR=`$Q $W $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part1.command.qsub -N freeze.$FREEZE.part1.chr.${i} -V -v CHROMOSOME=$i,DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN,USELD=$USELD $MODIFIED_PBS_ARGS`
        AFTEROK=$AFTEROK:$CURRENTCHR
        echo -e "$CURRENTCHR\tchr${i}"
    fi
  done
  PART2=`$Q $DIR_SCRIPT/pipeline.batch_exome_variant_calling_GATK_UG.HiSeq_part2.command.qsub -N freeze.$FREEZE.part2 -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN -W depend=afterok${AFTEROK} $MODIFIED_PBS_ARGS`
  echo -e "$PART2\tcombination/calibration/annotation"
}

#launch the qsub jobs
cd $DIR_RUN
if [[ $GAIIx == 1 ]]; then \
  QSUB=`$Q $DIR_SCRIPT/GAIIx.GATK.command.qsub -N freeze.$FREEZE -V -v DATE=$DATE,FREEZE=$FREEZE,DIR=$DIR_RUN`
  echo $QSUB
else \
  if [[ $ACCELERATE ]]; then accelerated; else normal; fi
fi
