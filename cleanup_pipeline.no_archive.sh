#!/bin/bash

rm -rf temp 
rm -v *trimmomaticPE.*paired.gz input*unpaired.fastq.gz .*.done *.TMP_MERGE.ba? *clean.ba? *clean.dedup.ba? *.intervals *.pre_recal.table *paired.aligned.v37.sam *.merged.sorted.bam* MT.*.clean.dedup.recal.ba? *.aligned.MT.merged.sorted.bai *.sam 
find -perm 600 -exec chmod 664 \{\} \;;
mkdir LOGS;
mv *.o[0-9]* *.out *.stderr *DataProcessingPipeline.*.txt launch_record* launch.txt *jobreport.pdf LOGS;
mkdir STATS;
mv TEQCreport DepthOfCoverage.* *.metrics *.flagstat *.vs.*.csv *.verifyBamID.* fastqc.* *.idxstats STATS;

