#!/bin/bash

#do to before:
#for sample in S*; do cd $sample/Ill*/GATK*; for file in $(cat expected.done.files.txt); do find $file; done; cd ../../../; done

#after:
##NOT FROM briaree3 tf-briaree!!!!
##
#for sample in S*; do cd $sample/Ill*/GATK*; /RQexec/dionnela/data/pipeline.svn/cleanup_pipeline.sh; cd ../../../; don

rm -rf temp
rm -v *trimmomaticPE.*paired.gz .*.done *.TMP_MERGE.bai
find -perm 600 -exec chmod 664 \{\} \;;
mkdir LOGS;
mv *.o[0-9]* *.out *.stderr *DataProcessingPipeline.*.txt launch_record* launch.txt *jobreport.pdf LOGS;
mkdir STATS;
mv TEQCreport DepthOfCoverage.* *.metrics *.idxstats *.flagstat *.vs.*.csv *.verifyBamID.* fastqc.* STATS;
tar czvf SMALL_FILES.tar.gz $(echo  $(find . -maxdepth 1 -type f -size -1024k)    $(find -name "MT.*.ba*")    $(find -name "*non_pri_aln.ba*")    LOGS STATS |    sed -e 's/ /\n/g' -e '/gz$/d' | sort | uniq);
rm -rfv $(echo    $(find . -maxdepth 1 -type f -size -1024k)    $(find -name "MT.*.ba*")    $(find -name "*non_pri_aln.ba*")    LOGS STATS temp |    sed -e 's/ /\n/g' | sed -e '/gz$/d' | sort | uniq);



