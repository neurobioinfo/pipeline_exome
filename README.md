Quick howto:

***
first set your own $PIPELINE_HOME variable as the base of this package i.e. where the package has been checked out
***

0) launch svn propset svn:ignore -F .svnignore . immediately after checking out the svn project (keep this file up to date with list of files to ignore)
1) Recompile bwa (2x) in $PIPELINE_HOME/soft/src/bwa* 
2) Recompile samtools in $PIPELINE_HOME/soft/src/samtools*
	cd $PIPELINE_HOME/soft/src/htslib
	./configure --prefix=$PIPELINE_HOME/soft/packages/htslib
	make
	make install
	cd $PIPELINE_HOME/soft/src/samtools
	./configure --prefix=$PIPELINE_HOME/soft/packages/samtools
	make
	make install
3) Recompile tabix in $PIPELINE_HOME/soft/src/tabix
4) Recompile verifyBAMId in $PIPELINE_HOME/soft/packages/verifyBamID/
	cd verifyBamID
	make cloneLib (in the case ../libStatGen does not exist)
	make
5) Recompile parsePileup in $PIPELINE_HOME/soft/src/statsFromRun	# do we still need this?
6) Recompile queueInterpreter in $PIPELINE_HOME/soft/src/queueInterpreter
7) Recompile iobuff in $PIPELINE_HOME/soft/src/iobuff
8) Check perl dependencies

if needed:
perl -MCPAN -e shell
o conf makepl_arg PREFIX='$PIPELINE_HOME/soft/packages/perl'
o conf commit
install <missing module>

9) manually create link to reference files and annovar databases (references: to $PIPELINE_HOME/data/reference; annovar dbs: to $PIPELINE_HOME/soft/packages/VarAnnot/humandb)
   (master copies of both directories exist here: https://verlaine.rqchp.qc.ca/pipeline/data/)
   (on briaree: folders can be found here: /RQexec/dionnela/data/pipeline.external_data)

 ex: >/pipeline.svn/data]$ ln -s /RQexec/GROUP/rouleau/COMMUN/data/pipeline.external_data/reference
     >/pipeline.svn/data]$ ln -s /RQexec/GROUP/rouleau/COMMUN/data/pipeline.external_data/varannot.humandb

*****
How to Launch pipeline:

export PIPELINE_HOME=/RQexec/dionnela/data/pipeline.svn
. $PIPELINE_HOME/data/templates/init_pipeline.sh $PIPELINE_HOME
. $PIPELINE_HOME/data/templates/function.library 
sample_init BA000VIL0001.test ./

sample_GATK_data_processing ...

>Example for exome, in the GATK dir:
/RQexec/dionnela/data/pipeline.svn/launch_pipeline.exome.fastqs.splits.sh -s S14268_A -d $PWD -r1f1 ../RAW/S14268_A_1.2012-12-06.fastq.gz -r1f2 ../RAW/S14268_A_2.2012-12-06.fastq.gz -r1rg run1 -r2f1 ../RAW/S14268_A_1.2013-01-15.fastq.gz -r2f2 ../RAW/S14268_A_2.2013-01-15.fastq.gz -r2rg run2 -x /RQexec/dionnela/data/pipeline.svn/data/templates/init_pipeline.briaree.exome.V4.pgx.sh --verbose > launch.out

*****
To verify that everything went well:
for sample in S*; do cd $sample/I*/G*; for file in $(cat expected.done.files.txt); do find $file; done; cd ../../../; done > /dev/null
