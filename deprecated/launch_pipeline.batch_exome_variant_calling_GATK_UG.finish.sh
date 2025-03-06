#in CURRENT FREEZE dir
for file in *.RESULTS.vcf; do echo $file; perl -ne 'if ($_ =~ /^#/ ) {next}; chomp; @data=split("\t",$_); if ( $previous_chr && $previous_chr !~ /$data[0]/ ) {$ref_chr=`grep -P \"^$previous_chr\t\" /RQexec/dionnela/data/pipeline.svn/data/reference/human_g1k_v37.fasta.fai | cut -f2`; chomp $ref_chr; $ratio=$previous_pos/$ref_chr; print "$previous_chr $previous_pos $ref_chr $ratio\n";} $previous_chr=$data[0]; $previous_pos=$data[1];' < $file; done

#in specific freeze dir:
FREEZE=RMGA.Various.batch1; for suffix in $(cat /RQexec/dionnela/data/pipeline.svn/launch_pipeline.batch_exome_variant_calling_GATK_UG.check_file); do find *${FREEZE}.${suffix} ; done
