#!/usr/bin/perl -w

use strict;
use warnings;

use Pod::Usage;
use Getopt::Long;

use File::Basename;
use File::Find;

our $VERSION = '$Revision: 1.4 $';
our $LAST_CHANGED_DATE = '$LastChangedDate: 2013-05-12 11:41 $';

our ($help, $man);
our ($flagstat, $rawBAM, $rawFASTQ, $rawFASTQC, $output);
#our ($logFile, $time);
#our ($rawDirectory, $runDirectory, $metricsFile);
our ($raw_reads, $filtered_reads, $mapped_reads, $percent_mapped, $paired_mapped_reads, $non_dup_reads, $dup_reads);
#our ($raw, $bwa, $mapped, $paired, $percentMapped, $nonDupMapped, $percentDup);
#our ($RAWFILE, @RUNFILES, @METFILES);
#our $SAMPLE;

processArguments();
parseMetrics();
printMetrics();

sub processArguments {
    my @command_line = @ARGV;       #command line argument
    GetOptions('help|h'=>\$help, 'man|m'=>\$man,'flagstat=s'=>\$flagstat, 'raw_fastqc=s'=>\$rawFASTQC, 'raw_bam=s'=>\$rawBAM, 'raw_fastq=s'=>\$rawFASTQ, 'output=s'=>\$output) or pod2usage ();

    $help and pod2usage (-verbose=>1, -exitval=>1, -output=>\*STDOUT);
    $man and pod2usage (-verbose=>2, -exitval=>1, -output=>\*STDOUT);
    @ARGV == 0 or pod2usage ("Syntax error\n" . join " ", @command_line);

    ($flagstat && -e $flagstat) or pod2usage("Error: Please specify the flagstat file...\n" . join " ", @command_line);
    if ($rawBAM) { ( -e $rawBAM) or pod2usage("Error: Invalid RAW BAM file...\n" . join " ", @command_line); }
    if ($rawFASTQ) { (-e $rawFASTQ) or pod2usage("Error: Invalid RAW FASTQ file...\n" . join " ", @command_line); }
    if ( !$rawBAM && !$rawFASTQ && !$rawFASTQC ) { pod2usage("Error: Provide either raw_bam, raw_fastq or fastqc ...\n" . join " ", @command_line); }
}

sub parseMetrics {
  parseFlagstat($flagstat);
  if ( $rawFASTQ ) { $raw_reads=rawReadsFromFastq($rawFASTQ); }
  elsif ( $rawFASTQC ) { $raw_reads=rawReadsFromFastqc($rawFASTQC); } 
  else { $raw_reads=rawReadsFromBam($rawBAM); }

}

sub parseFlagstat {
    my $file=shift;
    open(FLAGSTAT, "<$file") || die "Could not open file $file\n$!\n";
#202646138 + 0 in total (QC-passed reads + QC-failed reads)
#46881591 + 0 duplicates
#201781993 + 0 mapped (99.57%:-nan%)
#202646138 + 0 paired in sequencing
#101323069 + 0 read1
#101323069 + 0 read2
#201122846 + 0 properly paired (99.25%:-nan%)
#201638877 + 0 with itself and mate mapped
#143116 + 0 singletons (0.07%:-nan%)
#378027 + 0 with mate mapped to a different chr
#322084 + 0 with mate mapped to a different chr (mapQ>=5)

    my $line;
    $line = <FLAGSTAT>;
    $line =~ /(\d+) \+/;
    $filtered_reads = $1;

    $line = <FLAGSTAT>;
    $line =~ /(\d+) \+/;
    $dup_reads = $1;

    $line = <FLAGSTAT>;
    $line =~ /(\d+) \+/;
    $mapped_reads = $1;

    <FLAGSTAT>;
    <FLAGSTAT>;
    <FLAGSTAT>;

    $line = <FLAGSTAT>;
    $line =~ /(\d+) \+/;
    $paired_mapped_reads = $1;

    close(FLAGSTAT);
}

sub rawReadsFromFastq {
  my $file=shift;
  my $reads=`zcat $file | wc -l`;
  chomp $reads;
  $reads=$reads/2;
  return $reads;
}

sub rawReadsFromBam {
  my $file=shift;
  my $reads=`samtools view -c $file`;
  chomp $reads;
  return $reads;
}

sub rawReadsFromFastqc {
  my $file=shift;
  my $line;
  open(FASTQC, "<$file") || die "Could not open file $file\n$!\n";
###FastQC	0.10.0
#>>Basic Statistics	pass
##Measure	Value	
#Filename	S12386_1.fastq.gz	
#File type	Conventional base calls	
#Encoding	Illumina 1.5	
#Total Sequences	108352928	
#Filtered Sequences	0	
#Sequence length	100	
  <FASTQC>;
  <FASTQC>;
  <FASTQC>;
  <FASTQC>;
  <FASTQC>;
  <FASTQC>;
  $line = <FASTQC>;
  $line =~ /Total Sequences\t(\d+)/;
  my $reads=$1*2;
  return $reads;
}

sub printMetrics {
    open(OUT, ">$output") || die "Could not open file $output\n$!\n";
    print OUT "sequencer raw reads\t$raw_reads\n";
    print OUT "filtered paired reads\t$filtered_reads\n";
    print OUT "mapped reads\t$mapped_reads\n";
    print OUT "% mapped reads\t" . $mapped_reads/$filtered_reads . "\n";
    print OUT "paired mapped reads\t$paired_mapped_reads\n";
    print OUT "non-duplicated reads\t" . ($filtered_reads - $dup_reads) . "\n";
    print OUT "% duplication\t" . $dup_reads/$filtered_reads  . "\n";
#    print OUT "% mapped reads\t$percentMapped\n";
#    print OUT "non-duplicated mapped reads\t$nonDupMapped\n";
#    print OUT "% duplication\t$percentDup\n";
    close(OUT);
    print "\n\nMapping statistics have been writen in $output\n\n\n";
}

=head1 SYNOPSIS

 compileMappingStats.pl [arguments]

 File Arguments 
         --flagstat  <file>           specify the file containing the samtools flagstat results
         --output    <file>           specify the output file for this analysis
         --raw_bam   <file>           specify one or the other from raw_bam and raw_fastq (unfiltered reads)
         --raw_fastq <file>           specify one or the other from raw_bam and raw_fastq (unfiltered reads)

 Optional arguments : 
        -h, --help                    print help message
        -m, --man                     print complete documentation


 Function: Compile the mapping stats of an NGS run and print it in a text file.
           Used priorily to updateSamplePipelineStatus.pl which parses the mapping stat file in order to update the database.  

 Example: #download gene annotation database (for hg18 build) and save to humandb/ directory
      compileMappingStats.pl --raw ~/NGS/S12386/Illumina_HiSeq_Paired-IC-Exome_SS_50Mbp-2011_04_16/RAW --run ~/NGS/S12386/Illumina_HiSeq_Paired-IC-Exome_SS_50Mbp-2011_04_16/GATK_BWA.v37_50Mbp

 Version: $LastChangedDate: 2013-05-12 11:41 $

=head1 OPTIONS

=over 8

=item B<--help>

print a brief usage message and detailed explanation of options.

=item B<--man>

print the complete manual of the program.

=item B<--flagstat>

specify the file containing the samtools flagstat results

=item B<--output>

specify the output file for this analysis

=item B<--raw_bam>

specify one or the other from raw_bam and raw_fastq (unfiltered reads)

=item B<--raw_fastq>

specify one or the other from raw_bam and raw_fastq (unfiltered reads)

=back

-----------------------------------------

For questions or comments, please contact edouard.henrion@crchum.qc.ca.

=cut
