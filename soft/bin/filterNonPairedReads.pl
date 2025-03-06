#!/usr/bin/perl
use strict;
use Data::Dumper;

my $outName = $ARGV[0];
my $pair1 = $ARGV[1];
my $pair2 = $ARGV[2];
my $pair1FH;
my $pair2FH;
my $pair1outFH;
my $pair2outFH;
my $pairSingleFH;
#my $headChopChars = $ARGV[3];

open ($pair1FH, '-|', 'zcat', "$pair1") || die "Cannot Open pair 1 File";
open ($pair2FH, '-|', 'zcat', "$pair2") || die "Cannot Open pair 2 File";
my @pairedFiles = ($pair1FH,$pair2FH);

my %readCount;
my $totalRecords = 0;
for my $pairFH (@pairedFiles) {
  while (my $head = <$pairFH>) {
    <$pairFH>;
    <$pairFH>;
    <$pairFH>;
    $totalRecords++;
#    my $key = substr($head, 0, length($head)-3);
#    my $key = substr($head, 0, length($head)-$headChopChars-2);
    my $key = (split " ", $head)[0];
    my $value = $readCount{$key};
    if(defined $value) {
      $value++;
      $readCount{$key} = $value;
    }
    else {
      $readCount{$key} = 1;
    }
  }
  close($pairFH);
}

open ($pair1FH, '-|', 'zcat', "$pair1") || die "Cannot Open pair 1 File : $pair1 $!";
open ($pair2FH, '-|', 'zcat', "$pair2") || die "Cannot Open pair 2 File : $pair2 $!";
@pairedFiles = ($pair1FH,$pair2FH);
#open ($pair1outFH, "| gzip >$outName.pair1.fastq.gz") || die "Cannot Open pair out 1 File";
#open ($pair2outFH, "| gzip >$outName.pair2.fastq.gz") || die "Cannot Open pair out 2 File";
#open ($pairSingleFH, "| gzip >$outName.single.fastq.gz") || die "Cannot Open single File";
open ($pair1outFH, ">$outName.pair1.fastq") || die "Cannot Open pair out 1 File : $outName.pair1.fastq $!";
open ($pair2outFH, ">$outName.pair2.fastq") || die "Cannot Open pair out 2 File : $outName.pair2.fastq $!";
open ($pairSingleFH, ">$outName.single.fastq") || die "Cannot Open single File : $outName.single.fastq $!";

my @outputFiles = ($pair1outFH,$pair2outFH);

my $idx=0;
my $records=0;
my $pairs=0;
my $single=0;
print "Copying reads\n";
for my $pairFH (@pairedFiles) {
  my $outFH = $outputFiles[$idx];
  $idx++;

  while(my $head = <$pairFH>) {
    my $seq = <$pairFH>;
    my $head2 = <$pairFH>;
    my $qual = <$pairFH>;

#    my $key = substr($head, 0, length($head)-3);
#    my $key = substr($head, 0, length($head)-$headChopChars-2);
    my $key = (split " ", $head)[0];
    my $value = $readCount{$key};
    my $fh;
    $records++;
    if($value == 1) {
      $fh = $pairSingleFH;
      $single++;
    }
    elsif($value > 1) {
      $fh = $outFH;
      $pairs++;
    }
    else {
      die "No value found for: $key\n";
    }

    print $fh $head;
    print $fh $seq;
    print $fh $head2;
    print $fh $qual;
#    print "\rTotal: $totalRecords, Records: $records, Pairs: $pairs, Single: $single                      ";
  }
  close($outFH);
  close($pairFH);
}
close($pairSingleFH);
print "Total: $totalRecords, Records: $records, Pairs: $pairs, Single: $single\n";
