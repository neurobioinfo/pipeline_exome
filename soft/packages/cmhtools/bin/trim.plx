#! /usr/local/bin/perl -w
#
# $Id: trim.plx,v 1.5 2006/05/24 23:32:39 cmhall Exp $
#
# Actually cuts off the portions of reads that are trimmed according to lucy.
#
# example:
#
# cmhall@bio3:~/dev/cmhtools2/bin> cat test.fas
# >test.sequence 0 0 0 3 25
# AGCTATCAAAATTGGCACACACAGGTAGTTG
# cmhall@bio3:~/dev/cmhtools2/bin> cat test.qual
# >test.sequence
# 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 
# 20 20 20 20 20 20 20 20 20 20 20 20 20 
# cmhall@bio3:~/dev/cmhtools2/bin> ./trim.plx --seq=test.fas --qual=test.qual
# cmhall@bio3:~/dev/cmhtools2/bin> cat test.fas.t
# >test.sequence
# CTATCAAAATTGGCACACACAGG
# cmhall@bio3:~/dev/cmhtools2/bin> cat test.qual
# >test.sequence
# 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 20 
# 20 20 20 20 20 20 20 20 20 20 20 20 20 
# cmhall@bio3:~/dev/cmhtools2/bin>

use strict;
use Env qw(CMHTOOLS);
use lib "$CMHTOOLS/lib";
use DnaFasta;
use QualFasta;
use Getopt::Mixed;

#defaults
my $DEFAULT_SEQ_FILE  = '';
my $DEFAULT_QUAL_FILE = '';
my $DEFAULT_SEQ_WIDTH = 70;
my $DEFAULT_QUAL_WIDTH = 60;

$| = 1;

my $usage = <<END;
usage: $0
    [--help]
    --seq=<seq-file>
    --qual=<qual-file>

--help
    Print this help message.

--seq=<seq-file>
    The sequence file.

--qual=<qual-file>
    The quality file.

END

Getopt::Mixed::init(
    "help",
    "seq=s",
    "qual=s",
);
my %opts = ();
while((my $option, my $value) = Getopt::Mixed::nextOption()) {
    $opts{$option} = defined($value) ? $value : 1;
}
Getopt::Mixed::cleanup();

if(defined($opts{'help'})) {
    print "$usage\n";
    exit;
}

my $seq_file  = $opts{'seq'}  || $DEFAULT_SEQ_FILE;
my $qual_file = $opts{'qual'} || $DEFAULT_QUAL_FILE;

# check that both --seq and --qual were given.
if(!$seq_file || !$qual_file) {
    die "--seq and --qual are required parameters";
}

my $sf = new DnaFasta($seq_file);
my $qf = new QualFasta($qual_file);

my @seqs = $sf->get_seqs();
my @quals = $qf->get_seqs();

my %leftt  = ();
my %rightt = ();

open SEQ, ">$seq_file.t" or die;
foreach my $s (@seqs) {
    $s->get_title() =~ /^(.*) (.*) (.*) (.*) (.*) (.*)$/; 
    $s->set_title($1);
    $leftt{$1} = $5;
    $rightt{$1} = $6;
    $s->trim($5, $6);
    print SEQ $s->to_string($DEFAULT_SEQ_WIDTH);
}
close SEQ;

open QUAL, ">$qual_file.t" or die;
foreach my $s (@quals) {
    my $t = $s->get_title();
    $s->set_title($t);
    $s->trim($leftt{$t}, $rightt{$t});
    print QUAL $s->to_string($DEFAULT_QUAL_WIDTH);
}
close QUAL;
