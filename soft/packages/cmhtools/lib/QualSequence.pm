# $Id: QualSequence.pm,v 1.4 2006/05/24 23:32:48 cmhall Exp $

######################################################################

=head1 QualSequence.pm

Encapsulates a sequence of dna quality values.

=cut

######################################################################
package QualSequence;

use lib '.';
use base Sequence;
use Util;

######################################################################

=over 4

=item new($title, $sequence)

Constructs a new QualSequence object and returns a reference to it.
$title is the title of the new sequence and $sequence is the actual
sequence data.

=back

=cut

######################################################################
sub new {
    my ($proto, $title, $sequence) = @_;

    $title    ||= '';
    $sequence ||= '';

    my $class = ref($proto) || $proto;
    my $self = $class->SUPER::new($title, $sequence);
    return $self;
}

######################################################################

=over 4

=item to_string($width)

Returns a string representation of this QualSequence with sequence
data wrapped to $width.

=back

=cut

######################################################################
sub to_string {
    my ($self, $width) = @_;

    $width ||= '';

    my $qperline = 0;
    if($width) {
        $qperline = &Util::max(int($width / 3), 1);
    }
    my @seqs = split(/\s+/, $self->{'sequence'});
    my $count = 0;
    my $space = '';
    my $seq = '';
    foreach my $s (@seqs) {
        if(length($s) == 2) {
            $seq .= $space.$s;
        }
        elsif(length($s) == 1) {
            $seq .= $space." ".$s;
        }
        else {
            die "found a really strange qual value: $s";
        }

        if($s == 0) {
            print STDERR "found a '0' quality value in: ".
                         $self->{'title'}."\n";
        }

        $space = ' ';
        $count++;
        if($count == $qperline) {
            $count = 0;
            $space = "\n";
        }
    }

    return ">".$self->{'title'}."\n$seq\n";
}

######################################################################

=over 4

=item length()

Returns the length of this QualSequence.

=back

=cut

######################################################################
sub length {
    my ($self) = @_;

    my @s = split(/\s+/, $self->{'sequence'});

    return scalar(@s);
}

######################################################################

=over 4

=item trim($left, $right)

Trims the ends of this QualSequence.  $left and $right are 1-based
exclusive trimming indexes.  For example, if the sequence is: 21 22 23
24 25 26, and trim(2, 4) is called, the new sequence will be: 22 23
24.  Another way to think of it is that the 1-based indexes will be
'included' in the remaining sequence after trimming is done.

=back

=cut

######################################################################
sub trim {
    my ($self, $left, $right) = @_;

    $left ||= 0;
    $right ||= 9999999999;

    $left--;  # trimming coordinates are 1-based
    $right--;

    my @bases = split(/\s+/, $self->{'sequence'});
    my @trimmed_bases = splice(@bases, $left, $right - $left + 1);
    my $newseq = '';
    my $space = '';
    foreach my $s (@trimmed_bases) {
        $newseq .= "$space$s";
        $space = ' ';
    }

    $self->{'sequence'} = $newseq;
}

######################################################################

=over 4

=item print_avg_qual()

Returns a string representation of the average quality value in this
QualSequence.

=back

=cut

######################################################################
sub print_avg_qual {
    my ($self) = @_;

    my @seqs = split(/\s+/, $self->{'sequence'});
    my $count = 0;
    my $total = 0;
    my $avg = 0;
    my $s = "";

    foreach my $s (@seqs) {
        $total += $s;
        $count++;
    }

    if($count != 0) {
        $avg = $total / $count;
    }

    $s = sprintf("%-30s%3.3f\n", $self->{'title'}.": ", $avg);

    return $s;
}

1;
