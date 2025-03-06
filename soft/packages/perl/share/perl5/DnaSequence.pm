# $Id: DnaSequence.pm,v 1.3 2006/05/24 23:32:48 cmhall Exp $

######################################################################

=head1 DnaSequence.pm

Encapsulates a sequence of dna.

=cut

######################################################################
package DnaSequence;

use lib '.';
use base Sequence;
use Util;

######################################################################

=over 4

=item new($title, $sequence)

Constructs a new DnaSequence object and returns a reference to it.
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

Returns a string representation of this DnaSequence with sequence data
wrapped to $width.

=back

=cut

######################################################################
sub to_string {
    my ($self, $width) = @_;

    $width ||= '';

    my $seq = $self->{'sequence'};
    if($width) {
        $width = &Util::max($width, 1);
        &Util::wrap(\$seq, $width);
    }

    return ">".$self->{'title'}."\n$seq\n";
}

######################################################################

=over 4

=item length()

Returns the length of this DnaSequence.

=back

=cut

######################################################################
sub length {
    my ($self) = @_;

    return length($self->{'sequence'});
}

######################################################################

=over 4

=item trim($left, $right)

Trims the ends of this DnaSequence.  $left and $right are 1-based
exclusive trimming indexes.  For example, if the sequence is:
AAACCCGGGTTT, and trim(4, 10) is called, the new sequence will be:
CCCGGGT.  Another way to think of it is that the 1-based indexes will
be 'included' in the remaining sequence after trimming is done.

=back

=cut

######################################################################
sub trim {
    my ($self, $left, $right) = @_;

    $left  ||= 1;
    $right ||= 9999999999;

    $left--;  # trimming coordinates are 1-based
    $right--;

    $self->{'sequence'} = substr($self->{'sequence'}, $left,
                                 $right - $left + 1);
}

1;
