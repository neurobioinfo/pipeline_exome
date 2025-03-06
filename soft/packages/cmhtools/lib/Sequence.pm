# $Id: Sequence.pm,v 1.3 2006/05/24 23:32:48 cmhall Exp $

######################################################################

=head1 Sequence.pm

Encapsulates a generic sequence.

=cut

######################################################################
package Sequence;

use lib '.';

######################################################################

=over 4

=item new($title, $sequence)

Constructs a new generic sequence object and returns a reference to
it.  $title is the title of the sequence and $sequence is the actual
sequence data of the sequence.

=back

=cut
######################################################################
sub new {
    my ($proto, $title, $sequence) = @_;

    $title     ||= '';
    $sequence  ||= '';
    $formatter ||= '';

    my $class = ref($proto) || $proto;
    my $self = {};

    $self->{'title'}     = $title;
    $self->{'sequence'}  = $sequence;

    bless($self, $class);
    return $self;
}

######################################################################

=over 4

=item get_title()

Returns the entire title of this sequence.

=back

=cut

######################################################################
sub get_title {
    my ($self) = @_;

    return $self->{'title'};
}

######################################################################

=over4

=item get_title_leftmost()

Returns the leftmost space delimited element in the title line of this
sequence.

=back

=cut

######################################################################
sub get_title_leftmost {
    my ($self) = @_;

    my @title_elements = split(/\s+/, $self->{'title'});

    return $title_elements[0];
}

######################################################################

=over 4

=item set_title($title)

Sets the title of this sequence to $title.

=back

=cut

######################################################################
sub set_title {
    my ($self, $title) = @_;

    $title ||= '';

    $self->{'title'} = $title;
}

1;
