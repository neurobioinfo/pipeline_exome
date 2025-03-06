# $Id: Fasta.pm,v 1.3 2006/05/24 23:32:48 cmhall Exp $

######################################################################

=head1 Fasta.pm

Encapsulates a generic fasta file.  This module is intended to be used
as an abstract class.  It should be subclassed, and any functionality
pertaining to the type of sequence data should be implemented in the
subclass.

=cut

######################################################################
package Fasta;

use lib '.';
use Util;
use Sequence;

######################################################################

=over 4

=item new()

Constructs an new generic fasta object and returns a reference to it.

=back

=cut

######################################################################
sub new {
    my ($proto) = @_;

    my $class = ref($proto) || $proto;
    my $self = {};

    $self->{'seqs'} = [];

    bless($self, $class);
    return $self;
}

######################################################################

=over 4

=item get_seqs()

Returns an array containing all the Sequence objects contained in this
Fasta object.

=back

=cut

######################################################################
sub get_seqs {
    my ($self) = @_;

    return @{$self->{'seqs'}};
}

1;
