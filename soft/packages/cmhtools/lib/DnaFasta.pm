# $Id: DnaFasta.pm,v 1.3 2006/05/24 23:32:48 cmhall Exp $

######################################################################

=head1 DnaFasta.pm

Encapsulates a multifasta file of sequences.

=cut

######################################################################
package DnaFasta;

use lib '.';
use base Fasta;
use Util;
use DnaSequence;

######################################################################

=over 4

=item new($path)

Constructs a new DnaFasta object and returns a reference to it.  $path
is the path, full or relative, to a fasta file containing one or more
dna sequences.

=back

=cut

######################################################################
sub new {
    my ($proto, $path) = @_;

    $path ||= '';

    my $class = ref($proto) || $proto;
    my $self = $class->SUPER::new();

    my $seq = '';

    my $IN = 'IN';
    if($path) {
        open $IN, $path or die "failed to open $path";
    }
    else {
        $IN = 'STDIN';
    }
    while(<$IN>) {
        Util::trim(\$_);
        /^$/ && next;
        if(/^>(.*)$/) {
            if($title) {
                push(@{$self->{'seqs'}}, new DnaSequence($title, $seq));
            }
            $title = $1;
            $seq = '';
        }
        else {
            $seq .= $_;
        }
    }
    if($path) {
        close $IN;
    }

    if($title && $seq) {
        push(@{$self->{'seqs'}}, new DnaSequence($title, $seq));
    }

    return $self;
}

######################################################################

=over 4

=item print($width, $sorted)

Prints the fasta file represented by this object.  Each line of
sequence data is wrapped to $width.  If $sorted evaluates to true, the
sequences are sorted alphabetically by their title.

=back

=cut

######################################################################
sub print {
    my ($self, $width, $sorted) = @_;

    $width  ||= '';
    $sorted ||= '';

    if($width) {
        if($sorted) {
            foreach my $s (sort {$a->get_title() cmp $b->get_title()} @{$self->{'seqs'}}) {
                print $s->to_string($width);
            }
        }
        else {
            foreach my $s (@{$self->{'seqs'}}) {
                print $s->to_string($width);
            }
        }
    }
    else {
        if($sorted) {
            foreach my $s (sort {$a->get_title() cmp $b->get_title()} @{$self->{'seqs'}}) {
                print $s->to_string();
            }
        }
        else {
            foreach my $s (@{$self->{'seqs'}}) {
                print $s->to_string();
            }
        }
    }
}

1;
