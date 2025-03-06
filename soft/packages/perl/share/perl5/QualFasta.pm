# $Id: QualFasta.pm,v 1.4 2006/05/24 23:32:48 cmhall Exp $

######################################################################

=head1 QualFasta.pm

Encapsulates a multifasta file of quality values.

=cut

######################################################################
package QualFasta;

use lib '.';
use base Fasta;
use Util;
use QualSequence;

######################################################################

=over 4

=item new($path)

Constructs a new QualFasta object and returns a reference to it.  $path
is the path, full or relative, to a fasta file containing one or more
sequences of quality values.

=back

=cut

######################################################################
sub new {
    my ($proto, $path) = @_;

    $path      ||= '';

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
                Util::shrink_white_space(\$seq);
                push(@{$self->{'seqs'}}, new QualSequence($title, $seq));
            }
            $title = $1;
            $seq = '';
        }
        else {
            $seq .= " $_";
        }
    }
    if($path) {
        close $IN;
    }

    if($title && $seq) {
        Util::shrink_white_space(\$seq);
        push(@{$self->{'seqs'}}, new QualSequence($title, $seq));
    }

    return $self;
}

######################################################################

=over 4

=item print($width, $sorted)

prints the fasta file represented by this object.  Each line of
quality value data is wrapped to $width.  If $sorted evaluates to
true, the sequences are sorted alphabetically by their title.

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

######################################################################

=over 4

=item print_avg_qual($sorted)

For each sequence in this file, prints the average quality value.  If
$sorted evaluates to true, the results are sorted alphabetically by
the title of the sequence.

=back

=cut

######################################################################
sub print_avg_qual {
    my ($self, $sorted) = @_;

    $sorted ||= '';

    foreach my $s (sort {$a->get_title() cmp $b->get_title()} @{$self->{'seqs'}}) {
        print $s->print_avg_qual();
    }
}


1;
