# $Id: Util.pm,v 1.3 2006/05/24 23:32:48 cmhall Exp $

######################################################################

=head1 Util.pm

Utility functions that may be useful in many modules, but are not
found in perl itself.

=cut

######################################################################
package Util;

######################################################################

=over 4

=item wrap($string, $width)

Wraps the text pointed to by the scalar reference $string to $width by
inserting newlines.

=back

=cut

######################################################################
sub wrap {
    my ($string, $width) = @_;

    $string ||= '';
    $width  ||= '';

    if(!$string) {
        return;
    }

    if($width >= length($$string)) {
        return;
    }

    $$string = join("\n", split(/(.{$width})/, $$string));
    $$string =~ s/^\n//;
    $$string =~ s/\n\n/\n/g;
}

######################################################################

=over 4

=item trim($string)

Trims white space from the text pointed to by the scalar reference
$string.

=back

=cut

######################################################################
sub trim {
    my ($string) = @_;

    $string ||= '';

    if(!$string) {
        return;
    }

    $$string =~ s/^\s+|\s+$//g;
}

######################################################################

=over 4

=item shrink_white_space($string)

Trims the text pointed to by the scalar reference $string with the
trim() method of this module, then replaces all sequences of one or
more space, tab, or newline with a single space.

=back

=cut

######################################################################
sub shrink_white_space {
    my ($string) = @_;

    $string ||= '';

    if(!$string) {
        return;
    }

    &trim($string);
    $$string =~ s/\s+/ /g;
}

######################################################################

=over 4

=item max($a, $b)

Returns the maximum value of $a and $b.  If $a is equal to $b, returns
the value of $a.

=back

=cut

######################################################################
sub max {
    my ($a, $b) = @_;

    $a ||= 0;
    $b ||= 0;

    return $a >= $b ? $a : $b;
}

######################################################################

=over 4

=item min($a, $b)

Returns the minimum value of $a and $b.  If $a is equal to $b, returns
the value of $a.

=back

=cut

######################################################################
sub min {
    my ($a, $b) = @_;

    $a ||= 0;
    $b ||= 0;

    return $a <= $b ? $a : $b;
}

1;
