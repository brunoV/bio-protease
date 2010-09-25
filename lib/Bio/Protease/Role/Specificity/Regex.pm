package Bio::Protease::Role::Specificity::Regex;

# ABSTRACT: A role that implements a regex-based specificity

use Moose::Role;
use Bio::Protease::Types 'ProteaseRegex';

=attr regex

A C<ProteaseRegex>, which is basically an array reference of regular
expressions that describe the protease specificity. It can coerce from a
single regular expression into a single-element array of regexps.  Any
of the regexes in the array should match a given substrate for it to be
cleavable.

=cut

has regex => (
    is  => 'ro',
    isa => ProteaseRegex,
    coerce => 1,
);

sub _cuts {
    my ($self, $peptide) = @_;

    if ( grep { $peptide !~ /$_/ } @{$self->regex} ) {
        return;
    }

    return 'yes, it cuts';
}

=head1 SYNOPSIS

    package My::Protease;
    use Moose;

    with qw(Bio::ProteaseI Bio::Protease::Role::Specificity::Regex);

    package main;

    my $p = My::Protease->new( regex => qr/.{3}AC.{3}/ ); # coerces to [ qr/.../ ];

    my @products = $p->digest( 'AAAACCCC' );

    # @products: ('AAAA', 'CCCC')

=head1 DESCRIPTION

This role implements a regexp-based specificity for a class that also
consumes the L<Bio::ProteaseI> role. A peptide will be cleaved if any of
the regexes provided at construction time matches it. The regexes should
be tailored for 8-residue-long peptides, the cleavage site being between
the fourth and fifth residues.

For instance, if the specificity could be described as "cuts after
lysine or arginine", the appropriate regular expression would be
C<qr/.{3}[KR].{4}/>.

=cut

1;
