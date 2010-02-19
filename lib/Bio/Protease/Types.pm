package Bio::Protease::Types;

# ABSTRACT: Specific types for Bio::Protease

use MooseX::Types::Moose qw(Str ArrayRef);
use MooseX::Types -declare => [qw(ProteaseName ProteaseRegex)];
use namespace::autoclean;
use Carp qw(croak);

subtype ProteaseName, as Str;

subtype ProteaseRegex, as ArrayRef;

coerce ProteaseRegex, from Str, via { _str_to_prot_regex($_) };

coerce ProteaseName,
    from ArrayRef, via { 'custom' };

sub _str_to_prot_regex {
    my $specificity = shift;
    my $specificity_of = Bio::Protease->Specificities;

    croak "Not a known specificity\n"
        unless $specificity ~~ %$specificity_of;

    return $specificity_of->{$specificity};
}

__PACKAGE__->meta->make_immutable;

=head1 DESCRIPTION

This module defines specific types and type coercions to be used by
L<Bio::Protease>. It should not be used by end users or consumer of the
Bio::ProteaseI role.

=cut
