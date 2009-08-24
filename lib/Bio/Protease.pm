package Bio::Protease;
use Moose;
use MooseX::ClassAttribute;
use MooseX::Types::Moose qw(HashRef);
use namespace::autoclean;
use Bio::Protease::Types qw(ProteaseRegex ProteaseName);
extends 'Bio::ProteaseI';

class_has Specificities => (
    is      => 'ro',
    isa     => HashRef,
    lazy_build => 1,
);

sub BUILDARGS {
    my ($class, %args) = @_;

    $args{_regex} = $args{specificity};

    $class->SUPER::BUILDARGS(%args);
}

has _regex => (
    is  => 'ro',
    isa => ProteaseRegex,
    coerce => 1,
);

has specificity => (
    is  => 'ro',
    isa => ProteaseName,
    required => 1,
    coerce   => 1
);

augment _cuts => sub {
    my ($self, $peptide) = @_;

    if ( grep { $$peptide !~ /$_/ } @{$self->_regex} ) {
        return;
    }

    return 'yes, it cuts';

};

sub _build_Specificities {

    my %specificity_of = (
        'arg-c_proteinase'           => [ '.{3}R.{4}' ],
        'asp-n_endopeptidase'        => [ '.{4}D.{3}' ],
        'asp-n_endopeptidase_glu'    => [ '.{4}[DE].{3}' ],
        'bnps_skatole'               => [ '.{3}W.{4}' ],
        'caspase_1'                  => [ '[FWYL].[HAT]D[^PEDQKR].{3}' ],
        'caspase_2'                  => [ 'DVAD[^PEDQKR].{3}' ],
        'caspase_3'                  => [ 'DMQD[^PEDQKR].{3}' ],
        'caspase_4'                  => [ 'LEVD[^PEDQKR].{3}' ],
        'caspase_5'                  => [ '[LW]EHD.{4}' ],
        'caspase_6'                  => [ 'VE[HI]D[^PEDQKR].{3}' ],
        'caspase_7'                  => [ 'DEVD[^PEDQKR].{3}' ],
        'caspase_8'                  => [ '[IL]ETD[^PEDQKR].{3}' ],
        'caspase_9'                  => [ 'LEHD.{4}' ],
        'caspase_10'                 => [ 'IEAD.{4}' ],
        'chymotrypsin'               => [ '.{3}[FY][^P].{3}|.{3}W[^MP].{3}' ],
        'chymotrypsin_low'           => [ '.{3}[FLY][^P].{3}|.{3}W[^MP].{3}|.{3}M[^PY].{3}|.{3}H[^DMPW].{3}' ],
        'clostripain'                => [ '.{3}R.{4}' ],
        'cnbr'                       => [ '.{3}M.{4}' ],
        'enterokinase'               => [ '[DN][DN][DN]K.{4}' ],
        'factor_xa'                  => [ '[AFGILTVM][DE]GR.{4}' ],
        'formic_acid'                => [ '.{3}D.{4}' ],
        'glutamyl_endopeptidase'     => [ '.{3}E.{4}' ],
        'granzymeb'                  => [ 'IEPD.{4}' ],
        'hydroxylamine'              => [ '.{3}NG.{3}' ],
        'hcl'                        => [ '.{8}' ],
        'iodosobenzoic_acid'         => [ '.{3}W.{4}' ],
        'lysc'                       => [ '.{3}K.{4}' ],
        'lysn'                       => [ '.{4}K.{3}' ],
        'ntcb'                       => [ '.{4}C.{3}' ],
        'pepsin_ph1.3'               => [ '.[^HKR][^P][^R][FLWY][^P].{2}|.[^HKR][^P][FLWY].[^P].{2}' ],
        'pepsin'                     => [ '.[^HKR][^P][^R][FL][^P].{2}|.[^HKR][^P][FL].[^P].{2}' ],
        'proline_endopeptidase'      => [ '.{2}[HKR]P[^P].{3}' ],
        'proteinase_k'               => [ '.{3}[AFILTVWY].{4}' ],
        'staphylococcal_peptidase_i' => [ '.{2}[^E]E.{4}' ],
        'thermolysin'                => [ '.{3}[^XDE][AFILMV][^P].{2}' ],
        'thrombin'                   => [ '.{2}GRG.{3}|[AFGILTVM][AFGILTVWA]PR[^DE][^DE].{2}' ],
        'trypsin'                    => [ '.{2}(?!CKD).{6}', '.{2}(?!DKD).{6}', '.{2}(?!CKH).{6}', '.{2}(?!CKY).{6}', '.{2}(?!RRH).{6}', '.{2}(?!RRR).{6}', '.{2}(?!CRK).{6}',
                                        '.{3}[KR][^P].{3}|.{2}WKP.{3}|.{2}MRP.{3}' ]
    );

    return \%specificity_of;
}

__PACKAGE__->meta->make_immutable;

__END__

=pod

=head1 SYNOPSIS

    use Bio::Protease;
    my $protease = Bio::Protease->new(specificity => 'trypsin');

    my $protein = 'MRAERVIKP'; # Could also be a Bio::Seq object.

    # Perform a full digestion
    my @products = $protease->digest($protein);

    # products: ( 'MR', 'AER', 'VIKP' )

    # Get all the siscile bonds.
    my @sites = $protease->cleavage_sites($protein);

    # sites: ( 2, 5 )

    # Try to cut at a specific position.

    @products = $protease->cut($protein, 2);

    # products: ( 'MR', 'AERVIKP' )

=cut

=head1 DESCRIPTION

This module models the hydrolitic behaviour of a proteolytic enzyme.
Its main purpose is to predict the outcome of hydrolitic cleavage of a
peptidic substrate.

The enzyme specificity is currently modeled for 32 enzymes/reagents.
This models are somewhat simplistic as they are largely regex-based, and
do not take into account subtleties such as kinetic/temperature effects,
accessible solvent area, secondary or tertiary structure elements.
However, the module is flexible enough to allow the inclusion of any of
these effects via subclassing from the module's interface, Bio::ProteaseI.

=cut

=head1 Attributes And Methods

=head2 specificity

Set the enzyme's specificity. Required. Could be either of:

=over 4

=item * an enzyme name: e.g. 'enterokinase'

    my $enzyme = Bio::Protease->new(specificity => 'enterokinase');

There are currently definitions for 36 enzymes/reagents. See
C<Specificities>.

=item * an array reference of regular expressions:

    my $motif = ['MN[ED]K[^P].{3}'],

    my $enzyme = Bio::Protease->new(specificity => $motif);

The motif should always describe an 8-character long peptide. When a an
octapeptide matches the regex, its 4th peptidic bond (ie, between the
4th and 5th letter) will be marked for cleaving or reporting.

For example, the peptide AMQRNLAW is recognized as follows:

    .----..----.----..----. .-----.-----.-----.-----.
    | A  || M  | Q  || R  |*|  N  |  L  |  A  |  W  |
    |----||----|----||----|^|-----|-----|-----|-----|
    | P4 || P3 | P2 || P1 ||| P1' | P2' | P3' | P4' |
    '----''----'----''----'|'-----'-----'-----'-----'
                      cleavage site

Some specificity rules can only be described with more than one regular
expression (See the case for trypsin, for example). To account for those
cases, the array reference could contain an arbitrary number of regexes,
all of which should match the given octapeptide.

In the case your particular specificity rule requires an "or" clause,
you can use the "|" separator in a single regex.

=back

=cut

=head2 digest($substrate)

Performs a complete digestion of the peptide argument, returning a list
with possible products. It does not do partial digests (see method
C<cut> for that).

    my @products = $enzyme->digest($protein);

=head2 cut($substrate, $i)

Attempt to cleave $substrate at the C-terminal end of the $i-th residue
(ie, at the right). If the bond is indeed cleavable (determined by the
enzyme's specificity), then a list with the two products of the
hydrolysis will be returned. Otherwise, returns false.

    my @products = $enzyme->cut($peptide, $position);

=head2 cleavage_sites($protein)

Returns a list with siscile bonds (bonds susceptible to be cleaved as
determined by the enzyme's specificity). Bonds are numbered starting
from 1, from N to C-terminal.

=cut

=head1 Class Attributes

=head2 Specificities

A hash reference with all the available regexep-based specificities. The
keys are the specificity names, the value is an arrayref with the
regular expressions that define it.

    my @protease_pool = do {
        Bio::Protease->new(specificity => $_)
            for keys %{Bio::Protease->Specificities};
    }

As a rule, all specificity names are lower case. Currently, they include:

=over 2

=item * arg-cproteinase

=item * asp-n endopeptidase

=item * asp-n endopeptidase glu

=item * bnps skatole

=item * caspase 1

=item * caspase 2

=item * caspase 3

=item * caspase 4

=item * caspase 5

=item * caspase 6

=item * caspase 7

=item * caspase 8

=item * caspase 9

=item * caspase 10

=item * chymotrypsin

=item * chymotrypsin low

=item * clostripain

=item * cnbr

=item * enterokinase

=item * factor xa

=item * formic acid

=item * glutamyl endopeptidase

=item * granzymeb

=item * hydroxylamine

=item * iodosobenzoic acid

=item * lysc

=item * lysn

=item * ntcb

=item * pepsin ph1.3

=item * pepsin

=item * proline endopeptidase

=item * proteinase k

=item * staphylococcal peptidase i

=item * thermolysin

=item * thrombin

=item * trypsin

=back

For a complete description of their specificities, you can check out
<link here>, or look at the regular expressions of their definitions
in this same file.

=cut

=head1 SEE ALSO
L<Bio::Tools::SeqPattern>,  L<Bio::Seq>.

PeptideCutter This module's idea is largely based on Expasy's
PeptideCutter(L<http://www.expasy.ch/tools/peptidecutter/>). For more
information on the experimental evidence that supports both the
algorithm and the specificity definitions, check their page.

=head1 AUTHOR

Bruno Vecchi, C<< <vecchi.b at gmail.com> >>

=head1 BUGS

Please report any bugs or feature request via the github issue tracker
(L<http://github.com/brunoV/bio-protease/issues>).

=cut

#Please report any bugs or feature requests to C<bug-bio-protease at rt.cpan.org>, or through
#the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Protease>.  I will be notified, and then you'll
#automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::Protease


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-Protease>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-Protease>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-Protease>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-Protease/>

=back


=head1 COPYRIGHT & LICENSE

Copyright 2009 Bruno Vecchi, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.


=cut
