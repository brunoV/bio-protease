package Bio::Protease;
use 5.010_000;
use Moose;
use MooseX::ClassAttribute;
use Moose::Util::TypeConstraints;
use Carp;
use Memoize;
memoize ('cleavage_sites');
memoize ('is_substrate');
memoize ('_cuts');

my %specificity_of;

=head1 NAME

Bio::Protease - A class that describes a proteolytic enzyme.

=head1 VERSION

Version 0.01

=cut

class_has 'Specificities' => (
    is      => 'ro',
    isa     => 'ArrayRef',
    default => sub { [ keys %specificity_of ] },
);

subtype 'Bio::Protease::Specificity' => as 'CodeRef';
subtype 'Bio::Protease::SeqPattern'
    => as class_type('Bio::Tools::SeqPattern');

has _specif_code => (
   is         => 'ro',
   isa        => 'Bio::Protease::Specificity',
   init_arg   => 'specificity',
   required   => 1,
   coerce     => 1,
);

coerce 'Bio::Protease::Specificity',
    from 'Str',                       via { _str_to_specificity($_)     },
    from 'Bio::Protease::SeqPattern', via { _pattern_to_specificity($_) };

sub _cuts {
    my $self = shift;
    return $self->_specif_code->(@_);
}

sub _str_to_specificity {
    my $specificity = shift;

    croak "Not a known specificity\n"
        unless $specificity ~~ @{__PACKAGE__->Specificities};

    my @regexes = @{$specificity_of{$specificity}};

    my $coderef = _regex_to_coderef(@regexes);
    return $coderef;
}

sub _pattern_to_specificity {
    my $pattern_obj = shift;
    my $regex = $pattern_obj->str;

    my $coderef = _regex_to_coderef($regex);
    return $coderef;
}

sub _regex_to_coderef {
    my @regexes = @_;

    return sub {
        my $substrate = shift;
        my $length = length $substrate;

        if ( $length < 8 ) {
            if ( $length > 4 ) {
                $substrate .= 'X' x (8 - length $substrate);
            }
            else { return }
        }

        if ( grep { $substrate !~ /$_/ } @regexes ) {
            return;
        } else {
            return 1;
        }
    }
}

sub cut {
    my ( $self, $substrate, $pos ) = @_;

    unless ( defined $pos and $pos > 0 and $pos <= length $substrate ) {

        carp "Incorrect position.";
        return;
    }

    $substrate = uc $substrate;
    $substrate = 'XXX'. $substrate;
    $pos += 3;

    my $pep = substr($substrate, $pos - 4, 8);

    if ( $self->_cuts($pep) ) {
        my $product = substr($substrate, 0, $pos);
        substr($substrate, 0, $pos) = '';

        s/X//g for ($product, $substrate);

        return ($product, $substrate);
    }
    else { return }
}

sub digest {
    my ( $self, $substrate ) = @_;
    $substrate = uc $substrate;
    my @products;
    my ($i, $j) = (0, 0);

    $substrate = 'XXX' . $substrate;
    while ( my $pep = substr($substrate, $i, 8) ) {
        if ( $self->_cuts( $pep ) ) {
            my $product = substr($substrate, $j, $i + 4 - $j);
            push @products, $product;

            $j = $i + 4;
        }
        $i++;
    }
    push @products, substr($substrate, $j - length($substrate));

    s/X//g for @products[0, -1];

    return @products;
}

sub cleavage_sites {
    my ( $self, $substrate ) = @_;
    $substrate = uc $substrate;
    my @sites;
    my $i = 1;

    $substrate = 'XXX' . $substrate;
    while ( my $pep = substr($substrate, $i-1, 8 ) ) {
        if ( $self->_cuts( $pep ) ) { push @sites, $i };
        ++$i;
    }
    return @sites;
}

sub is_substrate {
    my ($self, $substrate) = @_;

    for my $pos (1 .. length $substrate) {
        return 1 if $self->cut($substrate, $pos);
    }

    return;
}



our $VERSION = '0.01';

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
these effects via user-defined specificity functions.

=cut

=head1 Attributes And Methods

=head2 specificity

Set the enzyme's specificity. Required. Could be either of:

=over 4

=item * an enzyme name: e.g. 'enterokinase'

    my $enzyme = Bio::Protease->new(specificity => 'enterokinase');

There are currently definitions for 36 enzymes/reagents. See C<Specificities>.

=item * a Bio::Tools::SeqPattern object.

    my $motif = Bio::Tools::SeqPattern->new(
        -SEQ  => 'MN[ED]K[^P].{3}',
        -TYPE => 'Amino',
    );

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


=item * a code reference.

    my $specificity = sub {
        my $peptide = shift;

        # ... some code that decides
        # ... if $peptide should be cut or not

        if ( peptide_should_be_cut ) { return 1 }
        else                         { return   }
    }

    my $enzyme = Bio::Protease->new(specificity => $coderef);

The code reference will be used by the methods C<digest>, C<cut> and
C<cleavage_sites>. It will always be passed a string with a length of 8
characters; if the coderef returns true, then the peptide bond between
the 4th and 5th residues will be marked as siscile, and the appropiate
action will be performed depending on which method was called.

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

A list with all the available regexep-based specificities.

    my @protease_pool = do {
        Bio::Protease->new(specificity => $_)
            for Bio::Protease->Specificities;
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

### Enzyme specificities

BEGIN {
    %specificity_of = (
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
}

no Moose;
__PACKAGE__->meta->make_immutable;
1;
