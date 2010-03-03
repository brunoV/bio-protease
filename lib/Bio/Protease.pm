package Bio::Protease;
use Moose;
use MooseX::ClassAttribute;
use MooseX::Types::Moose qw(HashRef);
use namespace::autoclean;
use Bio::Protease::Types qw(ProteaseRegex ProteaseName);
with 'Bio::ProteaseI';

# ABSTRACT: Digest your protein substrates with customizable specificity

has _regex => (
    is  => 'ro',
    isa => ProteaseRegex,
    init_arg => 'specificity',
    coerce => 1,
);

=attr specificity

Set the enzyme's specificity. Required. Could be either of:

=over 4

=item * an enzyme name: e.g. 'enterokinase'

    my $enzyme = Bio::Protease->new(specificity => 'enterokinase');

There are currently definitions for 36 enzymes/reagents. See
L<Specificities>.

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

has specificity => (
    is  => 'ro',
    isa => ProteaseName,
    required => 1,
    coerce   => 1
);

sub _cuts {
    my ($self, $peptide) = @_;

    if ( grep { $peptide !~ /$_/ } @{$self->_regex} ) {
        return;
    }

    return 'yes, it cuts';

};

=attr Specificities

This B<class attribute> contains a hash reference with all the available
regexep-based specificities. The keys are the specificity names, the
value is an arrayref with the regular expressions that define them.

    my @protease_pool = do {
        Bio::Protease->new(specificity => $_)
            for keys %{Bio::Protease->Specificities};
    }

As a rule, all specificity names are lower case. Currently, they include:

=over 2

=item * arg-cproteinase

=item * asp-n_endopeptidase

=item * asp-n_endopeptidase_glu

=item * bnps_skatole

=item * caspase_1

=item * caspase_2

=item * caspase_3

=item * caspase_4

=item * caspase_5

=item * caspase_6

=item * caspase_7

=item * caspase_8

=item * caspase_9

=item * caspase_10

=item * chymotrypsin

=item * chymotrypsin_low

=item * clostripain

=item * cnbr

=item * enterokinase

=item * factor_xa

=item * formic_acid

=item * glutamyl_endopeptidase

=item * granzymeb

=item * hydroxylamine

=item * iodosobenzoic_acid

=item * lysc

=item * lysn

=item * ntcb

=item * pepsin_ph1.3

=item * pepsin

=item * proline_endopeptidase

=item * proteinase_k

=item * staphylococcal_peptidase i

=item * thermolysin

=item * thrombin

=item * trypsin

=back

For a complete description of their specificities, you can check out
L<http://www.expasy.ch/tools/peptidecutter/peptidecutter_enzymes.html>,
or look at the regular expressions of their definitions in this same
file.

=cut

class_has Specificities => (
    is      => 'ro',
    isa     => HashRef,
    lazy_build => 1,
);

sub _build_Specificities {

    my %specificity_of = (
        'alcalase'                   => [ '.{3}[MYFLIVW].{4}'],
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

=head1 SYNOPSIS

    use Bio::Protease;
    my $protease = Bio::Protease->new(specificity => 'trypsin');

    my $protein = 'MRAERVIKP';

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

The enzyme specificity is currently modeled for 36 enzymes/reagents.
This models are somewhat simplistic as they are largely regex-based, and
do not take into account subtleties such as kinetic/temperature effects,
accessible solvent area, secondary or tertiary structure elements.
However, the module is flexible enough to allow the inclusion of any of
these effects by consuming the module's interface, L<Bio::ProteaseI>.
Alternatively, if your desired specificity can be correctly described by
a regular expression, you can pass it as a string to the specificity
attribute at construction time. See L<specificity> below.

=cut

=method digest

Performs a complete digestion of the peptide argument, returning a list
with possible products. It does not do partial digests (see method
C<cut> for that).

    my @products = $enzyme->digest($protein);

=method cut

Attempt to cleave C<$peptide> at the C-terminal end of the C<$i>-th
residue (ie, at the right). If the bond is indeed cleavable (determined
by the enzyme's specificity), then a list with the two products of the
hydrolysis will be returned. Otherwise, returns false.

    my @products = $enzyme->cut($peptide, $i);

=method cleavage_sites

Returns a list with siscile bonds (bonds susceptible to be cleaved as
determined by the enzyme's specificity). Bonds are numbered starting
from 1, from N to C-terminal. Takes a string with the protein sequence
as an argument:

    my @sites = $enzyme->cleavage_sites($peptide);

=method is_substrate

Returns true or false whether the peptide argument is a substrate or
not. Esentially, it's equivalent to calling C<cleavage_sites> in boolean
context, but with the difference that this method short-circuits when it
finds its first cleavable site. Thus, it's useful for CPU-intensive
tasks where the only information required is whether a polypeptide is a
substrate of a particular enzyme or not 

=cut


=head1 SEE ALSO

=over

=item * PeptideCutter

This module's idea is largely based on Expasy's
PeptideCutter (L<http://www.expasy.ch/tools/peptidecutter/>). For more
information on the experimental evidence that supports both the
algorithm and the specificity definitions, check their page.

=back

=cut
