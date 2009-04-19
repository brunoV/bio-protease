package Bio::Protease;
use Modern::Perl;
use Moose;
use MooseX::ClassAttribute;
use Moose::Util::TypeConstraints;
use Carp;

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

subtype 'Specificity' => as 'CodeRef';
subtype 'SeqPattern'
    => as class_type('Bio::Tools::SeqPattern');

has _specif_code => (
   is         => 'ro',
   isa        => 'Specificity',
   init_arg   => 'specificity',
   required   => 1,
   coerce     => 1,
);

coerce 'Specificity',
    from 'Str', via { &_str_to_specificity($_) },
    from 'SeqPattern', via { _pattern_to_specificity($_) };

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
}

sub _regex_to_coderef {
    my @regexes = @_;

    return sub {
        my $substrate = shift;
        unless ( length $substrate == 8 ) {
            $substrate .= 'X' x (8 - length $substrate);
        }

        if ( grep { $substrate !~ /$_/ } @regexes ) {
            return;
        } else {
            return 1;
        }
    }
}
sub _pattern_to_specificity {
    my $pattern_obj = shift;
    my $regex = $pattern_obj->str;

    my $coderef = _regex_to_coderef($regex);
    return $coderef;

}

sub cut {
    my ( $self, $substrate, $pos ) = @_;

    unless (defined $pos and $pos > 0 and $pos <= length $substrate) {
        carp "Incorrect position.";
        return;
    }

    say $pos;

    $substrate = uc $substrate;
    $substrate = 'XXXX'. $substrate;
    $pos += 4;

    my $pep = substr($substrate, $pos - 4, 8);

    if ($self->_cuts($pep)) {
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

    $substrate = 'XXXX' . $substrate;
    while ( my $pep = substr($substrate, $i, 8) ) {
        if ( $self->_cuts( $pep ) ) {
            my $product = substr($substrate, $j, $i + 4 - $j);
            push @products, $product;
            #    substr($substrate, 0, $i + 4) = '';
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
    my $i = 0;

    $substrate = 'XXXX' . $substrate;
    while ( my $pep = substr($substrate, $i, 8 ) ) {
        if ( $self->_cuts( $pep ) ) { push @sites, $i };
        ++$i;
    }
    return @sites;
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


=head1 DESCRIPTION

This module models the hydrolitic behaviour of a proteolytic enzyme.
Its main purpose is to model or predict the outcome of hydrolitic
cleavage of a peptidic substrate.

The enzyme specificity is currently modeled for 32 enzymes/reagents.
This models are somewhat simplistic as they are largely regex-based, and
do not take into account subtleties such as kinetic/temperature effects,
accessible solvent area or secondary or tertiary structure elements.
However, the module is flexible enough to account for any of these
effects via user-defined functions.

=head1 FUNCTIONS

=cut

=head1 AUTHOR

Bruno Vecchi, C<< <vecchi.b at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-protease at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-Protease>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




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


=head1 ACKNOWLEDGEMENTS


=head1 COPYRIGHT & LICENSE

Copyright 2009 Bruno Vecchi, all rights reserved.

This program is free software; you can redistribute it and/or modify it
under the same terms as Perl itself.


=cut

### Enzyme specificities

BEGIN {
    %specificity_of = (
        'arg-c proteinase'           => [ '.{3}R.{4}' ],
        'asp-n endopeptidase'        => [ '.{4}D.{3}' ],
        'asp-n endopeptidase glu'    => [ '.{4}[DE].{3}' ],
        'bnps skatole'               => [ '.{3}W.{4}' ],
        'caspase 1'                  => [ '[FWYL].[HAT]D[^PEDQKR].{3}' ],
        'caspase 2'                  => [ 'DVAD[^PEDQKR].{3}' ],
        'caspase 3'                  => [ 'DMQD[^PEDQKR].{3}' ],
        'caspase 4'                  => [ 'LEVD[^PEDQKR].{3}' ],
        'caspase 5'                  => [ '[LW]EHD.{4}' ],
        'caspase 6'                  => [ 'VE[HI]D[^PEDQKR].{3}' ],
        'caspase 7'                  => [ 'DEVD[^PEDQKR].{3}' ],
        'caspase 8'                  => [ '[IL]ETD[^PEDQKR].{3}' ],
        'caspase 9'                  => [ 'LEHD.{4}' ],
        'caspase 10'                 => [ 'IEAD.{4}' ],
        'chymotrypsin'               => [ '.{3}[FY][^P].{3}|.{3}W[^MP].{3}' ],
        'chymotrypsin low'           => [ '.{3}[FLY][^P].{3}|.{3}W[^MP].{3}|.{3}M[^PY].{3}|.{3}H[^DMPW].{3}' ],
        'clostripain'                => [ '.{3}R.{4}' ],
        'cnbr'                       => [ '.{3}M.{4}' ],
        'enterokinase'               => [ '[DN][DN][DN]K.{4}' ],
        'factor xa'                  => [ '[AFGILTVM][DE]GR.{4}' ],
        'formic acid'                => [ '.{3}D.{4}' ],
        'glutamyl endopeptidase'     => [ '.{3}E.{4}' ],
        'granzymeb'                  => [ 'IEPD.{4}' ],
        'hydroxylamine'              => [ '.{3}NG.{3}' ],
        'iodosobenzoic acid'         => [ '.{3}W.{4}' ],
        'lysc'                       => [ '.{3}K.{4}' ],
        'lysn'                       => [ '.{4}K.{3}' ],
        'ntcb'                       => [ '.{4}C.{3}' ],
        'pepsin ph1.3'               => [ '.[^HKR][^P][^R][FLWY][^P].{2}|.[^HKR][^P][FLWY].[^P].{2}' ],
        'pepsin'                     => [ '.[^HKR][^P][^R][FL][^P].{2}|.[^HKR][^P][FL].[^P].{2}' ],
        'proline endopeptidase'      => [ '.{2}[HKR]P[^P].{3}' ],
        'proteinase k'               => [ '.{3}[AFILTVWY].{4}' ],
        'staphylococcal peptidase i' => [ '.{2}[^E]E.{4}' ],
        'thermolysin'                => [ '.{3}[^XDE][AFILMV][^P].{2}' ],
        'thrombin'                   => [ '.{2}GRG.{3}|[AFGILTVM][AFGILTVWA]PR[^DE][^DE].{2}' ],
        'trypsin'                    => [ '.{2}(?!CKD).{6}', '.{2}(?!DKD).{6}', '.{2}(?!CKH).{6}', '.{2}(?!CKY).{6}', '.{2}(?!RRH).{6}', '.{2}(?!RRR).{6}', '.{2}(?!CRK).{6}',
                                        '.{3}[KR][^P].{3}|.{2}WKP.{3}|.{2}MRP.{3}' ]
    );
}

no Moose;
__PACKAGE__->meta->make_immutable;
1;
