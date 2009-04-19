package Bio::Protease;
use Modern::Perl;
use Moose;
use MooseX::ClassAttribute;
use Moose::Util::TypeConstraints;
use Carp;

my %specificity_of;

=head1 NAME

Bio::Protease - The great new Bio::Protease!

=head1 VERSION

Version 0.01

=cut

class_has 'Specificities' => (
    is      => 'ro',
    isa     => 'ArrayRef',
    default => sub { [ keys %specificity_of ] },
);

subtype 'Specificity',
    as 'Str',
    where { $_ ~~ @{__PACKAGE__->Specificities()} },
    message {
       "That's not a recognized specificity.
       Options are: @{[keys %specificity_of]}";
    };

has specificity => (
   is         => 'ro',
   isa        => 'Specificity',
   required   => 1,
   coerce     => 1,
);

has _cuts => (
    is          => 'ro',
    isa         => 'CodeRef',
    lazy_build  => 1,
);

sub _build__cuts {
    my $self = shift;
    return sub {
        my $substrate = shift;
        unless ( length $substrate == 8 ) {
            $substrate .= 'X' x (8 - length $substrate);
        }

        my @regexes = @{$specificity_of{$self->specificity}};

        if ( grep { $substrate !~ /$_/ } @regexes ) {
            return;
        } else {
            return 1;
        }
    }
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

    if ($self->_cuts->($pep)) {
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
        if ( $self->_cuts->( $pep ) ) {
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
        if ( $self->_cuts->( $pep ) ) { push @sites, $i };
        ++$i;
    }
    return @sites;
}

our $VERSION = '0.01';

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use Bio::Protease;

    my $foo = Bio::Protease->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 FUNCTIONS

=head2 function1

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
