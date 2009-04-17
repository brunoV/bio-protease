package Bio::Protease;
use Modern::Perl;
use Moose;
use MooseX::ClassAttribute;
use Moose::Util::TypeConstraints;

=head1 NAME

Bio::Protease - The great new Bio::Protease!

=head1 VERSION

Version 0.01

=cut

my %specificity_of = (
    'arg-c proteinase'           => \&_argc,
    'asp-n endopeptidase'        => \&_aspn,
    'asp-n endopeptidase glu'    => \&_aspn_glu,
    'bnps skatole'               => \&_bnps,
    'caspase 1'                  => \&_casp1,
    'caspase 2'                  => \&_casp2,
    'caspase 3'                  => \&_casp3,
    'caspase 4'                  => \&_casp4,
    'caspase 5'                  => \&_casp5,
    'caspase 6'                  => \&_casp6,
    'caspase 7'                  => \&_casp7,
    'caspase 8'                  => \&_casp8,
    'caspase 9'                  => \&_casp9,
    'caspase 10'                 => \&_casp2,
    'chymotrypsin'               => \&_chtryp_high,
    'chymotrypsin low'           => \&_chtryp_low,
    'clostripain'                => \&_clostripain,
    'cnbr'                       => \&_cnbr,
    'enterokinase'               => \&_enterokinase,
    'factor xa'                  => \&_factor_xa,
    'formic acid'                => \&_formic_acid,
    'glutamyl endopeptidase'     => \&_glutamil_endopeptidase,
    'granzymeb'                  => \&_granzyme_b,
    'hydroxylamine'              => \&_hydroxylamine,
    'iodosobenzoic acid'         => \&_iodobenzoic_acid,
    'lysc'                       => \&_lys_c,
    'lysn'                       => \&_lys_n,
    'ntcb'                       => \&_ntcb,
    'pepsin ph1.3'               => \&_pepsin_ph_1p3,
    'pepsin'                     => \&_pepsin,
    'proline endopeptidase'      => \&_proline_endopeptidase,
    'proteinase k'               => \&_proteinase_k,
    'staphylococcal peptidase i' => \&_staphilococcal_peptidase_i,
    'thermolysin'                => \&_thermolysin,
    'thrombin'                   => \&_thrombin,
    'trypsin'                    => \&_trypsin,
);

class_has 'Specificities' => (
    is      => 'ro',
    isa     => 'ArrayRef',
    default => sub { [ keys %specificity_of ] },
);

subtype 'Specificity',
    as 'Str',
    where { $_ ~~ @{__PACKAGE__->Specificities()} },
    message {
       "That's not a recognized specificity. Options are: ",
       keys %specificity_of;
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
    build => sub { $specificity_of{(shift)->specificity} },
);

sub _build__cuts {
     my $self = shift;

     return sub {
        my $substrate = shift;

        unless ( length($substrate) == 8 ) {
            my $tail_length = 8 - length($substrate);
            $substrate = $substrate . 'X' x $tail_length;
        }

        my @res = split('', $substrate);
        return $specificity_of{$self->specificity}->(@res);
    }
}

sub digest {
    my ( $self, $substrate ) = @_;

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

    return @products;
}

sub cleavage_sites {
    my ( $self, $substrate ) = @_;
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

sub _argc {
    my @res = @_;
    if ( uc $res[3] eq 'R' ) {
        return 1;
    } else {
        return;
    }
}

sub _aspn {
    my @res = @_;
    if ( uc $res[4] eq 'D' ) {
        return 1;
    } else {
        return;
    }
}

sub _aspn_glu {
    my @res = @_;
    if ( $res[4] =~ /[DE]/i ) {
        return 1;
    } else {
        return;
    }
}

sub _bnps {
    my @res = @_;

    if ( uc $res[3] eq 'W' ) {
        return 1;
    } else {
        return;
    }
}

sub _casp1 {
    my @res = @_;

    if (
        $res[0] !~ /[FWYL]  /ix or
        $res[2] !~ /[HAT]   /ix or
        $res[3] !~ /D       /ix or
        $res[4] =~ /[PEDQKR]/ix
       )   { return   }
    else { return 1 };
}

sub _casp2 {
    my @res = @_;

    if ( 
        $res[0] !~ /D       /ix or
        $res[1] !~ /V       /ix or
        $res[2] !~ /A       /ix or
        $res[3] !~ /D       /ix or
        $res[4] =~ /[PEDQKR]/ix
       )   { return   }
    else { return 1 };
}

sub _casp3 {
    my @res = @_;

    if ( 
        $res[0] !~ /D       /ix or
        $res[1] !~ /M       /ix or
        $res[2] !~ /Q       /ix or
        $res[3] !~ /D       /ix or
        $res[4] =~ /[PEDQKR]/ix
       )   { return   }
    else { return 1 };
}

sub _casp4 {
    my @res = @_;

    if ( 
        $res[0] !~ /L       /ix or
        $res[1] !~ /E       /ix or
        $res[2] !~ /V       /ix or
        $res[3] !~ /D       /ix or
        $res[4] =~ /[PEDQKR]/ix
       )   { return   }
    else { return 1 };
}

sub _casp5 {
    my @res = @_;

    if (
        $res[0] !~ /[LW]    /ix or
        $res[1] !~ /E       /ix or
        $res[2] !~ /H       /ix or
        $res[3] !~ /D       /ix
       )   { return   }
    else { return 1 }
}

sub _casp6 {
    my @res = @_;

    if ( 
        $res[0] !~ /V       /ix or
        $res[1] !~ /E       /ix or
        $res[2] !~ /[HI]    /ix or
        $res[3] !~ /D       /ix or
        $res[4] =~ /[PEDQKR]/ix
       )   { return   }
    else { return 1 };
}

sub _casp7 {
    my @res = @_;

    if (
        $res[0] !~ /D       /ix or
        $res[1] !~ /E       /ix or
        $res[2] !~ /V       /ix or
        $res[3] !~ /D       /ix or
        $res[4] =~ /[PEDQKR]/ix
       )   { return   }
    else { return 1 };
}

sub _casp8 {
    my @res = @_;

    if (
        $res[0] !~ /[IL]    /ix or
        $res[1] !~ /E       /ix or
        $res[2] !~ /T       /ix or
        $res[3] !~ /D       /ix or
        $res[4] =~ /[PEDQKR]/ix
       )   { return   }
    else { return 1 };
}

sub _casp9 {
    my @res = @_;

    if (
        $res[0] !~ /L       /ix or
        $res[1] !~ /E       /ix or
        $res[2] !~ /H       /ix or
        $res[3] !~ /D       /ix
       )   { return   }
    else { return 1 };
}

sub _casp10 {
    my @res = @_;

    if (
        $res[0] !~ /I       /ix or
        $res[1] !~ /E       /ix or
        $res[2] !~ /A       /ix or
        $res[3] !~ /D       /ix
       )   { return   }
    else   { return 1 };
}

sub _chtryp_high {
    my @res = @_;

    if (
        $res[3] =~ /[FY]/ix and
        $res[4] !~ /P   /ix
       )   { return 1 }
    elsif (
        $res[3] =~ /W   /ix and
        $res[4] !~ /[MP]/ix
    )      { return 1 }
    else   { return   };
}

sub _chtryp_low {
    my @res = @_;

    if (
        $res[3] =~ /[FLY]   /ix and
        $res[4] !~ /P       /ix
       )   { return 1 }
    elsif (
        $res[3] =~ /W       /ix and
        $res[4] !~ /[MP]    /ix
    )      { return 1 }
    elsif (
        $res[3] =~ /M       /ix and
        $res[4] !~ /[PY]    /ix
    )      { return 1 }
    elsif (
        $res[3] =~ /H       /ix and
        $res[4] !~ /[DMPW]  /ix
    )      { return 1 }
    else   { return   };
}

sub _clostripain {
    my @res = @_;

    if ($res[3] =~ /R/i)   { return 1 }
    else                   { return   };
}

sub _cnbr {
    my @res = @_;

    if ($res[3] =~ /M/i)   { return 1 }
    else                   { return   };
}

sub _enterokinase {
    my @res = @_;

    if (
        $res[0] =~ /[DN]   /ix and
        $res[1] =~ /[DN]   /ix and
        $res[2] =~ /[DN]   /ix and
        $res[3] =~ /K      /ix
    )      { return 1 }
    else   { return   };
}

sub _factor_xa {
    my @res = @_;

    if (
        $res[3] =~ /R          /ix and
        $res[2] =~ /G          /ix and
        $res[1] =~ /[DE]       /ix and
        $res[0] =~ /[AFGILTVM] /ix
    )      { return 1 }
    else   { return   };
}

sub _formic_acid { # wtf?
    my @res = @_;

    if ($res[3] =~ /D/ix) { return 1 }
    else                  { return   };
}

sub _glutamil_endopeptidase {
    my @res = @_;

    if ($res[3] =~ /E/ix) { return 1 }
    else                  { return   };
}

sub _granzyme_b {
    my @res = @_;

    if (
        $res[0] =~ /I /ix and
        $res[1] =~ /E /ix and
        $res[2] =~ /P /ix and
        $res[3] =~ /D /ix
    )      { return 1 }
    else   { return   };
}

sub _hydroxylamine {
    my @res = @_;

    if (
        $res[3] =~ /N /ix and
        $res[4] =~ /G /ix
    )      { return 1 }
    else   { return   };
}

sub _iodobenzoic_acid {
    my @res = @_;

    if ($res[3] =~ /W/ix) { return 1 }
    else                  { return   };
}

sub _lys_c {
    my @res = @_;

    if ($res[3] =~ /K/ix) { return 1 }
    else                  { return   };
}

sub _lys_n {
    my @res = @_;

    if ($res[4] =~ /K/ix) { return 1 }
    else                  { return   };
}

sub _ntcb {
    my @res = @_;

    if ($res[4] =~ /C/ix) { return 1 }
    else                  { return   };
}

sub _pepsin_ph_1p3 {
    my @res = @_;

    if (
        $res[1] !~ /[HKR] /ix and
        $res[2] !~ /P     /ix and
        $res[3] !~ /R     /ix and
        $res[4] =~ /[FLWY]/ix and
        $res[5] !~ /P     /ix
       )   { return 1 }
    elsif (
        $res[1] !~ /[HKR] /ix and
        $res[2] !~ /P     /ix and
        $res[3] =~ /[FLWY]/ix and
        $res[5] !~ /P     /ix
    )      { return 1 }
    else   { return   };
}

sub _pepsin { # pH more than 2
    my @res = @_;

    if (
        $res[1] !~ /[HKR] /ix and
        $res[2] !~ /P     /ix and
        $res[3] !~ /R     /ix and
        $res[4] =~ /[FL]  /ix and
        $res[5] !~ /P     /ix
       )   { return 1 }
    elsif (
        $res[1] !~ /[HKR] /ix and
        $res[2] !~ /P     /ix and
        $res[3] =~ /[FL]  /ix and
        $res[5] !~ /P     /ix
    )      { return 1 }
    else   { return   };
}

sub _proline_endopeptidase {
    my @res = @_;

    if (
        $res[2] =~ /[HKR]/ix and
        $res[3] =~ /P    /ix and
        $res[4] !~ /P    /ix
    )      { return 1 }
    else   { return   };
}

sub _proteinase_k {
    my @res = @_;

    if ($res[3] =~ /[AFILTVWY]/ix) { return 1 }
    else                            { return   };
}

sub _staphilococcal_peptidase_i {
    my @res = @_;

    if (
        $res[3] =~ /E/i and
        $res[2] !~ /E/i
    )      { return 1 }
    else   { return   };
}

sub _thermolysin {
    my @res = @_;

    if (
        $res[4] =~ /[AFILMV]/ix and
        $res[3] =~ /[^DEX]  /ix and
        $res[5] !~ /P       /ix

    )      { return 1 }
    else   { return   };
}

sub _thrombin {
    my @res = @_;

    if (
        $res[2] =~ /G          /ix and
        $res[3] =~ /R          /ix and
        $res[4] =~ /G          /ix
       )   { return 1 }
    elsif       (
        $res[2] =~ /P          /ix and
        $res[3] =~ /R          /ix and
        $res[0] =~ /[AFGILTVM] /ix and
        $res[1] =~ /[AFGILTVWA]/ix and
        $res[4] !~ /[DE]       /ix and
        $res[5] !~ /[DE]       /ix
    )      { return 1 }
    else   { return   };
}

sub _trypsin {

    my @res = @_;

    if (
        $res[2] =~ /[CD] /ix and
        $res[3] =~ /K    /ix and
        $res[4] =~ /D    /ix
       )   { return }
    elsif       (
        $res[2] =~ /C    /ix and
        $res[3] =~ /K    /ix and
        $res[4] =~ /[HY] /ix
    )      { return }
    elsif       (
        $res[2] =~ /C    /ix and
        $res[3] =~ /R    /ix and
        $res[4] =~ /K    /ix
    )      { return }
    elsif       (
        $res[2] =~ /R    /ix and
        $res[3] =~ /R    /ix and
        $res[4] =~ /[HR] /ix
    )      { return }
    elsif       (
        $res[3] =~ /[KR] /ix and
        $res[4] !~ /P    /ix
    )      { return 1 }
    elsif       (
        $res[2] =~ /W    /ix and
        $res[3] =~ /K    /ix and
        $res[4] =~ /P    /ix
    )      { return 1 }
    elsif       (
        $res[2] =~ /M   /ix and
        $res[3] =~ /R   /ix and
        $res[4] =~ /P   /ix
    )      { return 1 }
    else   { return   };
}

no Moose;
__PACKAGE__->meta->make_immutable;
1;
