#!/usr/bin/perl
use Modern::Perl;
use YAML::Any;

my $cut_probability = {
    A => { A => 0.5 },
};

my $seq = 'AAAAAAAAAA';

my @products = digest($seq);

say Dump(@products);

sub digest {
    my ( $substrate ) = @_;
    my @monomers = split('', $substrate);
    my @products;
    my $i = 0;

    while ( exists $monomers[$i+1] ) {
        if ( _cut ( $monomers[$i], $monomers[$i+1] ) ) {
            push @products, splice(@monomers, 0, $i+1);
            $i = 0;
        }
        else { ++$i }
    }

    return @products;
}

sub _cut {
    my ($left, $right) = @_;
    my $prob = $cut_probability->{$left}{$right};

    if ( rand() >= $prob ) {
        return 1;
    }

    else { return }
}


