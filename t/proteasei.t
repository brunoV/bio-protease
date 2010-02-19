use Test::More;
use Modern::Perl;

{
    package My::Protease;
    use Moose;
    with qw(Bio::ProteaseI);

    sub _cuts {
        my ( $self, $substrate ) = @_;

        if ( $substrate eq 'MAELVIKP' ) { return 1 }
        else                            { return   }
    };

}

my $protease = My::Protease->new;

isa_ok( $protease, 'My::Protease' );

can_ok( $protease, qw(cut is_substrate digest cleavage_sites) );

ok !$protease->cut('foo', -1), "no cutting in senseless positions";

my @products = $protease->digest( 'AAAAMAELVIKPYYYYYYY' );

is_deeply(
    \@products, ["AAAAMAEL", "VIKPYYYYYYY"],
    "Subclassing works as expected"
);

done_testing();
