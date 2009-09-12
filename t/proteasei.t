use Test::More;
use Modern::Perl;

{
    package My::Protease;
    use Moose;
    extends qw(Bio::ProteaseI);

    augment _cuts => sub {
        my ( $self, $subs_ref ) = @_;

        if ( $$subs_ref eq 'MAELVIKP' ) { return 1 }
        else                            { return   }
    };

}

my $protease = My::Protease->new;

isa_ok( $protease, 'My::Protease' );

my @products = $protease->digest( 'AAAAMAELVIKPYYYYYYY' );

is_deeply(
    \@products, ["AAAAMAEL", "VIKPYYYYYYY"],
    "Subclassing works as expected"
);

done_testing();
