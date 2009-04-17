#!perl -T

use Test::More tests => 1;

BEGIN {
	use_ok( 'Bio::Protease' );
}

diag( "Testing Bio::Protease $Bio::Protease::VERSION, Perl $], $^X" );
