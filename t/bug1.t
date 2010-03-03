#!/usr/bin/env perl
use Test::More;
use Test::Exception;

use Bio::Protease;

lives_ok { $protease = Bio::Protease->new( specificity => qr/.{3}[GA].{4}/ ) };

done_testing();
