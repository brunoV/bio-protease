use Modern::Perl;

use lib qw(../lib);
use Test::More qw(no_plan);
use Test::Exception;

use ok 'Bio::Protease';

my $enzyme;

lives_ok { $enzyme = Bio::Protease->new(specificity => 'trypsin') };

isa_ok($enzyme, 'Bio::Protease');

#is($enzyme->specificity, 'trypsin');

dies_ok { $enzyme = Bio::Protease->new() };
dies_ok { $enzyme = Bio::Protease->new(specificity => 'ooo') };
