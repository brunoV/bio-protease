use Modern::Perl;

use lib qw(../lib);
use Test::More qw(no_plan);
use Test::Exception;

use ok 'Bio::Protease';

my $enzyme;

lives_ok { $enzyme = Bio::Protease->new(specificity => 'trypsin') };

is $enzyme->specificity, 'trypsin';

isa_ok($enzyme, 'Bio::Protease');

dies_ok { $enzyme = Bio::Protease->new() };
dies_ok { $enzyme = Bio::Protease->new(specificity => 'ooo') };

my $pattern = ['AGGAL[^P]'];

lives_ok { $enzyme = Bio::Protease->new(specificity => $pattern) };

is $enzyme->specificity, 'custom';
