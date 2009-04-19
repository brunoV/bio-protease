use Modern::Perl;

use lib qw(../lib);
use Test::More qw(no_plan);
use Test::Exception;

use ok 'Bio::Protease';

my $enzyme;

lives_ok { $enzyme = Bio::Protease->new(specificity => 'trypsin') };

isa_ok($enzyme, 'Bio::Protease');

dies_ok { $enzyme = Bio::Protease->new() };
dies_ok { $enzyme = Bio::Protease->new(specificity => 'ooo') };

lives_ok { $enzyme = Bio::Protease->new(specificity => sub { return 1 }) };

use Bio::Tools::SeqPattern;
my $pattern = Bio::Tools::SeqPattern->new(
    -SEQ => 'AGGAL[^P]', -TYPE => 'Amino'
);

lives_ok { $enzyme = Bio::Protease->new(specificity => $pattern) };
