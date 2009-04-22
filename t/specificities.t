use Modern::Perl;

use lib qw(../lib);
use Test::More qw(no_plan);
use Test::Exception;
use Test::Warn;
use YAML::Any;

use ok 'Bio::Protease';

my $test_seq = <<EOL
mattsfpsmlfyfcifllfhgsmaqlfgqsstpwqssrqgglrgcrfdrlqafeplrqvr
sqagiteyfdeqneqfrctgvsvirrviepqglvlpqyhnapalvyilqgrgftgltfpg
cpatfqqqfqpfdqsqfaqgqsqsqtikdehqrvqrfkqgdvvalpagivhwcyndgdap
ivaiyvfdvnnnanqleprqkkfllagnnkfllagnnanqleprqkefllagnnkreqqs
gnnifsglsvqllsealgisqqaaqgsksndqrgrvirvsqglqflkpivsqqvpveqqv
yqpiqtqdvqatqyqvgqstqyqvgkstpyqggqssqyqagqswdqsfngleenfcslea
rknienpqhadtynpragritrlnsknfpilnivqmsatrvnlyqnailspfwninahsv
iymiqgharvqvvnnngqtvfsdilhrgqllivpqhfvvlknaeregcqyisfktnpnsm
vshiagktsilralpidvlanayrisrqearnlknnrgeefgaftpkltqtgfqsyqdie
easssavraseMVNSNQNQNGNSNGHDDDFPQDSITEPEHMRKLFIGGLDYRTTDENLKA
VMKDPRTKRSRGFGFITYSHSSMIDEAQKSRPHKIDGRVEPKRAVPRQDIDSPNAGATVK
KLFVGALKDDHDEQSIRDYFQHFGNIVDNIVIDKETGKKRGFAFVEFDDYDPVDKVVLQK
QHQLNGKMVDVKKALPKNDQQGGGGGRGGPGGRAGGNRGNMGGGNYGNQNGGGNWNNGGN
NWGNNRGNDNWGNNSFGGGGGGGGGYGGGNNSWGNNNPWDNGNGGGNFGGGGNNWNGGND
FGGYQQNYGGGPQRGGGNFNNNRMQPYQGGGGFKAGGGNQGNYGNNQGFNNGGNNRRYHE
KWGNIVDVVMVNSNQNQNGNSNGHDDDFPQDSITEPEHMRKLFIGGLDYRTTDENLKAHE
VMKDPTSTSTSTSTSTSTSTSTMIDEAQKSRPHKIDGRVEPKRAVPRQDIDSPNAGATVK
KLFVGALKDDHDEQSIRDYFQHLLLLLLLDLLLLDLLLLDLLLFVEFDDYDPVDKVVLQK
QHQLNGKMVDVKKALPKNDQQGGGGGRGGPGGRAGGNRGNMGGGNYGNQNGGGNWNNGGN
NWGNNRGNDNWGNNSFGGGGGGGGGYGGGNNSWGNNNPWDNGNGGGNFGGGGNNWNGGND
FGGYQQNYGGGPQRGGGNFNNNRMQPYQGGGGFKAGGGNQGNYGNNQGFNNGGNNRRYKW
GNIVDVV
EOL
;

$test_seq =~ s/\n//g;
#say length $test_seq;

open( my $fh, '<', 't/specificities.txt' )
    or die "Couldn't open test data file specificities.txt: $!\n";
my $data = join('', <$fh>);
my $true_values = Load $data;

my $results;
my @products;
foreach my $specificity ( @{Bio::Protease->Specificities} ) {
    my $protease = Bio::Protease->new(specificity => $specificity);
    my @cleavage_sites = $protease->cleavage_sites($test_seq);
    $results->{$specificity} =  [scalar @cleavage_sites, [@cleavage_sites] ];
}

is_deeply($results, $true_values);

# -Test cut
my $seq = 'AARAGQTVRFSDAAA';
my $protease = Bio::Protease->new(specificity => 'trypsin');

warning_like { !$protease->cut($seq) } qr/Incorrect position/;

ok !$protease->cut($seq, 1);

is_deeply([ $protease->cut($seq, 3) ], [ 'AAR', 'AGQTVRFSDAAA' ]);
is_deeply([ $protease->cut($seq, 9) ], [ 'AARAGQTVR', 'FSDAAA' ]);

#is_deeply(\@digest_products, \@cut_products);

# test digest
@products = $protease->digest($seq);

is_deeply( [@products], ['AAR', 'AGQTVR', 'FSDAAA'] );
