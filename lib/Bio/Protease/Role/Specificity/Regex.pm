package Bio::Protease::Role::Specificity::Regex;
use Moose::Role;
use Bio::Protease::Types 'ProteaseRegex';

has regex => (
    is  => 'ro',
    isa => ProteaseRegex,
    coerce => 1,
);

sub _cuts {
    my ($self, $peptide) = @_;

    if ( grep { $peptide !~ /$_/ } @{$self->regex} ) {
        return;
    }

    return 'yes, it cuts';
}

1;
