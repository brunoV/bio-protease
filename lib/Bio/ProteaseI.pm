package Bio::ProteaseI;
use Moose;
use Carp;
use Memoize qw(memoize flush_cache);
use namespace::autoclean;

memoize ('cleavage_sites');
memoize ('is_substrate');
memoize ('digest');

sub cut {
    my ( $self, $substrate, $pos ) = @_;

    unless ( defined $pos and $pos > 0 and $pos <= length $substrate ) {

        carp "Incorrect position.";
        return;
    }

    $substrate = uc $substrate;
    $substrate = 'XXX'. $substrate;
    $pos += 3;

    my $pep = substr($substrate, $pos - 4, 8);

    if ( $self->_cuts(\$pep) ) {
        my $product = substr($substrate, 0, $pos);
        substr($substrate, 0, $pos) = '';

        s/X//g for ($product, $substrate);

        return ($product, $substrate);
    }

    else { return }
}

sub digest {
    my ( $self, $substrate ) = @_;
    $substrate = uc $substrate;
    my @products;
    my ($i, $j) = (0, 0);

    $substrate = 'XXX' . $substrate;
    while ( my $pep = substr($substrate, $i, 8) ) {
        if ( $self->_cuts( \$pep ) ) {
            my $product = substr($substrate, $j, $i + 4 - $j);
            push @products, $product;

            $j = $i + 4;
        }
        $i++;
    }
    push @products, substr($substrate, $j - length($substrate));

    s/X//g for @products[0, -1];

    return @products;
}

sub is_substrate {
    my ($self, $substrate) = @_;

    for my $pos (1 .. length $substrate) {
        return 1 if $self->cut($substrate, $pos);
    }

    return;
}

sub _cuts {

    # Substrate needs to be passed by reference in order to modify the
    # argument for the inner() call. Otherwise, all modifications to
    # $substrate would only be local to this sub, and a fresh unmodified
    # copy would be given to inner(), giving unwanted results.

    my ($self, $substrate) = @_;
    my $length = length $$substrate;
    if ( $length < 8 ) {
        if ( $length > 4 ) {
            $$substrate .= 'X' x (8 - $length);
        }
        else { return }
    }

    inner();

}

sub cleavage_sites {
    my ( $self, $substrate ) = @_;
    $substrate = uc $substrate;
    my @sites;
    my $i = 1;

    $substrate = 'XXX' . $substrate;
    while ( my $pep = substr($substrate, $i-1, 8 ) ) {
        if ( $self->_cuts( \$pep ) ) { push @sites, $i };
        ++$i;
    }
    return @sites;
}

sub DEMOLISH {
    flush_cache('digest');
    flush_cache('cleavage_sites');
    flush_cache('is_substrate');
}

__PACKAGE__->meta->make_immutable;

__END__

=pod

=head1 NAME

Bio::ProteaseI - A base class to build your customized Protease

=head1 SYNOPSIS

    package My::Protease;
    use Moose;
    extends qw(Bio::ProteaseI);

    augment _cuts => sub {
        my $substrate = shift;

        # some code that decides
        # if $peptide should be cut or not

        if ( $peptide_should_be_cut ) { return 1 }
        else                          { return   }
    };

=head1 DESCRIPTION

This module describes the interface for Bio::Protease. You only need to
use this if you want to build your custom specificity protease and
regular expressions won't do; otherwise look at Bio::Protease instead.

=head1 METHODS

All of the methods provided in Bio::Protease (namely, cut, digest,
is_substrate and cleavage_sites) are defined here, incluiding a stub of
the specificity-determining one, '_cuts'. It has to be completed by the
subclass with an 'augment' call.

=head1 HOW TO SUBCLASS

This subroutine will be used by the methods C<digest>, C<cut>,
C<cleavage_sites> and C<is_substrate> It will always be passed a string
with a length of 8 characters; if the subroutine returns true, then the
peptide bond between the 4th and 5th residues will be marked as siscile,
and the appropiate action will be performed depending on which method
was called

TODO

=cut
