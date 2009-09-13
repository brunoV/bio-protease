package Bio::ProteaseI;

# ABSTRACT: A base class to build your customized Protease

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

=head1 SYNOPSIS

    package My::Protease;
    use Moose;
    extends qw(Bio::ProteaseI);

    augment _cuts => sub {
        my ($self, $substrate) = @_;

        # some code that decides
        # if $peptide should be cut or not

        if ( $peptide_should_be_cut ) { return 1 }
        else                          { return   }
    };

=head1 DESCRIPTION

This module describes the interface for L<Bio::Protease>. You only need
to use this if you want to build your custom specificity protease and
regular expressions won't do; otherwise look at L<Bio::Protease>
instead.

=head1 METHODS

All of the methods provided in Bio::Protease (namely, C<cut>, C<digest>,
C<is_substrate> and C<cleavage_sites>) are defined here, incluiding a
stub of the specificity-determining one, C<_cuts>. It has to be
completed by the subclass with an C<augment> call.

=head1 HOW TO SUBCLASS

=head2 Step 1: create a child class.

    package My::Protease;
    use Moose;
    extends qw(Bio::ProteaseI);

    1;

Simply create a new Moose class, and inherit from the Bio::ProteaseI
interfase using C<extends>.

=head2 Step 2: augment _cuts()

The C<_cuts> subroutine will be used by the methods C<digest>, C<cut>,
C<cleavage_sites> and C<is_substrate>. It will B<always> be passed a
reference to a string of 8 characters; if the subroutine returns true,
then the peptide bond between the 4th and 5th residues will be marked as
siscile, and the appropiate action will be performed depending on which
method was called.

Your specificity logic should only be concerned in deciding whether the
8-residue long peptide passed to it as an argument should be cut between
the 4th and 5th residues. This is done by using the C<augment> method
modifier (for more information on Method Modifiers, please read up on
L<Moose::Manual::MethodModifiers>), like so:

    augment _cuts => sub {
        my ( $self, $peptide ) = @_;

        # some code that decides
        # if $peptide should be cut or not

        if ( $peptide_should_be_cut ) { return 1 }
        else                          { return   }
    };

And that's it. Your class will inherit all the methods mentioned above,
and will work according to the specificity logic that you define in your
C<_cuts()> subroutine.

=head2 Example: a ridiculously specific protease.

Suppose you want to model a protease that only cleaves the sequence
C<MAEL^VIKP>. Your Protease class would be like this:

    package My::Ridiculously::Specific::Protease;
    use Moose;
    extends qw(Bio::ProteaseI);

    augment _cuts => sub {
        my ( $self, $substrate_ref ) = @_;

        if ( $$substrate_ref eq 'MAELVIKP' ) { return 1 }
        else                                 { return   }
    };

    1;

Then you can use your class easily in your application:

    #!/usr/bin/env perl
    use Modern::Perl;

    use My::Ridiculously::Specific::Protease;

    my $protease = My::Ridiculously::Specific::Protease->new;
    my @products = $protease->digest( 'AAAAMAELVIKPYYYYYYY' );

    say for @products; # ["AAAAMAEL", "VIKPYYYYYYY"]

=cut
