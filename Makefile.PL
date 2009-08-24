use inc::Module::Install;

name     'Bio-Protease';
all_from 'lib/Bio/Protease.pm';
author   'Bruno Vecchi <vecchi.b@gmail.com>';
license  'perl';

build_requires 'Test::More';
build_requires 'Test::Exception';
requires Moose                   => 0.72;
requires MooseX::ClassAttribute  => 0.08;
requires namespace::autoclean    => 0;
requires Carp                    => 0;
requires Memoize                 => 0;
requires MooseX::Types::Moose     => 0;
requires Modern::Perl            => 0;
requires YAML::Any               => 0;

auto_install;

WriteAll;
