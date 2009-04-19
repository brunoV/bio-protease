use inc::Module::Install;

name     'Bio-Protease';
all_from 'lib/Bio/Protease.pm';
author   'Bruno Vecchi <vecchi.b@gmail.com>';
license  'perl';

build_requires 'Test::More';
requires Moose                     => 0.72;
requires 'MooseX::ClassAttribute'  => 0.08;

auto_install;

WriteAll;
