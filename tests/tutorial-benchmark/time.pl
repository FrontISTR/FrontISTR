#!/usr/bin/perl
use strict;
use warnings;
use Time::HiRes qw(time);

my $output = shift @ARGV;
my $start = time();
my $status = system { $ARGV[0] } @ARGV;
my $elapsed = time() - $start;
open my $file, '>', $output or die "cannot write $output: $!\n";
print {$file} "$elapsed\n";
close $file;
exit($status == -1 ? 127 : $status & 127 ? 128 + ($status & 127) : $status >> 8);
