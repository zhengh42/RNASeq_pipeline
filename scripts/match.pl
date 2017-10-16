#! /usr/bin/perl -w
use strict;

die "perl $0 file1 col1 file2 col2" unless @ARGV >= 4;

my $n1 = $ARGV[1] - 1 ;
my $n2 = $ARGV[3] - 1;

my %eval;
open(IN,$ARGV[0]) or die;
while(<IN>){
        chomp;
        my @t = split /\t/;
        $eval{$t[$n1]} = $_;
}
close IN;

open(INA,$ARGV[2]) or die;
while(<INA>){
        chomp;
        my @t = split /\t/;
        if(exists($eval{$t[$n2]})){print "$eval{$t[$n2]}\n";}
        else{print "$t[$n2]\tna\n";}
}
close INA;

