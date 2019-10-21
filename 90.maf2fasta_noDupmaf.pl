#!/usr/bin/perl -w

=head1 Description

This script is used to transform the noDup maf from Cactus to fasta format.

=head1 Contact & Version

  Author: Qi Fang, fangqi@genomics.cn
  Version: 1.0,  Date: 2019-09-09

=head1 Command-line Option

  --sp <str>    list of species (the same as the species name in guide tree), if not, the guidetree would be used and don't output the ancestors
  --input <str> input MAF format file

  --help        output help information to screen

=head1 Usage Exmples

  90.maf2fasta_nofiltering.pl -i test.maf -s species.list > test.fasta

=cut

use strict;
use Getopt::Long;

my ($Input,$Name,$Help);
GetOptions(
        "sp:s"=>\$Name,
        "input:s"=>\$Input,
        "help"=>\$Help
);

die `pod2text $0` if ($Help or !$Input);

my %hash;
if($Name){
    open IN,$Name;
    while(<IN>){
    	chomp;
    	my @a = split /\s+/;
    	$hash{$a[0]}=1;
    }close IN;
}else{
    open IN,$Input;
    while(<IN>){
        next, unless($_=~/^# hal/);
        my @a = split /\(|\)|:|,|;|\s+/;
        foreach my $s (@a){
            if($s and $s!~/[0-9]/ and $s ne "hal" and $s ne "#"){
                $hash{$s}=1;
            }
        }
        last;
    }
    close IN;
}

my (%seq,%len);my $last_len=0;
open IN,$Input;
$/ = "\n\n";
while(<IN>){
	next, if($_=~/#/);
	chomp;
	my @a = split /\s+/;
    my $start;
    for(my $i=1;$i<=$#a;$i++){
        if($a[$i] eq "s"){
            $start = $i+1;
            last;
        }
    }
    foreach my $key (sort keys %hash){
        my $d=0;
        for(my $i=$start;$i<=$#a;$i+=7){
            my $name = $1, if($a[$i]=~/(\w+)\.(\S+)/);
            my $chr=$2;
            if($name eq $key){
                if(!exists $seq{$name}){
                    for(my $j=1;$j<=$last_len;$j++){
                        $seq{$name} .= "-";
                    }
                    $seq{$name} .= $a[$i+5];
                }else{$seq{$name} .= $a[$i+5];}

                $len{$name} += $a[$i+2];
                $d=1;
                last;
            }
        }
        if($d == 0){
            my $gap = length($a[7]);
            for(my $j=1;$j<=$gap;$j++){
                last, if(!exists $seq{$key});
                $seq{$key} .= "-";
            }
        }
    }
    $last_len+=length($a[7]);
}
close IN;

foreach my $key (sort keys %hash){
	next, unless(exists $len{$key} and $len{$key}>0);
	print ">$key\n$seq{$key}\n";
}
