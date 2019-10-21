#!/usr/bin/perl -w

=head1 Description

This script is used to exchange the sequences name from cactus to some unique ones.

=head1 Contact & Version

  Author: Qi Fang, fangqi@genomics.cn
  Version: 1.0,  Date: 2019-03-29

=head1 Command-line Option

  --name <str>      list of name list to correct the sequence names, default: 363birds_correct_name.table.gz
  --input <str>     input MAF format file

  --help            output help information to screen

=head1 Usage Exmples

  00.maf_correct_rename.pl -i Anc003refChr182.bed.maf -n 363birds_correct_name.table.gz > Anc003refChr182.bed.maf.N

=cut

use strict;
use Getopt::Long;

my ($Input,$Name,$Help);
GetOptions(
        "name:s"=>\$Name,
        "input:s"=>\$Input,
        "help"=>\$Help
);

die `pod2text $0` if ($Help);

$Name ||= "/hwfssz5/ST_DIVERSITY/B10K/PUB/USER/fangqi/share/Cactus_Alignments_Tools/363birds_correct_name.table.gz";

open IN,"gzip -dc $Name | ";
my %list;
while(<IN>){
    chomp;
    my @a = split /\t+/;
    $list{$a[1]}=$a[0];
}close IN;

open IN,$Input;
$/ = "\n\n";
while(<IN>){
    chomp;
    if($_=~/^#/){
        print "$_\n\n";
        next;
    }
    $_=~s/^a .*/a/;
    my @a = split /\s+/;
    print "$a[0]\n";
    for(my $i=2;$i<$#a;$i+=7){
        if(exists $list{$a[$i]}){
            $a[$i]=$list{$a[$i]};
        }
        print "$a[$i-1]\t$a[$i]\t$a[$i+1]\t$a[$i+2]\t$a[$i+3]\t$a[$i+4]\t$a[$i+5]\n";
    }
    print "\n";
}close IN;
