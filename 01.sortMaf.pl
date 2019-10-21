#!/usr/bin/perl -w

=head1 Description

This script is used to sort the input MAF format multiple alignments by sequence names and their strands.

=head1 Contact & Version

  Author: Qi Fang, fangqi@genomics.cn
  Version: 1.0,  Date: 2018-12-20
  Version: 1.1,  Date: 2019-02-03 Add the support to multiple hits of reference

=head1 Command-line Option

  --ref <str>       reference of input multiple alignments
  --input <str>     input MAF format file

  --help            output help information to screen

=head1 Usage Exmples

  01.sortMaf.pl -i Anc001.maf -r Anc001 > Anc001.bed.maf.sort

=cut

use strict;
use Getopt::Long;

my ($Input,$Ref,$Help);
GetOptions(
        "ref:s"=>\$Ref,
        "input:s"=>\$Input,
        "help"=>\$Help
);

die `pod2text $0` if ($Help);

open IN,$Input;
$/ = "\n\n";
while(<IN>){
    chomp;
    if($_=~/^#/){
        print "$_\n\n";
        next;
    }
    $_=~s/^a .*/a/;
    my %hash;
    my @a = split /\s+/;
###make reference be the first 
    $hash{"000".$a[2]}=$a[2]."\t".$a[3]."\t".$a[4]."\t".$a[5]."\t".$a[6]."\t".$a[7];
    for(my $i=9; $i<=$#a; $i+=7){
        my $new_name=&get_strand($a[$i],$a[$i+3]);
        $new_name.=$a[$i+1];
        $hash{$new_name}=$a[$i]."\t".$a[$i+1]."\t".$a[$i+2]."\t".$a[$i+3]."\t".$a[$i+4]."\t".$a[$i+5];
    }

    my $reftag=0;
    foreach my $key (sort keys %hash){
        if($key=~/$Ref/){
            if($reftag==0){
                print "a\ns\t$hash{$key}\n";
                $reftag=1;
            }else{print "s\t$hash{$key}\n";}
        }
    }

    foreach my $key (sort keys %hash){
        unless($key=~/$Ref/){
            print "s\t$hash{$key}\n";
        }
    }

    print "\n";
}close IN;


######################################################################
##consider about the strand, idea from lifang
######################################################################
sub get_strand{
    my ($name,$strand) = @_;
    if($strand eq "+"){
        return $name."_P";
    }else{
        return $name."_M";
    }
}
