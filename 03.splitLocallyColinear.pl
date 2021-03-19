#!/usr/bin/perl -w
=head1 Description

This script is used to make sure the homolog group coordinates with the order of the left to right, then split the crossed sequences to different homolog group so that make several locally colinear blocks.

=head1 Contact & Version

  Author: Qi Fang, fangqi@genomics.cn
  Version: 1.0,  Date: 2019-02-14
  Version: 2.0,  Date: 2019-03-29, Modified: change the @c to @{$hash{$split_tag}} so that the coordinates would be correct among more than one splits.
  Version: 2.1,  Date: 2019-03-30, Modified: change the line45 from $seq to $seg{$split_tag} so that there no wrong dump with the sequence splits.
  Version: 2.2,  Date: 2019-12-10, Modified: Debug the coodinates redundancy of merging neighboring blocks when split the INDEL out which was reported by YuanDeng.
  Version: 2.3,  Date: 2020-02-28, Modified: Debug the split length ($s) error when the alignments could be merged.
  Version: 3.0,  Date: 2020-09-22, Modified: Further debug the split length ($s) error when the alignments could be merged.

=head1 Usage Exmples

  03.splitLocallyColinear.pl Anc001.maf.sort.all.fasta > Anc001.maf.sort.all.s.fasta

=cut

use strict;
my $input = shift;
open IN,$input;
$/=">";<IN>;
while(<IN>){
    my (%hash,%seg);
    my ($split_tag,$tag,$s)=(1,0,0);
    chomp;
    my @a = split /\s+/;
    my $seq = $a[3];
    my @c = split(/,|;/,$a[1]);
    $seg{$split_tag}=$seq;
    if($#c==1){
        push @{$hash{$split_tag}}, ($c[0],$c[1]);
    }else{
        for(my $i=0;$i<=$#c-2;$i+=2){
            if($c[$i+2]>$c[$i+1]){
                push @{$hash{$split_tag}}, ($c[$i],$c[$i+1]);
                $s+=$c[$i+1]-$c[$i];
            }elsif($c[$i+2]==$c[$i+1]){
                push @{$hash{$split_tag}}, ($c[$i],$c[$i+3]);
                $s+=$c[$i+3]-$c[$i];$i+=2;
            }else{
                $tag++;
                my $sT=0; my $sE=0;
                unshift @{$hash{$split_tag+$tag}}, ($c[$i],$c[$i+1]);
                $sT+=$c[$i+1]-$c[$i];
                for(my $j=$#{$hash{$split_tag}};$j>0;$j-=2){
                    if($hash{$split_tag}[$j]<$c[$i+2]){
                        last;
                    }elsif($hash{$split_tag}[$j]==$c[$i+2]){ 
                        $hash{$split_tag}[-1]=$c[$i+3];
                        $sE=$c[$i+3]-$c[$i+2]; ##debug by fangqi at 20200228; debug by fangqi at 20200922: if add $s here it would disturb the sequences substr in this run
                        $i+=2; ##debug by fangqi at 20191210
                        last;
                    }else{
                        unshift @{$hash{$split_tag+$tag}}, ($hash{$split_tag}[$j-1],$hash{$split_tag}[$j]);
                        $sT+=$hash{$split_tag}[$j]-$hash{$split_tag}[$j-1];
                        $s-=$hash{$split_tag}[$j]-$hash{$split_tag}[$j-1];
                        pop @{$hash{$split_tag}};
                        pop @{$hash{$split_tag}};
                    }
                }
                my $sN=0;
                for(my $k=0;$k<length($seq);$k++){
                    my $nu = substr($seg{$split_tag},$k,1);
                    if($nu ne "-"){
                        $sN++;
                        if($sN>$s and $sN<=$s+$sT){
                            $seg{$split_tag+$tag}.=$nu;
                            substr($seg{$split_tag}, $k, 1, "-");
                        }else{
                            $seg{$split_tag+$tag}.="-";
                        }
                    }else{
                        $seg{$split_tag+$tag}.="-";
                   }
                }
                $s+=$sE;
            }
        }
        push @{$hash{$split_tag}}, ($c[-2],$c[-1]);
    }
    for(my $k=0;$k<=$tag;$k++){
        if($tag>0){
            my $tmpTag = $split_tag+$k;
            print ">$a[0]_p$tmpTag\t";
        }else{
            print ">$a[0]\t";
        }
        my ($pos,$len);
        for(my $i=0;$i<=$#{$hash{$split_tag+$k}};$i+=2){
            $pos.=$hash{$split_tag+$k}[$i].",".$hash{$split_tag+$k}[$i+1].";";
            $len.=($hash{$split_tag+$k}[$i+1]-$hash{$split_tag+$k}[$i]).",";
        }
        if(!defined $pos){print STDERR "$a[0]\n";}
        chop $pos;
        chop $len;
        print "$pos\t$len\n$seg{$split_tag+$k}\n";
    }
}close IN;
