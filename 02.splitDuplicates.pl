#!/usr/bin/perl -w

=head1 Description

This script is used to connect the homology blocks to homology groups:
a)  The aligned blocks from the same scaffold and same strand would be homology blocks;
b)  The former and later homology blocks having the smallest pairwise distance would be in the same homology groups;
c)  At last, we would get many connected homology blocks (result as one homology group for one connected homology blocks) which all aligned to the reference.

=head1 Contact & Version

  Author: Qi Fang, fangqi@genomics.cn
  Version: 1.0,  Date: 2018-12-20
  Version: 1.1,  Date: 2019-02-03 Add the support to multiple hits of reference
  Version: 2.0,  Date: 2019-09-11 Modified: add the auto-detecting with species list from hal2maf output (the second line: guide tree)

=head1 Command-line Option

  --species_list <str>      list of species needed in results (the same as the species name in guide tree), if not, the guidetree would be used to extracted without ancestors. Notice: species name should not include any number!
  --input <str>             input MAF format file, should be sorted by 01.sortMaf.pl

  --help                    output help information to screen

=head1 Usage Exmples

  02.splitDuplicates.pl -i Anc001.maf.sort -s hal_species.list > Anc001.maf.sort.all.fasta

=cut

use strict;
use Getopt::Long;

my ($Input,$List,$Help);
GetOptions(
        "species_list:s"=>\$List,
        "input:s"=>\$Input,
        "help"=>\$Help
);

die `pod2text $0` if ($Help);

my %hash;
if($List){
##only get the listed species
    open IN,$List;
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

##get para number
my (%para,$last,$last_tag);
open IN,$Input;
while(<IN>){
    next, unless($_=~/^s/);
    chomp;
    my @a = split /\s+/;

#consider about the strand
    my $name_strand = &get_strand($a[1],$a[4]);

#get hash for all paralog-like
    if(!$last or $last ne $name_strand){
        $last_tag=1;$last=$name_strand;
    }
    elsif($last eq $name_strand and $last_tag >= 1){
        $last_tag++; 
        $para{$name_strand}=1, unless(exists $para{$name_strand}); 
        $para{$name_strand} = $last_tag>$para{$name_strand} ? $last_tag : $para{$name_strand};
        $last=$name_strand;
    }
}close IN;

##main
my (%seq,%len,%pos,%max_para);
my $last_len=0;
open IN,$Input;
$/ = "\n\n";
while(<IN>){
	next, if($_=~/#/);
	chomp;
	my @a = split /\s+/;
##get each seq
    for(my $i=2;$i<=$#a;$i+=7){
        if($i==2){
            $a[$i]=$a[$i]."_A";
        }else{
            $a[$i]=&get_strand($a[$i],$a[$i+3]);
        }
        my ($name,$chr) = ($1,$2), if($a[$i]=~/(\w+)\.(\S+)/);
        if(exists $hash{$name}){
##para having
            if(exists $para{$a[$i]}){
##init pos
                if(!exists $pos{$name}{$chr}{1}){
                    push @{$pos{$name}{$chr}{1}}, ($a[$i+1],$a[$i+1]+$a[$i+2]);
                    for(my $j=1;$j<=$last_len;$j++){$seq{$name}{$chr."_p1"} .= "-";}
                    $seq{$name}{$chr."_p1"} .= $a[$i+5];
                    $max_para{$a[$i]}++;
                    my $a_i = $a[$i];
                    $a[$i]=$name.".".$chr."_p1"; ##rename the chr
                    my $time3=0;
                    for(my $time=2;$time<=$para{$a_i};$time++){
                        my $time2 = ($time-1) * 7;
                        if($a[$i+$time2]){
                            my $tmp_name = &get_strand($a[$i+$time2],$a[$i+3+$time2]);
                            if($a_i eq $tmp_name){
                                push @{$pos{$name}{$chr}{$time}}, ($a[$i+1+$time2],$a[$i+1+$time2]+$a[$i+2+$time2]);
                                if(!exists $seq{$name}{$chr."_p".$time}){
                                    for(my $j=1;$j<=$last_len;$j++){
                                        $seq{$name}{$chr."_p".$time} .= "-";
                                    }
                                }
                                $seq{$name}{$chr."_p".$time} .= $a[$i+5+$time2];
                                $max_para{$a_i}++;
                                $a[$i+$time2]=$name."\.".$chr."_p".$time; ##rename the chr
                                $time3 = $time2;
                            }
                        }
                    }
                    $i+=$time3;
                }else{
                    my ($time3,$time);my $a_i=$a[$i]; my %pos2;
                    for($time=1;$time<=$para{$a_i};$time++){
                        my $time2 = ($time-1) * 7;
                        if($a[$i+$time2]){
                            my $tmp_name = &get_strand($a[$i+$time2],$a[$i+3+$time2]);
                            if($time2 == 0){
                                $pos2{$name}{$chr}{$time}=$a[$i+1+$time2];
                                $time3 = $time2;
                            }elsif($time2 > 0 and $a_i eq $tmp_name){
                                $a[$i+$time2]=$tmp_name;
                                $pos2{$name}{$chr}{$time}=$a[$i+1+$time2];
                                $time3 = $time2;
                            }else{last;}
                        }else{last;}
                    }
                    $time--;
                    my %min_dis_pa;
                    my $min_dis_p=&compare_dis(\%pos,$name,$chr,$time,\%pos2); ##present compare to former
                    %min_dis_pa=%$min_dis_p;
                    for(my $time4=1;$time4<=$time;$time4++){
                        my $time2 = ($time4-1) * 7;
                        if(exists $min_dis_pa{$time4}){
                            push @{$pos{$name}{$chr}{$min_dis_pa{$time4}}}, ($a[$i+1+$time2],$a[$i+1+$time2]+$a[$i+2+$time2]);
                            $seq{$name}{$chr."_p".$min_dis_pa{$time4}} .= $a[$i+5+$time2];
                            $a[$i+$time2]=$name."\.".$chr."_p".$min_dis_pa{$time4}; ##rename the chr
                        }else{
                            my $time_tag = $max_para{$a_i}+1;
                            push @{$pos{$name}{$chr}{$time_tag}}, ($a[$i+1+$time2],$a[$i+1+$time2]+$a[$i+2+$time2]);
                            if(!exists $seq{$name}{$chr."_p".$time_tag}){
                                for(my $j=1;$j<=$last_len;$j++){
                                    $seq{$name}{$chr."_p".$time_tag} .= "-";
                                }
                            }
                            $seq{$name}{$chr."_p".$time_tag} .= $a[$i+5+$time2];
                            $max_para{$a_i}++;
                            $a[$i+$time2]=$name."\.".$chr."_p".$time_tag; ##rename the chr
                        }
                    }
                    $i+=$time3;
                }
            }
##no para having
            else{
                if(!exists $seq{$name}{$chr}){
                    for(my $j=1;$j<=$last_len;$j++){
                        $seq{$name}{$chr} .= "-";
                    }
                }
                $seq{$name}{$chr} .= $a[$i+5];
###v2 start: add pos
                push @{$pos{$name}{$chr}{1}}, ($a[$i+1],$a[$i+1]+$a[$i+2]);
###v2 end: add pos
            }
##just to filter non aligned species at last
            $len{$name} += $a[4];
        }
    }

##get non aligned seq and add by gaps
    my $gap = length($a[7]);
    foreach my $key (sort keys %hash){
CIRL:
        foreach my $tag (keys %{$seq{$key}}){
            for(my $k=2;$k<=$#a;$k+=7){
                my ($name,$chr) = ($1,$2), if($a[$k]=~/(\w+)\.(\S+)/);
                if($key eq $name and $tag eq $chr){next CIRL;}
            }
            for(my $j=1;$j<=$gap;$j++){
                $seq{$key}{$tag} .= "-";
            }
        }
    }

    $last_len+=length($a[7]);
}
close IN;

my ($sp_num,$total_len);
foreach my $key (sort keys %hash){
	next, unless(exists $len{$key} and $len{$key}>0);
	foreach my $tag (sort keys %{$seq{$key}}){
###v2 start: add pos
        my ($chr,$num);
        if($tag=~/(\S+_[APM])_p([0-9]+)/){
            ($chr,$num)=($1,$2);
        }elsif($tag=~/(\S+_[APM])/){
            ($chr,$num)=($1,1);
        }
        my ($re_len,$re_pos,$last_len,$last_pos);
        for(my $i=0; $i<=$#{$pos{$key}{$chr}{$num}}; $i+=2){
            if($i==0){
                $re_pos=$pos{$key}{$chr}{$num}[$i];
                $last_pos=$pos{$key}{$chr}{$num}[$i+1];
                $last_len=$pos{$key}{$chr}{$num}[$i+1]-$pos{$key}{$chr}{$num}[$i];
                if($i+1==$#{$pos{$key}{$chr}{$num}}){
                    $re_pos.=",".$pos{$key}{$chr}{$num}[$i+1];
                    $re_len=$last_len;
                }
            }elsif($i+1==$#{$pos{$key}{$chr}{$num}}){
                if($last_pos==$pos{$key}{$chr}{$num}[$i]){
                    $re_pos.=",".$pos{$key}{$chr}{$num}[$i+1];
                    $last_len+=$pos{$key}{$chr}{$num}[$i+1]-$pos{$key}{$chr}{$num}[$i];
                    $re_len.=$last_len;
                }else{
                    $re_pos.=",".$last_pos.";".$pos{$key}{$chr}{$num}[$i].",".$pos{$key}{$chr}{$num}[$i+1];
                    $re_len.=$last_len.",";
                    $last_len=$pos{$key}{$chr}{$num}[$i+1]-$pos{$key}{$chr}{$num}[$i];
                    $re_len.=$last_len;
                }
            }else{
                if($last_pos==$pos{$key}{$chr}{$num}[$i]){
                    $last_pos=$pos{$key}{$chr}{$num}[$i+1];
                    $last_len+=$pos{$key}{$chr}{$num}[$i+1]-$pos{$key}{$chr}{$num}[$i];
                }else{
                    $re_pos.=",".$last_pos.";".$pos{$key}{$chr}{$num}[$i];
                    $last_pos=$pos{$key}{$chr}{$num}[$i+1];
                    $re_len.=$last_len.",";
                    $last_len=$pos{$key}{$chr}{$num}[$i+1]-$pos{$key}{$chr}{$num}[$i];
                }
            }
        }
###v2 end: add pos
		print ">$key.$tag\t$re_pos\t$re_len\n$seq{$key}{$tag}\n";
	}
}


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


######################################################################
##get the minimized pairwise distance to connect sequences
######################################################################
sub compare_dis{
    my ($pos,$name,$chr,$time,$pos2) = @_;
    my (%last,%last1,%last2,%last1D,%last2D,);
    my $dis=9999999999;

##get the pairwise distance
    for(my $i=1; $i<=$time; $i++){
        foreach my $j (sort keys %{$pos->{$name}->{$chr}}){
            if($pos->{$name}->{$chr}->{$j}[-1] > $pos2->{$name}->{$chr}->{$i}){$last1{$i}{$j}=$dis;}
            else{$last1{$i}{$j} = $pos2->{$name}->{$chr}->{$i} - $pos->{$name}->{$chr}->{$j}[-1];}
        }
    }

    foreach my $j (sort keys %{$pos->{$name}->{$chr}}){
        for(my $i=1; $i<=$time; $i++){
            if($pos->{$name}->{$chr}->{$j}[-1] > $pos2->{$name}->{$chr}->{$i}){$last2{$j}{$i}=$dis;}
            else{$last2{$j}{$i} = $pos2->{$name}->{$chr}->{$i} - $pos->{$name}->{$chr}->{$j}[-1]}
        }
    }

##compare the pairwise distance
    foreach my $i (sort keys %last1){
        my $c1=1;
        foreach my $j (sort {$last1{$i}{$a}<=>$last1{$i}{$b}} keys %{$last1{$i}}){
            next, if(exists $last2D{$j});
            if($c1==1){
                my $c2=1;my $lastc2;
                foreach my $is (sort {$last2{$j}{$a}<=>$last2{$j}{$b}} keys %{$last2{$j}}){
                    next, if(exists $last1D{$is});
                    if($c2==1){
##if the min pair is same, they are the best choice
                        if($last1{$i}{$j} == $last2{$j}{$is} && $i == $is && $last1{$i}{$j} != $dis){
                            $last{$i}=$j;
                            $last1D{$i}=1;$last2D{$j}=1;
                            $c1++;
                        }
                        $lastc2=$last2{$j}{$is};
                        $c2++;
                    }
                    else{
                        if($lastc2 == $last2{$j}{$is} && $lastc2 != $dis){
                            print STDERR "There are the (1st) same position $lastc2 aligned to reference at $pos->{$name}->{$chr}->{$j}[-2] $name/$chr/p$j!\n";
                        }
                        $c2=1;last;
                    }
                }
            }else{
                $c1=1;last;
            }
        }
    }
    return \%last;
}
