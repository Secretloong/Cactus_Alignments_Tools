#!/usr/bin/perl -w

=head1 Description

This script is used to merge the best homology groups:
    a)  Method 1 (seg) == discarded
        i.Segment (homology group) based method;
        ii.Iterate to map the longest homology groups back to reference, only the next non-overlapped longest homology group would be merged to the first longest homology group;
        iii.When the next non-overlapped longest homology group was merged into the gaps of the former longest homology groups, the gaps region should be longer than 10bp. Otherwise this next homology group should be discarded. This will reduce chimerism of different sequences coming from the very small homology group.
    b)  Method 2 (base) == discarded
        i.Base (each bp) based method; Step is one bp, but window is _len_ (default is 10bp);
        ii.Select the base from the longest homology group of all the homology groups for each base in the multiple alignments;
        iii.To reduce chimerism of different sequences, the continuous bases from the same homology group should be longer than 10 bp and the break points (in other words, perhaps these break points are the inserts in reference, the deletion in targets or just the error in alignments) should less than 2, then them will be merged into the result. NEW, the break points have been also retained.
    c)  Method 3 (new) == default
        i.like Base-based method;
        ii.Select the base from the best homology group of all the homology groups for each base in the multiple alignments; the best one means one having the highest identity in each length ranks; all homology groups are seperated to different length ranks by 100bp.
        iii.the break points have been also retained, because the indels in ref have been considerated in identity calculation.

=head1 Contact & Version

  Author: Qi Fang, fangqi@genomics.cn
  Version: 1.0,  Date: 2018-12-20
  Version: 2.0,  Date: 2019-01-30  Modified: using the identity and length as the treshold to filter the overlaped bases.
  Version: 2.1,  Date: 2019-06-27  Debug: 1. get_pos() when the $LenV longer than the blocks, make sure the end ordinates would be correct 
                                             when the righest block only alined partially.
                                          2. initail the $last_i (=-1) when the homology groups selection went back to the first rank groups.
  Version: 3.0,  Date: 2019-09-11 Modified: discard Method 1 and 2; add --ortholog to output the longest and most identical locally colinear alignments
  Version: 3.1,  Date: 2019-09-17 Modified: add test to sequences only with N; add "n" to "N" identification (some alignments having masked n); change $/ back to \n
  Version: 3.2,  Date: 2020-05-07 Modified: debug the coordinates within the total "n" or "N" blocks, don't affect the alignments
  Version: 3.3,  Date: 2020-09-15 Modified: add test to all N blocks in unique target sequences, so that the subsequent analysis could treat total "n" or "N" blocks
                                            for all blocks uniformly, don't affect the alignments
  Version: 3.4,  Date: 2020-09-18 Modified: add $len_con_base_InterN and the N (in one block) spliting in get_pos function
  Version: 4.0,  Date: 2020-09-21 Modified: rewrite the N spliting, so that the initial/internal N could all be deal with
  Version: 4.1,  Date: 2020-10-04 Debug: line 473

=head1 Command-line Option

  --input <str>     input duplicates splited FASTA file
  --reference <str> the species name which is based on to extract these alignments, default Anc001
  --len <num>       length to reduce influence of little sequences, default 10
  --pos <str>       whether calculate the position information, "y" or "n", default is yes ("y"), it would be much slower
  --ortholog        suppress the merge activities, will only output the longest and most identical locally colinear alignments

  --help            output help information to screen

=head1 Usage Exmples

  04.mergeColinearBlocks.pl -i Anc001.maf.sort.all.s.fasta -r Anc001 > Anc001.maf.sort.all.s.fasta.merge

=cut

use strict;
use Getopt::Long;

my ($Input,$LenV,$PosInfo,$Ref,$Ortholog,$Help);
GetOptions(
        "input:s"=>\$Input,
        "len:n"=>\$LenV,
        "pos:s"=>\$PosInfo,
        "reference:s"=>\$Ref,
        "ortholog!"=>\$Ortholog,
        "help"=>\$Help
);

die `pod2text $0` if ($Help);

$Ref ||= "Anc001";
$LenV ||= 10;
$PosInfo ||= "y";

my (%segment,%seq,%len,%out,$Ref_seq,$Ref_segment,$Ref_len);
open IN,$Input;
$/ = ">";
<IN>;
while(<IN>){
    chomp;
    my @a = split /\s+/;
    my ($sp,$seg)=($1,$2), if($a[0]=~/(\w+)\.(\S+)/);
    $seq{$sp}{$seg}=$a[3];

    my @b = split(/,/,$a[2]);
    my $tmp_len;
    for(my $i=0;$i<=$#b;$i++){$tmp_len+=$b[$i];}
    $len{$sp}{$seg}=$tmp_len;

    my @c = split(/;/,$a[1]);
    for(my $i=0;$i<=$#c;$i++){$segment{$sp}{$seg}.=$seg.":".$c[$i].";";}
    chop($segment{$sp}{$seg});

    if($sp eq $Ref){$Ref_seq=$seq{$sp}{$seg};$Ref_segment=$segment{$sp}{$seg};$Ref_len=$len{$sp}{$seg};}

}close IN;
$/ = "\n";

##################################################
#### Merge duplicates through window, 
#### and consider identity to reference
##################################################

##print Reference
if($PosInfo eq "y"){
    $out{$Ref}=$Ref_segment."\n".$Ref_seq;
}else{
    $out{$Ref}=$Ref_seq;
}

##merge for each species
foreach my $sp (sort keys %segment){
    next, if($sp eq $Ref);

###sort the segments by length
    my @len_seg = sort {$len{$sp}{$b}<=>$len{$sp}{$a}} keys %{$segment{$sp}};

    my (%match,%mismatch,%gap1,%gap2);
###Parsing the segment positions
    my (%pos_site,%pos_name);
    for(my $i=0;$i<=$#len_seg;$i++){
        my @a = split(/;/,$segment{$sp}{$len_seg[$i]});
        foreach my $a (@a){
            my @b = split(/:/,$a);
            push @{$pos_name{$i}}, $b[0];
            my @c = split(/,/,$b[1]);
####calculate break point for target
            if(exists $pos_site{$i}){
                if($c[0]-$pos_site{$i}[-1]>0){
                    $gap1{$i}++;
                }
            }
            push @{$pos_site{$i}}, ($c[0],$c[1]);
        }
    }

    if($#len_seg==0){
###add test to all N sequences
        my $tmp_out_test=$seq{$sp}{$len_seg[0]};
        $tmp_out_test=~s/[-Nn]//g;
        if($tmp_out_test=~/^$/){next;}

        if($PosInfo eq "y"){
###add test to all N blocks, 20200915
            my $tmp_t_seg;
            my ($tj,$tk,$tj3,$tag_tn)=(0,0,0,0);
            for(my $tn=0;$tn<length($seq{$sp}{$len_seg[0]});$tn++){
                my $ts = substr($seq{$sp}{$len_seg[0]},$tn,1);
                if($ts eq "-" or $ts eq "N" or $ts eq "n"){
                    if($tk>0){
                        for(my $j=$tj3; $j<$tj; $j+=2){
                            my $tjn=$j/2;
                            $tmp_t_seg.=$pos_name{0}[$tjn].":".$pos_site{0}[$j].",".$pos_site{0}[$j+1].";";
                        }
                        my $tjn=$tj/2;
                        my $new_end=$pos_site{0}[$tj]+$tk;
                        $tmp_t_seg.=$pos_name{0}[$tjn].":".$pos_site{0}[$tj].",".$new_end.";";
                        if($pos_site{0}[$tj]+$tk == $pos_site{0}[$tj+1]){
                            $tj+=2;
                        }else{
                            $pos_site{0}[$tj]=$pos_site{0}[$tj]+$tk;
                        }
                        $tj3=$tj;$tk=0;
                    }
                    if($ts eq "N" or $ts eq "n"){
                        if($pos_site{0}[$tj+1]-$pos_site{0}[$tj] < $tk+1 and $tj<$#{$pos_site{0}}-1 ){$tk=0;$tj+=2;$pos_site{0}[$tj]++;}
                        elsif($pos_site{0}[$tj+1]-$pos_site{0}[$tj] == $tk+1 and $tj<$#{$pos_site{0}}-1){$tk=0;$pos_site{0}[$tj]++;$tj+=2;}
                        else{$pos_site{0}[$tj]++;} ##Identity contains N, but here is not containing N
                        $tag_tn=1;
                    }else{$tag_tn=0;}
                }else{
                    $tk++;
                    if($pos_site{0}[$tj+1]-$pos_site{0}[$tj] < $tk and $tj<$#{$pos_site{0}}-1){$tk=1;$tj+=2;}
                }
            }
            if($tk>0 or $tag_tn==1){ ##$tag_tn==1 make sure the last N blocks could be emitted to seg, than deleted by merge_pos
                for(my $j=$tj3; $j<=$tj; $j+=2){
                    my $jn=$j/2;
                    $tmp_t_seg.=$pos_name{0}[$jn].":".$pos_site{0}[$j].",".$pos_site{0}[$j+1].";";
                }
                $tj3=$tj;
                $tag_tn=0;
            }
            my $tmp_t_seg2=&merge_pos($tmp_t_seg);
            $out{$sp}=$tmp_t_seg2."\n".$seq{$sp}{$len_seg[0]};
        }else{
            $out{$sp}=$seq{$sp}{$len_seg[0]};
        }
        next;
    }

###identity calculating
    my %gt2;
    for(my $n=0;$n<$Ref_len;$n++){
        my $sr = uc(substr($Ref_seq,$n,1));
        for(my $i=0;$i<=$#len_seg;$i++){
            if(!exists $gt2{$i}){$gt2{$i}=0;}
            my $st = uc(substr($seq{$sp}{$len_seg[$i]},$n,1));
            if($sr eq $st and $sr ne "-"){
                $match{$i}++;
                $gt2{$i} = 0;
            }elsif($sr ne $st and $sr ne "-" and $st ne "-"){
                $mismatch{$i}++;
                $gt2{$i} = 0;
            }elsif($st eq "-" and $sr ne "-"){
                if($gt2{$i} == 0){$gap2{$i}++;$gt2{$i}=1;}
                else{$gt2{$i}=1;}
            }
        }
    }
###initializing with the first segment (longest one or the highest identity one) and reorder the hash of all homology groups
    my (%ident,%identG);
    my ($tmp_out,$tmp_out_seg);
    my $chosen_i = 0;
    my @new_order;
    my ($newi,$tmp_i)=(0,0);my %tmp_order;
    for(my $i=0;$i<=$#len_seg;$i++){
        if(!exists $gap1{$i}){$gap1{$i}=0;}
        if(!exists $gap2{$i}){$gap2{$i}=0;}
        if(!exists $match{$i}){$match{$i}=0;}
        if(!exists $mismatch{$i}){$mismatch{$i}=0;}
        $identG{$i}=$match{$i}/($match{$i}+$mismatch{$i}+$gap1{$i}+$gap2{$i});
        
        $tmp_order{$i}=$identG{$i};
        if($i+1<=$#len_seg and ($len{$sp}{$len_seg[$tmp_i]}-$len{$sp}{$len_seg[$i+1]}) > 100 and ($len{$sp}{$len_seg[$tmp_i]}-$len{$sp}{$len_seg[$i+1]}) > (length($seq{$sp}{$len_seg[0]})/100)){
            push @new_order, (sort {$tmp_order{$b}<=>$tmp_order{$a}} keys %tmp_order);
            %tmp_order=();$tmp_i=$i+1;
        }
    }
    push @new_order, (sort {$tmp_order{$b}<=>$tmp_order{$a}} keys %tmp_order);
    $chosen_i = $new_order[0];

    $tmp_out = $seq{$sp}{$len_seg[$chosen_i]};

    if($Ortholog){
        $out{$sp} = $segment{$sp}{$len_seg[$chosen_i]}."\n".$tmp_out;
        next;
    }

###$j1: $pos; $k1: extending length; $j3: last $j1 (next initial);$tag_j: last base was aligned;$tag_n: last base was N;
    my ($j1,$k1,$j3,$tag_j,$tag_n)=(0,0,0,0,0);
    my ($len_con_base,$len_con_base_InterGap,$len_con_base_InterN,$last_i)=(0,0,0,-1);

    for(my $n=0;$n<length($seq{$sp}{$len_seg[$chosen_i]});$n++){
        my $s1 = substr($tmp_out,$n,1);
        if($s1 eq "-" or $s1 eq "N" or $s1 eq "n"){
###non-aligned region
####last base was aligned, record that and initialize the next cycle
            if($PosInfo eq "y"){
                if($tag_j==1){
                    for(my $j=$j3; $j<$j1; $j+=2){
                        my $jn=$j/2;
                        $tmp_out_seg.=$pos_name{$chosen_i}[$jn].":".$pos_site{$chosen_i}[$j].",".$pos_site{$chosen_i}[$j+1].";";
                    }
                    my $jn=$j1/2;
                    my $new_end=$pos_site{$chosen_i}[$j1]+$k1;
                    $tmp_out_seg.=$pos_name{$chosen_i}[$jn].":".$pos_site{$chosen_i}[$j1].",".$new_end.";";
                    if($pos_site{$chosen_i}[$j1]+$k1 == $pos_site{$chosen_i}[$j1+1]){
                        $j1+=2;
                    }else{
                        $pos_site{$chosen_i}[$j1]=$pos_site{$chosen_i}[$j1]+$k1;
                    }
                    $j3=$j1;$k1=0;
                }
                if($s1 eq "N" or $s1 eq "n"){
                    if($pos_site{$chosen_i}[$j1+1]-$pos_site{$chosen_i}[$j1] < $k1+1 and $j1<$#{$pos_site{$chosen_i}}-1 ){$k1=0;$j1+=2;$pos_site{$chosen_i}[$j1]++;}
#debug the only N blocks in coordinates, fangqi 20200507 0:50, `$pos_site{$chosen_i}[$j1]++` before `$j1+=2;`
                    elsif($pos_site{$chosen_i}[$j1+1]-$pos_site{$chosen_i}[$j1] == $k1+1 and $j1<$#{$pos_site{$chosen_i}}-1){$k1=0;$pos_site{$chosen_i}[$j1]++;$j1+=2;}
                    else{$pos_site{$chosen_i}[$j1]++;} ##Identity contains N, but here is not containing N 
                    $tag_n=1;
                }else{$tag_n=0;}
                $tag_j=0;
            }

####select the best segment
            for(my $i=1;$i<=$#new_order;$i++){
                my $s2 = substr($seq{$sp}{$len_seg[$new_order[$i]]},$n,1);
                if($s2 eq "-" or $s2 eq "N" or $s2 eq "n"){
####add judgment to the second level segments about the inter gaps(-/N)
                    if($last_i == $i){
                        $len_con_base_InterGap++;
                        if($s2 eq "N" or $s2 eq "n"){
                            $len_con_base_InterN++;
                        }
                    }
                    next;
                }else{
                    if($last_i == $i){
                        $len_con_base++;
                    }else{
                        $len_con_base=1;
                        $len_con_base_InterGap=0;
                        $len_con_base_InterN=0;
                    }
                    $last_i=$i;
                    if($len_con_base>$LenV){
                        if($len_con_base==$LenV+1){
                            if($PosInfo eq "y"){
                                my $pos_tmp = &get_pos($seq{$sp}{$len_seg[$new_order[$i]]},$segment{$sp}{$len_seg[$new_order[$i]]},$n,$LenV,$len_con_base_InterN);
                                my @pos_all = @$pos_tmp;
                                for(my $p=0;$p<=$#pos_all-2;$p+=3){
                                    $tmp_out_seg.=$pos_all[$p].":".$pos_all[$p+1].",".$pos_all[$p+2].";";
                                }
                            }
                            $s2 = substr($seq{$sp}{$len_seg[$new_order[$i]]},$n-$len_con_base-$len_con_base_InterGap+1, $len_con_base+$len_con_base_InterGap);
                            substr($tmp_out, $n-$len_con_base-$len_con_base_InterGap+1, $len_con_base+$len_con_base_InterGap, $s2);
                        }else{
                            if($PosInfo eq "y"){
                                my $pos_tmp = &get_pos($seq{$sp}{$len_seg[$new_order[$i]]},$segment{$sp}{$len_seg[$new_order[$i]]},$n,0,0);
                                my @pos_all = @$pos_tmp;
                                for(my $p=0;$p<=$#pos_all-2;$p+=3){$tmp_out_seg.=$pos_all[$p].":".$pos_all[$p+1].",".$pos_all[$p+2].";";}
                            }
                            substr($tmp_out, $n, 1, $s2);
                        }
                        last;
                    }else{last;}
                }
            }
        }else{
###aligned region
            if($PosInfo eq "y"){
                $tag_j=1;
                $tag_n=0;
                $k1++;
                $last_i=-1;
                if($pos_site{$chosen_i}[$j1+1]-$pos_site{$chosen_i}[$j1] < $k1 and $j1<$#{$pos_site{$chosen_i}}-1){$k1=1;$j1+=2;}
            }
        }
    }

###add test to all N sequences
    my $tmp_out_test=$tmp_out;
    $tmp_out_test=~s/[-Nn]//g;
    if($tmp_out_test=~/^$/){next;}
###
    if($PosInfo eq "y"){
        if($tag_j==1 or $tag_n==1){
            for(my $j=$j3; $j<=$j1; $j+=2){
                my $jn=$j/2;
                $tmp_out_seg.=$pos_name{$chosen_i}[$jn].":".$pos_site{$chosen_i}[$j].",".$pos_site{$chosen_i}[$j+1].";";
            }
            $j3=$j1;
            $tag_n=0;
        }
        my $tmp_out_seg2=&merge_pos($tmp_out_seg);
        if($tmp_out_seg2=~/;$/){chop($tmp_out_seg2);}
        $out{$sp}=$tmp_out_seg2."\n".$tmp_out;
    }else{
        $out{$sp}=$tmp_out;
    }

}

foreach my $sp (sort keys %out){
    if($PosInfo eq "y"){
        print ">$sp\t$out{$sp}\n";
    }else{
        print ">$sp\n$out{$sp}\n";
    }
}


##################################################
#### subfunction for per base method
##################################################
sub get_pos{
    my ($seq,$seg,$n,$len,$intN)=@_;
    my $tmpN = $intN; ##used to present the dynamic stats of N in several blocks
    my ($tmpG,$tmpG0) = (0,0);
    my (@pos_site,@pos_name);
    my @a = split(/;/,$seg);
    my $len0=$len;
    foreach my $a (@a){
        my @b = split(/:/,$a);
        push @pos_name, $b[0];
        my @c = split(/,/,$b[1]);
        push @pos_site, ($c[0],$c[1]);
    }

    my @tmp_pos;
    my ($j1,$k1,$k2)=(0,0,0);
    for(my $i=0; $i<=$n; $i++){
        my $s = substr($seq,$i,1);
        if($s eq "-"){
            next;
        }else{
            $k1++;
            if($s ne "N" and $s ne "n"){$k2++;}
            if($pos_site[$j1+1]-$pos_site[$j1] < $k1 and $j1<$#pos_site-1){
                $k1=1;
                if($s ne "N" and $s ne "n"){$k2=1;}else{$k2=0;}
                $j1+=2;
            }
        }
    }
    if($k2>=$len+1){ ##locate in one block
        $len=$len+1;
        my ($tmp_no_N_len,$tmp_have_N_len,$firN,$intialN);
        ($tmp_no_N_len,$tmp_have_N_len,$firN,$tmpN,$tmpG,$intialN)= &test_len($seq,$n,$k1,$len,0,$tmpN);
        $tmpG0+=$tmpG;
        $tmp_have_N_len=$tmp_have_N_len-$intialN;
        my $tmp_start = $pos_site[$j1]+$k1-$len-$firN;
        if($firN>0){
            my $tmp_pos_s = &block_split(\@pos_name,\@pos_site,$seq,$tmp_start,$n,$j1,$len,$tmp_have_N_len,$firN); ##non-N total length,including-Ntotal length,N-length)
            my @tmp_pos_a = @$tmp_pos_s;
            push @tmp_pos, @tmp_pos_a;
        }else{
            push @tmp_pos, ($pos_name[$j1/2],$tmp_start,$pos_site[$j1]+$k1);
        }
    }else{ ##locate in several blocks
        for(my $j2=$j1;$j2>=0;$j2=$j2-2){
            if($j2==$j1){ ##first block
                $len=$len+1-$k2;
                my ($tmp_no_N_len,$tmp_have_N_len,$firN,$intialN);
                ($tmp_no_N_len,$tmp_have_N_len,$firN,$tmpN,$tmpG,$intialN)= &test_len($seq,$n,$k1,$k2,0,$tmpN);
                $tmpG0+=$tmpG;
                my $tmp_start = $pos_site[$j2]+$k1-$k2-$firN;
                $tmp_have_N_len=$tmp_have_N_len-$intialN;
                if($k1==$pos_site[$j2+1]-$pos_site[$j2]){
                    if($firN>0){
                        my $tmp_pos_s = &block_split(\@pos_name,\@pos_site,$seq,$tmp_start,$n,$j2,$k2,$tmp_have_N_len,$firN);
                        my @tmp_pos_a = @$tmp_pos_s;
                        unshift @tmp_pos, @tmp_pos_a;
                    }else{
                        unshift @tmp_pos, ($pos_name[$j2/2],$tmp_start,$pos_site[$j2+1]);
                    }
                }else{
                    if($firN>0){
                        my $tmp_pos_s = &block_split(\@pos_name,\@pos_site,$seq,$tmp_start,$n,$j2,$k2,$tmp_have_N_len,$firN);
                        my @tmp_pos_a = @$tmp_pos_s;
                        unshift @tmp_pos, @tmp_pos_a;
                    }else{
                        unshift @tmp_pos, ($pos_name[$j2/2],$tmp_start,$tmp_start+$k2+$firN);
                    }
                }
            }else{
                if($len<=($pos_site[$j2+1]-$pos_site[$j2])){ ## end block or multiple-N blocks
                    my $tmp_n = $n-$len0-1-$intN-$tmpG0+$len+$tmpN; ##$len0 need to add 1, so minus 1 here
                    $tmpG=&test_intelGap($seq,$tmp_n);
                    $tmpG0+=$tmpG;
                    $tmp_n = $n-$len0-1-$intN-$tmpG0+$len+$tmpN;
                    my $tmp_len = $pos_site[$j2+1]-$pos_site[$j2];
                    my ($tmp_no_N_len,$tmp_have_N_len,$firN,$intialN);
                    ($tmp_no_N_len,$tmp_have_N_len,$firN,$tmpN,$tmpG,$intialN)= &test_len($seq,$tmp_n,0,$len,$tmp_len,$tmpN);
                    $tmpG0+=$tmpG;
                    $tmp_have_N_len=$tmp_have_N_len-$intialN;
                    if($firN>0 and $tmp_len==$firN){
                        next;
                    }elsif($firN>0){
                        my $tmp_start = $pos_site[$j2+1]-$tmp_no_N_len-$firN;
                        $tmp_len = $firN+$tmp_no_N_len;
                        my $tmp_pos_s = &block_split(\@pos_name,\@pos_site,$seq,$tmp_start,$tmp_n,$j2,$tmp_no_N_len,$tmp_have_N_len,$firN);
                        my @tmp_pos_a = @$tmp_pos_s;
                        unshift @tmp_pos, @tmp_pos_a;
                    }else{
                        unshift @tmp_pos, ($pos_name[$j2/2],$pos_site[$j2+1]-$tmp_no_N_len,$pos_site[$j2+1]);
                    }
                    $len=$len-$tmp_no_N_len;
                    if($len==0){last;}elsif($len<0){print STDERR "Wrong with \$len in $pos_name[$j2/2]:$n:$Input\n";}
                }else{ ## intel blocks
                    my $tmp_len = $pos_site[$j2+1]-$pos_site[$j2];
                    my $tmp_n = $n-($len0+1+$intN+$tmpG0-$len-$tmpN);
                    $tmpG=&test_intelGap($seq,$tmp_n);
                    $tmpG0+=$tmpG;
                    $tmp_n = $n-($len0+1+$intN+$tmpG0-$len-$tmpN);
                    my ($tmp_no_N_len,$tmp_have_N_len,$firN,$intialN);
                    ($tmp_no_N_len,$tmp_have_N_len,$firN,$tmpN,$tmpG,$intialN)= &test_len($seq,$tmp_n,$tmp_len,0,$tmp_len,$tmpN);
                    $tmpG0+=$tmpG;
                    my $tmp_start = $pos_site[$j2]+$intialN;
                    if($firN>0){
                        my $tmp_pos_s = &block_split(\@pos_name,\@pos_site,$seq,$tmp_start,$tmp_n,$j2,$tmp_no_N_len,$tmp_len,$firN);
                        my @tmp_pos_a = @$tmp_pos_s;
                        unshift @tmp_pos, @tmp_pos_a;
                    }else{
                        unshift @tmp_pos, ($pos_name[$j2/2],$tmp_start,$pos_site[$j2+1]);
                    }
                    $len=$len-$tmp_no_N_len; #debug $len=$len-($pos_site[$j2+1]-$pos_site[$j2]); at 20201004
                }
            }
        }
    }
    return(\@tmp_pos);
}

sub test_intelGap{
    my ($seq,$n)=@_;
    my $gaps=0;
    for(my $i=$n; $i>=0; $i--){
        my $s = substr($seq,$i,1);
        if($s eq '-'){
            $gaps++;
        }else{
            last;
        }
    }
    return($gaps);
}

sub test_len{
    my ($seq,$n,$block_len,$tmp_no_N_len,$tmp_have_N_len,$tmpN)=@_;
    my $tag;
    my ($intialN_tag,$intialN)=(0,0);
    my ($firN,$gaps)=(0,0);
    if($tmp_no_N_len==0){$tag=0;}else{$tag=1;}
    if($block_len==0){
        my ($tmp_no_N_len0,$tmp_have_N_len0) = ($tmp_no_N_len,$tmp_have_N_len);
        ($tmp_no_N_len,$tmp_have_N_len)=(0,0);
        for(my $i=$n; $i>=0; $i--){
            my $s = substr($seq,$i,1);
            if($s ne '-'){
                $tmp_have_N_len++; 
                if($s ne 'N' and $s ne 'n'){
                    $tmp_no_N_len++;
                }elsif($s eq 'N' or $s eq 'n'){
                    $firN++;$tmpN--;
                }
            }else{$gaps++;}
            if($tmp_no_N_len0==$tmp_no_N_len or $tmp_have_N_len0==$tmp_have_N_len){last;}
        }
    }else{
        my $hill=0;
        for(my $i=$n; $i>=0; $i--){
            my $s = substr($seq,$i,1);
            if($intialN_tag>0){
                if($s eq 'n' or $s eq 'N'){$tmpN--;$tmp_have_N_len++;$intialN++;}
                elsif($s ne '-'){last;}
                else{$gaps++;}
            }else{
                if($s ne '-'){
                    $tmp_have_N_len++, if($tag==1);
                    $hill++, if($tag==0);
                    if($s ne 'N' and $s ne 'n'){
                        $tmp_no_N_len++, if($tag==0);
                        $hill++, if($tag==1);
                    }elsif($s eq 'N' or $s eq 'n'){
                        $firN++;$tmpN--;
                    }
                }else{$gaps++;}
            }
            if($tmp_have_N_len!=$block_len and (($tag==1 and $hill==$tmp_no_N_len) or ($tag==0 and $hill==$tmp_have_N_len))){$intialN_tag++;}
            if($tmp_have_N_len==$block_len and (($tag==1 and $hill==$tmp_no_N_len) or ($tag==0 and $hill==$tmp_have_N_len))){last;}
        }
    }
    return($tmp_no_N_len,$tmp_have_N_len,$firN,$tmpN,$gaps,$intialN);
}

sub block_split{
    my ($pos_name,$pos_site,$seq,$tmp_start,$n,$j1,$len,$lenWN,$intN)=@_;
    my @pos_name = @$pos_name;
    my @pos_site = @$pos_site;
    my @tmp_pos;
    my ($tmp_start0,$k2)=($tmp_start,0);

    ##get the no gaps start
    my ($effect_len,$aln_len,$tag_N)=(0,0,0);
    my $tmpN = $intN;
    for(my $i=$n; $i>=0; $i--){
        $aln_len++;
        my $s = substr($seq,$i,1);
        if($s eq "N" or $s eq "n"){
            if($tmpN <= 0){print STDERR "Wrong with \$tmpN in $pos_name[$j1/2]:$tmpN:$n:$i:$tmp_start:$Input\n";last;}
            else{
                $tag_N++;
                $tmpN--;
                $effect_len++;
            }
        }elsif($s ne "-"){
            $effect_len++;
        }

        if($effect_len==$lenWN){last;}
    }

    if($tag_N==0){
        push @tmp_pos, ($pos_name[$j1/2],$tmp_start0,$tmp_start0+$len);
    }else{
        my $n2=0; ##count the sequential Ns
        for(my $i=$n-$aln_len+1; $i<=$n; $i++){ ##$i=$n-$len-1-$intN is wrong, because the count of $len_con_base is based on 1, the id of them is based on 0
            my $s = substr($seq,$i,1);
            if($s eq "N" or $s eq "n"){
                push @tmp_pos, ($pos_name[$j1/2],$tmp_start0,$tmp_start0+$k2);
                if($k2==0){$tmp_start0++;}
                else{$tmp_start0=$tmp_start0+$k2+1;}
                $k2=0;
            }elsif($s eq "-"){
                next;
            }else{
                $k2++;
            }
        }
        push @tmp_pos, ($pos_name[$j1/2],$tmp_start0,$tmp_start0+$k2);
    }
    return(\@tmp_pos);
}

sub merge_pos{
    my ($seg)=@_;
    my (@pos_site,@pos_name,$seg_new);
    my @a = split(/;/,$seg);
    foreach my $a (@a){
        my @b = split(/:/,$a);
        push @pos_name, $b[0];
        my @c = split(/,/,$b[1]);
        push @pos_site, ($c[0],$c[1]);
    }

    if($#pos_site>1){
        my $i;my $jn;
        for($i=0;$i<=$#pos_site-2;$i+=2){
            $jn = $i/2;
            if($pos_name[$jn] eq $pos_name[$jn+1]){
                if($pos_site[$i+1] == $pos_site[$i+2]){
                    $pos_site[$i+2]=$pos_site[$i];
                }elsif($pos_site[$i+1] == $pos_site[$i]){ ##debug the only N blocks in coordinates, fangqi 20200507 0:50
                    next;
                }else{
                    $seg_new.=$pos_name[$jn].":".$pos_site[$i].",".$pos_site[$i+1].";";
                }
            }else{
                $seg_new.=$pos_name[$jn].":".$pos_site[$i].",".$pos_site[$i+1].";";
            }
        }
        $jn = $i/2;
        unless($pos_site[$i+1] == $pos_site[$i]){ ##debug the only N blocks in coordinates, fangqi 20200507 14:40
            $seg_new.=$pos_name[$jn].":".$pos_site[$i].",".$pos_site[$i+1].";";
        }
    }else{
        $seg_new=$pos_name[0].":".$pos_site[0].",".$pos_site[1].";";
    }

    return($seg_new);
}

##################################################
#### subfunction for segmental method
#### find the overlaps, excluding '-' and 'N'
##################################################
sub find_overlap{
    my ($seg1,$seg2,$pos1,$pos2,$LenV) = @_;
    my (@tag,$pos); my ($start,$len_con_gap)=(-1,0);

###Parsing the segment positions
    my (@pos1_site,@pos1_name,@pos2_site,@pos2_name);
    my @a = split(/;/,$pos1);
    foreach my $a (@a){
        my @b = split(/:/,$a);
        push @pos1_name, $b[0];
        my @c = split(/,/,$b[1]);
        push @pos1_site, ($c[0],$c[1]);
    }
    my @a2 = split(/;/,$pos2);
    foreach my $a (@a2){
        my @b = split(/:/,$a);
        push @pos2_name, $b[0];
        my @c = split(/,/,$b[1]);
        push @pos2_site, ($c[0],$c[1]);
    }

###main
    my ($j1,$k1,$j2,$k2,$j3,$j4)=(0,0,0,0,0,0); #for segment position
    for(my $i=0;$i<length($seg1);$i++){
        my $s1 = substr($seg1,$i,1);
        my $s2 = substr($seg2,$i,1);

###discard overlapped segments
        if(($s1 ne "-" and $s1 ne "N" and $s1 ne "n") and ($s2 ne "-" and $s2 ne "N" and $s2 ne "n")){
            @tag=();$pos="OVERLAPPING!";
            last;
        }

###complement the gaps in the first sequence
        if($s1 eq "-"){
            $len_con_gap++;
            if($s2 ne "-"){
                $k2++;
                if($pos2_site[$j2+1]-$pos2_site[$j2] < $k2 and $j2<$#pos2_site-1){$k2=1;$j2+=2;}
                if($s2 ne "N" and $s2 ne "n" and $start == -1){$start = $i;}
            }
        }elsif($s1 eq "N" or $s1 eq "n"){
            $len_con_gap++;
            if($s2 ne "-"){
                $k2++;
                if($pos2_site[$j2+1]-$pos2_site[$j2] < $k2 and $j2<$#pos2_site-1){$k2=1;$j2+=2;}
                if($s2 ne "N" and $s2 ne "n" and $start == -1){$start = $i;}
            }
            if($pos1_site[$j1+1]-$pos1_site[$j1] < $k1+1 and $j1<$#pos1_site-1 ){$k1=0;$j1+=2;$pos1_site[$j1]++;}
            else{$pos1_site[$j1]++;}
        }else{
            if($len_con_gap>$LenV){
###zero-based to 1-based
                push @tag, ($start, $i);
###get the position of seq1
                if($pos1_site[$j1+1]-$pos1_site[$j1]==$k1){
                    for(my $j=$j3; $j<=$j1; $j+=2){
                        my $jn=$j/2;
                        $pos.=$pos1_name[$jn].":".$pos1_site[$j].",".$pos1_site[$j+1].";";
                    }
                    $j3=$j1+2;
                }else{
                    for(my $j=$j3; $j<$j1; $j+=2){
                        my $jn=$j/2;
                        $pos.=$pos1_name[$jn].":".$pos1_site[$j].",".$pos1_site[$j+1].";";
                    }
                    my $jn=$j1/2;
                    my $new_end=$pos1_site[$j1]+$k1;
                    $pos.=$pos1_name[$jn].":".$pos1_site[$j1].",".$new_end.";";
                    $pos1_site[$j1]=$pos1_site[$j1]+$k1;
                    $j3=$j1;
                }
###get the position of seq2
                if($pos2_site[$j2+1]-$pos2_site[$j2]==$k2){
                    for(my $j=$j4; $j<=$j2; $j+=2){
                        my $jn=$j/2;
                        $pos.=$pos2_name[$jn].":".$pos2_site[$j].",".$pos2_site[$j+1].";";
                    }
                    $j4=$j2+2;
                }else{
                    for(my $j=$j4; $j<$j2; $j+=2){
                        my $jn=$j/2;
                        $pos.=$pos2_name[$jn].":".$pos2_site[$j].",".$pos2_site[$j+1].";";
                    }
                    my $jn=$j2/2;
                    my $new_end=$pos2_site[$j2]+$k2;
                    $pos.=$pos2_name[$jn].":".$pos2_site[$j2].",".$new_end.";";
                    $pos2_site[$j2]=$pos2_site[$j2]+$k2;
                    $j4=$j2;
                }
            }else{
                $j4=$j2;
            }
            $len_con_gap=0;$start=-1;
            $k1++;
            if($pos1_site[$j1+1]-$pos1_site[$j1] < $k1 and $j1<$#pos1_site-1){$k1=1;$j1+=2;}
        }
    }

    if($j1<$#pos1_site){
        for(my $j=$j3; $j<=$j1; $j+=2){
            my $jn=$j/2;
            $pos.=$pos1_name[$jn].":".$pos1_site[$j].",".$pos1_site[$j+1].";";
        }
    }

    return (\@tag,$pos);
}
