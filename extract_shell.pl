#!/usr/bin/perl -w

=head1 Description

This script is used to create the work shell which is used to extract MAF from HAL.

=head1 Contact & Version

  Author: Qi Fang, fangqi@genomics.cn
  Version: 1.0,  Date: 2019-01-17
  Version: 1.1,  Date: 2019-02-18, 1. --unique of hal2maf will lost some columns, deleted; 2. add 03.splitLocallyColinear.pl to split some inversions on the same scaffold
  Version: 2.0,  Date: 2019-09-11 Modified: discard Method 1 and 2; add --ortholog to output the longest and most identical locally colinear alignments and --species_list
  Version: 2.1,  Date: 2019-10-21 Modified: add hal2maf detection form environment

=head1 Command-line Option

  --hal <str>       HAL alignments file
  --input <str>     input genome length BED format file
  --pos <str>       whether show the position information, default is y(/n);
  --ref <str>       reference of multiple alignments
  --split_dir <str> output directory
  --window <num>    window length, default 100k; 1 means the original BED information will be used, don't split to windows
  --ortholog        suppress the merge activities, will only output the longest and most identical locally colinear alignments
  --species_list <str>      list of species needed in results (the same as the species name in guide tree), if not, the guidetree would be used to extracted without ancestors. Notice: species name should not include any number!

  --help            output help information to screen

=head1 Usage Exmples

  perl extract_shell.pl -i Anc001.fasta.len.bed -r Anc001 --split_dir Anc001_extract --hal allspecies.hal

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin);

my ($Input,$Ref,$Split_dir,$Win,$Hal,$Help,$Ortholog,$Pos,$SpeciesList);
GetOptions(
        "ref:s"=>\$Ref,
        "input:s"=>\$Input,
        "split_dir:s"=>\$Split_dir,
        "species_list:s"=>\$SpeciesList,
        "window:n"=>\$Win,
        "hal:s"=>\$Hal,
        "ortholog!"=>\$Ortholog,
        "pos:s"=>\$Pos,
        "help"=>\$Help
);

die `pod2text $0` if ($Help);

$Win ||= 100000;
$Pos ||= "y";

if($Ortholog){
    $Ortholog="--ortholog";
}else{
    $Ortholog="";
}

my $speciesList;
if($SpeciesList){
    $speciesList="--species_list ".$SpeciesList;
}else{
    $speciesList="";
}

if($Win == 1){
    `cp $Input $Input.split`;
}else{
    `python $Bin/split_bed.py $Input $Win > $Input.split`;
}

open IN,"<$Input.split";
mkdir "$Split_dir" unless (-d "$Split_dir");
my $parentdir="00";
my $subdir = "000";
my $parentloop=0;
my $loop = 0;
my $cmd;

while(<IN>){
    if($loop % 100 == 0){
        if($parentloop % 100 ==0){
            $parentdir++;
            mkdir ("$Split_dir/$parentdir");
            $subdir="000";
        }
        $subdir++;
        mkdir("$Split_dir/$parentdir/$subdir");
        $parentloop++;
    }
    chomp;
    my @a = split(/\s+/,$_);
    my $qr_file = "$Split_dir/$parentdir/$subdir/$a[3].bed";
    open OUT, ">$qr_file" || die "fail creat $qr_file";
    print OUT "$a[0]\t$a[1]\t$a[2]\t$a[3]\n";
    close OUT;
    my $Hal2Maf = `which hal2maf`;
    chomp $Hal2Maf;
    $cmd .= "$Hal2Maf --onlyOrthologs --refGenome $Ref --refTargets $qr_file $Hal $qr_file.0.maf && perl $Bin/91.correctMafName.pl -i $qr_file.0.maf > $qr_file.maf && perl $Bin/01.sortMaf.pl -i $qr_file.maf -r $Ref > $qr_file.maf.sort && perl $Bin/02.splitDuplicates.pl -i $qr_file.maf.sort $speciesList > $qr_file.maf.sort.all.fasta && perl $Bin/03.splitLocallyColinear.pl $qr_file.maf.sort.all.fasta > $qr_file.maf.sort.all.s.fasta && perl $Bin/04.mergeColinearBlocks.pl -i $qr_file.maf.sort.all.s.fasta $Ortholog -p $Pos -r $Ref > $qr_file.maf.sort.all.fasta.merge && rm -rf $qr_file.0.maf $qr_file.maf.sort $qr_file.maf.sort.all.fasta && gzip $qr_file.maf && gzip $qr_file.maf.sort.all.s.fasta\n";
    $loop++;
}

open OUT, ">split_hal2fa_$Ref.sh" || die "fail creat split_hal2fa_$Ref.sh";
print OUT $cmd;
close OUT;
