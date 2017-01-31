#!/usr/bin/perl
# author: Astride Audirac
# date: 15/05/15
# script to dehomopolymerize sequences CCCCC-->CCC 
# to use this script type: perl dehomopolymerize.pl [inputfile] [outputfile] \n

use strict;
use warnings;

my $j=0;
my $k=0;
my $l=0;
my $m=0;

my $file_name = $ARGV[0];
chomp $file_name;
my $output= $ARGV[1];
open (IN,"<",$file_name) or die "error in: $!";
my @row=<IN>;
open (OUT,">","$file_name.DH") or die "error in: $!";

for (my $p=0; $p < scalar @row; $p ++)
{
	if (defined $row[$p]){
		while($row[$p]=~ m/AAAA/ or $row[$p]=~m/CCCC/ or $row[$p]=~m/TTTT/ or $row[$p]=~m/GGGG/){ 
			if ($row[$p]=~ m/AAAA/) {
			$row[$p]=~s/AAAA/AAA/;
			$j=$j+1;
			
			}
			elsif ($row[$p]=~ m/TTTT/){
			$row[$p]=~s/TTTT/TTT/;
			$k=$k+1;
			
			}
			elsif ($row[$p]=~ m/GGGG/){
			$row[$p]=~s/GGGG/GGG/;
			$l=$l+1;
			}
			elsif ($row[$p]=~ m/CCCC/){
			$row[$p]=~s/CCCC/CCC/;
			$m=$m+1;
			}
			 
		}	#end while
	print OUT $row[$p];
	

	}
}
close OUT;
close IN;
my $datestring = localtime();
open (LOG,">","$file_name log.txt") or die "unable do log";
print LOG "$file_name dehomopolymerized \n$datestring:\n";
print LOG " $j A where removed \n";
print LOG " $k T where removed \n";
print LOG " $l G where removed \n";
print LOG " $m C where removed \n";
my $r=$j+$k+$l+$m;
print LOG "bases removed $r.\n";
close LOG;
exit;

