use warnings;
use strict;

my $dir="/data/ted/COAD/";
my $out="esophagael-de.txt";
my $input="esophagael-sorted.txt";
open OUT,">$dir$out" or die "Can't open output $dir$out";
open IN, "$dir$input" or die "Can't open input $dir$input";

my $old;
my $new;

while(<IN>){
	chomp;
	next, if $old and ($old eq $_);
	$old=$_;
	print OUT $old, "\n";
}
