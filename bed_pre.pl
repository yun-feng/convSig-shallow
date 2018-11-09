use warnings;
use strict;

my $chrom_file ="/data/ted/multi/Ubiquitous.bed";
my $out_file =">/data/ted/multi/Ubiquitous_s.bed";

my $ochr;
my $os;
my $oe;
open CHROM,$chrom_file or die "Can't open mutation file $chrom_file";
open OUT,$out_file or die "Can't open mutation file $out_file";
while(<CHROM>){
	chomp;
	($ochr,$os,$oe)=(split "\t", $_);
	$ochr =~ /^chr(.*)/;
	$ochr=$1;
	$ochr=23, if $ochr eq 'X';
	print OUT (join "\t", ($ochr,$os,$oe)),"\n";
}
