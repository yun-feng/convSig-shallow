use warnings;
use strict;

my $dir="/data/ted/COAD/";
my $out="colorectal.txt";
my $input="simple_somatic_mutation.open.COCA-CN.tsv";
open OUT,">$dir$out" or die "Can't open output $dir$out";
my @line;
my $mut;

sub Base_to_number($){
        $_=shift;
        return 0 if $_ eq "A";
        return 1 if $_ eq "G";
        return 2 if $_ eq "C";
        return 3 if $_ eq "T";
        return 4;
}


open FH, "$dir$input" or die "Can't open input $dir$input";

while(<FH>){
	chomp;
	@line=(split "\t", $_);
	print "error:",$line[12],"\n", if $line[12]	ne "GRCh37";
	next , if $line[12]	ne "GRCh37";
	next , if $line[8] eq 'X';
	next , if $line[8] eq 'Y';
	next , if Base_to_number($line[15])>3;
	next , if Base_to_number($line[16])>3;
#	print (join "\t",@line[4,8,9,10,15,16]);
#	print "\n";
#	print $line[33],"\n";
	next , if $line[33] ne 'WGS';
	print OUT (join "\t",@line[4,8,9,10,15,16]);
	print OUT "\n";
}
