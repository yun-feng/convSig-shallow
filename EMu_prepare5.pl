use warnings;
use strict;

my $fasta_file="/well/htseq/Genomes/BWA/GRCh37";
my $mut_file="/data/ted/COAD/colorectal-de.txt";

sub Base_to_number($){
	$_=shift;
	return 0 if $_ eq "A";
	return 1 if $_ eq "G";
	return 2 if $_ eq "C";
	return 3 if $_ eq "T";
	return 4;
}

sub Feature_to_table($$$$$){
	my ($pos1,$pos2,$pos3,$pos4,$pos5)=@_;
	my ($n1,$n2,$n3,$n4,$n5)=(Base_to_number($pos1),Base_to_number($pos2),Base_to_number($pos3),Base_to_number($pos4),Base_to_number($pos5));
	
	return () if $n1>3;
	return () if $n2>3;
	return () if $n3>3;
	return () if $n4>3;
	return () if $n5>3;
	
	if($n3<2){
		($n1,$n2,$n3,$n4,$n5)=(3-$n1,3-$n2,$n3,3-$n4,3-$n5);
	}
	else{
		($n1,$n2,$n3,$n4,$n5)=($n1,$n2,3-$n3,$n4,$n5);
	}
	my $value=$n1*128*3+$n2*32*3+$n3*16*3+$n4*4*3+$n5*3;
	return ($value,$value+1,$value+2);
}

my $sample=-1;
my %sample_hash;
my @sample_table;
my $sample_name;
my @chrom_table;
my @pos_table;
my @mut_type;
my @ori_type;
my ($chrom,$pos,$ori,$mut);
open MUT,$mut_file or die "Can't open mutation file $mut_file";
while(<MUT>){
        local $/="\n";
        chomp;
        ($sample_name,$chrom,$pos,$ori,$mut)=(split "\t", $_);
         $mut=Base_to_number($mut);
        next, if Base_to_number($ori)>3;
        next, if $mut>3;
        next, if ($chrom eq 'X' or  $chrom eq 'Y');
        push @chrom_table, $chrom;
        push @pos_table, $pos;
        push @ori_type, $ori;
        $mut=3-$mut, if Base_to_number($ori)<2;
        $mut=2 if $mut>2;
        push @mut_type, $mut;
        $sample_hash{$sample_name}=++$sample, if(not exists($sample_hash{$sample_name}));
        push @sample_table, $sample_name;
}


my @wt_table;
my %mut_hash;
my @line;
my ($base1,$base2,$base3,$base4,$base5);

for (keys %sample_hash){
	$sample=$_;
	$mut_hash{$sample}=[];
	push @{$mut_hash{$sample}}, 0, for(1..128*4*3);
}


push @wt_table,0, for (1..128*4*3);

$chrom=0;
$pos=-2;

open FASTA, $fasta_file or die "Can't open fasta file";

my $temp_chrom;
my $temp_pos;
my $temp_sample;
my $temp_mut;
my $temp_ori;
my @temp_array;

$temp_chrom=shift @chrom_table;
$temp_pos=shift @pos_table;

while(<FASTA>){
	chomp;
	if(/^>/){
		$chrom++;
		last if $chrom>24;
		$pos=-2;
		$base1="N";
		$base2="N";
		$base3="N";
		$base4="N";
		next;
	}
	@line=split //, $_;
	while(@line){
		$base5=shift @line;
		@temp_array=Feature_to_table($base1,$base2,$base3,$base4,$base5);
		$_++ for(@wt_table[@temp_array]);
		
		$pos++;
		
		
		if ($temp_chrom and not ($temp_chrom>$chrom || $temp_pos>$pos)){
			while(1){
					$temp_sample=shift @sample_table;
					$temp_mut=shift @mut_type;
					$temp_ori=shift @ori_type;
					print "Here: $temp_chrom	$temp_pos	$temp_ori dosn't match $base1$base2$base3$base4$base5" if $temp_ori ne $base3;
					${$mut_hash{$temp_sample}}[$temp_array[$temp_mut]]++;
					
					$temp_chrom=shift @chrom_table;
					$temp_pos=shift @pos_table;
					last if(not $temp_chrom or ($temp_chrom>$chrom || $temp_pos>$pos));
				
			}
		}
		$base1=$base2;
		$base2=$base3;
		$base3=$base4;
		$base4=$base5;
		
	}
	print $chrom," ",$pos,"\n";
	
}

my $dir="/data/ted/COAD/";
open OUT, ">${dir}}coca_5.mut.txt" or die "Can't open output file ${dir}XRXS.mut.txt";
for (keys %mut_hash){
	print OUT (join "\t", @{$mut_hash{$_}});
	print OUT "\n";
}

open OUT, ">${dir}coca_5.wt.txt" or die "Can't open output file ${dir}XRXS.wt.txt";
print OUT (join "\t", @wt_table);
print OUT "\n";


open OUT, ">${dir}coca_5.sample.txt" or die "Can't open output file ${dir}XRXS.sample.txt";
for (keys %sample_hash){
	print OUT $_,"\t",$sample_hash{$_},"\n";
}

