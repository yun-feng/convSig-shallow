use warnings;
use strict;

my $fasta_file="/data/ted/COAD/test-prepare/GRCh37_head";
my $mut_file="/data/ted/COAD/test-prepare/mut_head";

sub Base_to_number($){
	$_=shift;
	return 0 if $_ eq "A";
	return 1 if $_ eq "G";
	return 2 if $_ eq "C";
	return 3 if $_ eq "T";
	return 4;
}


sub Feature_to_table($$$){
	my ($pos1,$pos2,$pos3)=@_;
	my ($n1,$n2,$n3)=(Base_to_number($pos1),Base_to_number($pos2),Base_to_number($pos3));
	
	return () if $n1>3;
	return () if $n2>3;
	return () if $n3>3;
	
	if($n2<2){
		($n1,$n2,$n3)=(3-$n3,$n2,3-$n1);
	}
	else{
		($n1,$n2,$n3)=($n1,3-$n2,$n3);
	}
	my $value=$n1*8*3+$n2*4*3+$n3*3;
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
my ($base1,$base2,$base3);

for (keys %sample_hash){
	$sample=$_;
	$mut_hash{$sample}=[];
	push @{$mut_hash{$sample}}, 0, for(1..96);
}


push @wt_table,0, for (1..96);

$chrom=0;
$pos=-1;

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
	if(/^>/){   #if the next line is the start of a new chromosone
		$chrom++; # Count how many chromosomes have visited
		last if $chrom>24; #if this chromosome if between 1-22 and X,Y
		$pos=-1;  #the position points to the base before the first base of this chromosome 
		$base1="N"; #no base are there before the fist base of this chromosome
		$base2="N";
		next;
	}
	@line=split //, $_;
	while(@line){
		$base3=shift @line;  #the base that we are currently considering
		@temp_array=Feature_to_table($base1,$base2,$base3); #the numerical place for fragments of type ($base1$base2$base3)
		$_++ for(@wt_table[@temp_array]); 
		
		$pos++;
		
		
		if ($temp_chrom and not ($temp_chrom > $chrom || $temp_pos>$pos)){ #if there are still mutations from the mutation file and the position being considered correspond to the position of the mutation being considered ( this is why we need to sort the mutation file)
			while(1){ #since there are multiple samples in the mutation file, it is possible there are multiple mutations at the same position
					$temp_sample=shift @sample_table; 
					$temp_mut=shift @mut_type;
					$temp_ori=shift @ori_type;
					print "Here: $temp_chrom	$temp_pos	$temp_ori dosn't match $base1$base2$base3\n" if $temp_ori ne $base2;
					next ,if $temp_ori ne $base2; #check if the midlle base is the same as the middle base from mutation file...just to make sure we are using the same genome reference as the one used to call the mutations 
					${$mut_hash{$temp_sample}}[$temp_array[$temp_mut]]++; #add one to the mutation table
					
					$temp_chrom=shift @chrom_table; #consider the next mutation
					$temp_pos=shift @pos_table;
					last if(not $temp_chrom or ($temp_chrom > $chrom || $temp_pos>$pos)); # if the next mutation is not at this position, break 
				
			}
		}
		$base1=$base2; #slide the position by 1
		$base2=$base3;
		
	}
	print $chrom," ",$pos,"\n";
	
}

my $dir="/data/ted/COAD/test-prepare/";
open OUT, ">${dir}eso.mut.txt" or die "Can't open output file ${dir}COAD.mut.txt";
for (keys %mut_hash){
	print OUT (join "\t", @{$mut_hash{$_}});
	print OUT "\n";
}

open OUT, ">${dir}eso.wt.txt" or die "Can't open output file ${dir}COAD.wt.txt";
print OUT (join "\t", @wt_table);
print OUT "\n";


open OUT, ">${dir}eso.sample.txt" or die "Can't open output file ${dir}COAD.sample.txt";
for (keys %sample_hash){
	print OUT $_,"\t",$sample_hash{$_},"\n";
}
