#!/usr/bin/env perl
# Purpose: Merging SNP calls of Identical Taxas from Hapmap file.
# Usage: MergeIdenticalTaxa.perl <input.hmp> <output.hmp>
# Contact: s.sivasubramani@cgia.org
# Date: 13-03-2019

# Checking for the Aurgument inputs
if($ARGV[1] eq ""){
	die "Usage: MergeIdenticalTaxa.perl <input.hmp> <output.hmp>\n";
}

# Hash for IUPAC to BI
%iupac2bi = (
    A => 'AA',
    T => 'TT',
    C => 'CC',
    G => 'GG',
    M => 'AC',
    R => 'AG',
    W => 'AT',
    S => 'CG',
    Y => 'CT',
    K => 'GT',
    B => 'CGT',
    D => 'AGT',
    H => 'ACT',
    V => 'ACG',
    N => 'NN'
);

# Hash for BI to IUPAC
%iupac = (
    A    => 'A',
    T    => 'T',
    C    => 'C',
    G    => 'G',
    AC   => 'M',
    AG   => 'R',
    AT   => 'W',
    CG   => 'S',
    CT   => 'Y',
    GT   => 'K',
    CGT   => 'B',
    AGT   => 'D',
    ACT   => 'H',
    ACG   => 'V',
    ATCG   => 'N',
    N    => 'N'
);

# Input Hapmap File
$inHMP = $ARGV[0];
# Output Hapmap File
$outHMP = $ARGV[1];

# Open input hapmap file with file handle FA || exit from script
open(FA, $inHMP) || die "Can not open the input Hapmap file.\n";
# open output hapmap file with file handle OUT || exit from script
open(OUT, ">$outHMP") || die "Can not create the output Hapmap file.\n";
# Reading first lind of FA
chomp($head = <FA>);
# Spliting the Hapmap header to read Taxa Names
@H = split ("\t", $head);
# Counting how many taxas are duplicated (Identical) and making a Hash of index with its incremental count
for ($i=11; $i<=$#H; $i++){
	push(@{ $taxa{$H[$i]} }, $i);
}

# Output first line concatenation
$headLine = join("\t", @H[0..10]);
foreach $key(sort keys %taxa){
	$headLine .= "\t".$key;
}
print OUT $headLine,"\n";

# Reading from second line (line-by-line)
while (<FA>) {
	chomp();
	# Number of lines processed.
	$lineNo++;
	if($lineNo % 1000 == 0){
		print "Processed $lineNo lines...\n";
	}

	# Split the hapmap line into array
	@arr = split("\t", $_); @seq = ();
	# first 12 columns of meta info into outline
	$outLine = join("\t", @arr[0..10]);
	# Reading all the taxa allele calls in the same order we have written out head. Refer line number 71
	foreach $key(sort keys %taxa){
		# If the taxa is non duplicated just add that to the output allele calles array
		if(scalar @{$taxa{$key}} == 1 ){
			$idx = $taxa{$key}[0];
			push(@seq, $arr[$idx]);
		}
		# else (for duplicated/identical taxa) of any number
		else{
			$base = "";
			# Concatenating all the allele calls of given(#84) duplicated Taxa
			foreach $idx(sort @{$taxa{$key}}){
				$base .= $iupac2bi{$arr[$idx]};
			}
			# Removing Ns (missing calls)
			$base =~ s/N//ig;
			# Removing duplicated characters
			$base =~ s/(.)(?=.*?\1)//g;
			# if the base allele call is empty, means its a missing call in all the given duplicated taxa, so considering N as allele call.
			if($base eq ""){
				$base = "N";
			}
			# For multi allele call, to read from the IUPAC has, we need to have all the permuations of calls. So we get a ARRAY of shuffled elements. Eg: GT => (GT, TG)
			@finalCall = shuffle($base);
			# Mapped allele call to the IUPAC hash to derive whats the merges allele call
			foreach $call(@finalCall){
				if($iupac{$call}){
					# Adding merges allele call to the output array
					push(@seq, $arr[$idx]);
				}
			}
		}
	}
	# printing output line to Hapmap file handle
	print OUT $outLine,"\t", join("\t", @seq),"\n";
}

# Subroutine to shuffle charectors (upto length of 3) of a string and returns an arrayt containing all possible strings.
sub shuffle (){
	$string = shift @_;
	@out = ();
	if(length($string) == 1){
		push(@out, $string);
	}
	if(length($string) == 2){
		push(@out, $string);
		$rev = reverse $string;
		push(@out, $rev );
	}
	if(length($string) == 3){
		@calls = split(//, $string);
		push(@out, "$calls[0]$calls[1]$calls[2]");
		push(@out, "$calls[0]$calls[2]$calls[1]");
		push(@out, "$calls[1]$calls[0]$calls[2]");
		push(@out, "$calls[1]$calls[2]$calls[0]");
		push(@out, "$calls[2]$calls[0]$calls[1]");
		push(@out, "$calls[2]$calls[1]$calls[0]");
	}
	return @out;
}