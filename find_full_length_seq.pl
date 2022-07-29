#!/usr/bin/perl -w                                                                                                                                                                                                                                                                                                                        

use strict;
use warnings;

###takes protein seq fasta and removes sequences without a starting M and ending *###
#Input 1: Transdecoder protein fasta final output
open (my $fasta_in, "<", "/path/to/transdecoder/output/protein_fasta.alltranscripts.fasta.transdecoder.pep");
#Output 1: protein fasta with only full length proteins
open (my $fasta_out, ">", "/path/to/fulllength/protein/coding/output_fulllength_seqs.pep");


my $counter = 0;
my %seq_list;
my $name;
my $seq;
my $good = 0;
while (my $fasta = <$fasta_in>){
	chomp $fasta;
	if ($fasta =~ />/){
		$counter = 1;
		if ($good == 1){
			$seq_list{$name} = $seq;
			$good = 0;
		} else {
		}
		$name = $fasta;
	} else {
		if ($counter == 1){
			$seq = $fasta;
			$counter++;
			$good = 1;
		} elsif ($counter >= 1){
			$seq = $seq . $fasta;
		}
	}
}

my $num;
my $no_M = 0;
my $no_end = 0;
foreach my $name (sort keys %seq_list){
	my $seq = $seq_list{$name};
	my $first = substr($seq, 0, 1);
	my $last = substr($seq, -1);
	if ($first eq "M"){ #start codon
		if ($last eq "*"){ #stop codon
			print $fasta_out "$name\n$seq\n";
			$num++;
		} else {
			$no_end++;
		}
	} else {
		$no_M++;
	}
}

print "number of full length sequences = $num\nnumber without start = $no_M\nnumber without stop = $no_end\n";
close($fasta_in);
close($fasta_out);