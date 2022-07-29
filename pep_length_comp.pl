#!/usr/bin/perl -w                                                                                                                                                                                                                                                                                                                        

use strict;
use warnings;

open (my $trans_list, "<", "/Users/bhazzard/Documents/PacBioPaperData/Bloodstage_fulllength_seqs.pep");
open (my $ref_list, "<", "/Users/bhazzard//Downloads/PlasmoDB-54_PvivaxP01_AnnotatedProteins.fasta");
open (my $len_out, ">", "/Users/bhazzard/Documents/PacBioPaperData/Bloodstage_fulllength_pep_lens.txt");

print $len_out "Peptide_ID\tSource\tLength\n";
while (my $list_line = <$trans_list>){
	chomp $list_line;
	if ($list_line =~ ">"){
		my @fields = split / /, $list_line;
		my $name = substr($fields[0], 1);
		my $len = substr($fields[5], 4);
		print $len_out "$name\tPacBio\t$len\n";
	}
}
close($trans_list);
while (my $line = <$ref_list>){
	chomp $line;
	if ($line =~ ">"){
		my @fields = split /\|/, $line;
		my $name = substr($fields[0], 1);
		my $len = substr($fields[7], 16);
		print $len_out "$name\tRef\t$len\n";
	}
}
close($ref_list);
close($len_out);