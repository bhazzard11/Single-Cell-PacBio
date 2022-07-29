#!/usr/bin/perl -w                                                                                                                                                                                                                                                                                                                        

use strict;
use warnings;

open (my $gene_list, "<", "/Users/bhazzard/Documents/Blast_v54/Bloodstage_vs_ref.tab");
open (my $pep_seqs, "<", "/Users/bhazzard/Documents/PacBioPaperData/Bloodstage_fulllength_seqs.pep");
open (my $ninty_seqs, ">", "/Users/bhazzard/Documents/Blast_v54/Bloodstage_fifty_ninty_fulllength_pep_seqs.fasta");
open (my $outnum, ">", "/Users/bhazzard/Documents/Blast_v54/Bloodstage_per.txt");

###get all transcripts that are less than 90% cov and greater than 30###
my %blast_geneid;
my %over_ninty;
my %ref_match;
while (my $line = <$gene_list>){
	chomp $line;
	my @fields = split /\t/, $line;
	my $trans_id = $fields[0];
	my $cov = $fields[2];
	my $match = $fields[1];
	if ($cov > 90){
		if (defined $over_ninty{$trans_id}){
			if ($over_ninty{$trans_id} > $cov){
				$over_ninty{$trans_id} = $cov;
				$ref_match{$trans_id} = $match;
			} else {
				$ref_match{$trans_id} = $match;
				$over_ninty{$trans_id} = $cov;
			}
		} else {
			$ref_match{$trans_id} = $match;
			$over_ninty{$trans_id} = $cov;
		}
	} elsif (defined $over_ninty{$trans_id}){
	#} elsif ($cov >= 0){
		} elsif (defined $blast_geneid{$trans_id}){
			if ($blast_geneid{$trans_id} >= $cov){
				$blast_geneid{$trans_id} = $cov;
				$ref_match{$trans_id} = $match;
			}
		} else {
			$blast_geneid{$trans_id} = $cov;
			$ref_match{$trans_id} = $match;
		}
	#}
}
close($gene_list);

my $good = 0;
while (my $line2 = <$pep_seqs>){
	chomp $line2;
	if ($line2 =~ ">"){
	#if ($counter == 0){
		my @names = split / /, $line2;
		my $trans = substr ($names[0], 1);
		if (defined $blast_geneid{$trans}){
			if ($blast_geneid{$trans} >= 50){
				print $ninty_seqs "$line2\n";
				#$good = 1;
			}
			$good = 1;
		} else {
			$good = 0;
		}
	} else {
		if ($good == 1){
			print $ninty_seqs "$line2\n";
		}
	}
}
close($ninty_seqs);
foreach my $trans (sort keys %blast_geneid){
	print $outnum "$trans\t$blast_geneid{$trans}\t$ref_match{$trans}\n";
}
foreach my $ninty_trans (sort keys %over_ninty){
	print $outnum "$ninty_trans\t$over_ninty{$ninty_trans}\t$ref_match{$ninty_trans}\n";
}
close($outnum);