#!/usr/bin/perl -w                                                                                                                                                                                                                                                                                                                        

use strict;
use warnings;

##make fasta files of 5' and 3' UTRs and do kmer analaysis##
open (my $full_length, "<", "/Users/bhazzard/Documents/Bloodstage_combined.coding.gtf");
open (my $gff, "<", "/Users/bhazzard/Documents/Bloodstage_combined.alltranscripts.fasta.transdecoder.gff3");
open (my $fasta_in, "<", "/Users/bhazzard/Documents/PacBioPaperData/Bloodstage_combined.alltranscripts.fasta");
open (my $iso_list, "<", "/Users/bhazzard/Documents/PacBioPaperData/Bloodstage_fulllength_geneisoformlist.txt");
my %isoforms;
while (my $iso_line = <$iso_list>){
	chomp $iso_line;
	my @names = split /\t/, $iso_line;
	my $isos = $names[1];
	$isoforms{$isos} = 1;
}

my %gene_list;
while (my $list = <$full_length>){
	chomp $list;
	my @fields = split /\t/, $list;
	my $type = $fields[2];
	if ($type eq "mRNA"){
		my @info = split /;/, $fields[8];
		my $fullname = substr($info[0], 3);
		if (defined $isoforms{$fullname}){
		} else {
			my @name = split /\./, $fullname;
			my $trans = join(".", $name[0], $name[1], $name[2]);
			$gene_list{$trans} = 1;
		}
	}
}
close($full_length);

my %three_cords;
my %five_cords;
while (my $locations = <$gff>){
	chomp $locations;
	my @fields = split /\t/, $locations;
	if (scalar @fields >= 8){
		my $trans = $fields[0];
		my $type = $fields[2];
		my $start = $fields[3];
		my $end = $fields[4];
		my $strand = $fields[6];
		#print "$trans\n$strand\n";
		if ($strand eq "+"){
			if (defined $gene_list{$trans}){
				if ($type eq "five_prime_UTR"){
					$five_cords{$trans} = $end;
				} elsif ($type eq "three_prime_UTR"){
					$three_cords{$trans} = $start;
				}	
			}
		}
	}
}
close($gff);

open (my $three_fasta, ">", "/Users/bhazzard/Documents/PacBioPaperData/Bloodstage_3UTR.fasta");
open (my $five_fasta, ">", "/Users/bhazzard/Documents/PacBioPaperData/Bloodstage_5UTR.fasta");
my %three_seqs;
my %five_seqs;
my $counter = 0;
my $trans_name;
while (my $sequences = <$fasta_in>){
	chomp $sequences;
	if ($counter == 0){
		my @names = split / /, $sequences;
		$trans_name = substr($names[0], 1);
		$counter++;
	} else {
		my $seq = $sequences;
		if (defined $three_cords{$trans_name}){
			my $three_start = $three_cords{$trans_name};
			my $three = substr($seq, $three_start);
			$three_seqs{$trans_name} = $three;
			if (length($three) >= 10){
				print $three_fasta ">$trans_name\n$three\n";
			}
		}
		if (defined $five_cords{$trans_name}){
			my $five_end = $five_cords{$trans_name};
			my $five = substr($seq, 0, $five_end);
			$five_seqs{$trans_name} = $five;
			#print "$trans_name\n$seq\n$five_end\n$five\n";
			print $five_fasta ">$trans_name\n$five\n";
		}
		$counter = 0;
	}	
}
close($fasta_in);
close($three_fasta);
close($five_fasta);

open (my $three_kmers_out, ">", "/Users/bhazzard/Documents/PacBioPaperData/Bloodstage_3UTRlast100bp_kmers.txt");
open (my $five_kmers_out, ">", "/Users/bhazzard/Documents/PacBioPaperData/Bloodstage_5UTRlast100bp_kmers.txt");

my %three_kmers;
foreach my $three_trans (sort keys %three_seqs){
	my $seq = $three_seqs{$three_trans};
	my $len = length($seq);
	if ($len >= 100){
		$seq = substr($seq, -100);
		$len = length($seq);
	}
	my $loop = 0;
	while ($loop < ($len - 5)){
		my $kmer = substr($seq, $loop, 5);
		if (defined $three_kmers{$kmer}){
			$three_kmers{$kmer}++;
		} else {
			$three_kmers{$kmer} = 1;
		}
		$loop++;
	}
}
my %five_kmers;
foreach my $five_trans (sort keys %five_seqs){
	my $seq = $five_seqs{$five_trans};
	my $len = length($seq);
	if ($len >= 100){
		$seq = substr($seq, -100);
		$len = length($seq);
	}
	my $loop = 0;
	while ($loop < ($len - 5)){
		my $kmer = substr($seq, $loop, 5);
		if (defined $five_kmers{$kmer}){
			$five_kmers{$kmer}++;
		} else {
			$five_kmers{$kmer} = 1;
		}
		$loop++;
	}
}

foreach my $three_ks (sort keys %three_kmers){
	print $three_kmers_out "$three_ks\t$three_kmers{$three_ks}\n";
}
close($three_kmers_out);
foreach my $five_ks (sort keys %five_kmers){
	print $five_kmers_out "$five_ks\t$five_kmers{$five_ks}\n";
}
close($five_kmers_out);

