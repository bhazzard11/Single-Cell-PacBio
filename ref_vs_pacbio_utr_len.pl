#!/usr/bin/perl -w                                                                                                                                                                                                                                                                                                                        

use strict;
use warnings;

###Finds UTR lengths of reference and prints to file with comparison lengths from predictions###
open (my $ref_trans, "<", "/Users/bhazzard/Downloads/PlasmoDB-54_PvivaxP01.gff");
open (my $coding_trans, "<", "/Users/bhazzard/Documents/Spz_combined.alltranscripts.fasta.transdecoder.gff3");
open (my $ref_match, "<", "/Users/bhazzard/Documents/Blast_v54/Spz_per_id_over_fifty.txt");
open (my $output, ">", "/Users/bhazzard/Documents/Ref_vs_spz_UTR_lengths.txt");

my %ref_five;
my %ref_three;
my %trans_list;
while (my $line1 = <$ref_trans>){
	chomp $line1;
	my @fields = split /\t/, $line1;
	if (scalar(@fields) >= 8){
		my $ident = $fields[2];
		my $start = $fields[3];
		my $end = $fields[4];
		my $length = int($end - $start);
		my @name_mess = split /;/, $fields[8];
		if ($ident eq "five_prime_UTR"){
			my $name = substr($name_mess[1], 7);
			if (defined $ref_five{$name}){
				my $new_len = $length + $ref_five{$name};
				$ref_five{$name} = $length;
			} else {
				$ref_five{$name} = $length;
			}
		} elsif ($ident eq "three_prime_UTR"){
			my $name = substr($name_mess[1], 7);
			if (defined $ref_three{$name}){
				my $new_len = $length + $ref_three{$name};
				$ref_three{$name} = $length;
			} else {
				$ref_three{$name} = $length;
			}
			$ref_three{$name} = $length;
		} elsif ($ident eq "CDS"){
			my $name = substr($name_mess[1], 7);
			$trans_list{$name} = "REF";
		}
	}
}	
close($ref_trans);

my %check_list;
while (my $fasta = <$ref_match>){
	chomp $fasta;
	my @fields = split /\t/, $fasta;
	my $mine = $fields[0];
	my @ref_full = split /\-/, $fields[2];
	my $ref = $ref_full[0];
	my $per_id = $fields[1];
	if ($per_id >= 90){
		$check_list{$mine} = $ref;
	}
}

my %five_list;
my %three_list;
while (my $line2 = <$coding_trans>){
	chomp $line2;
	my @fields = split /\t/, $line2;
	if (scalar(@fields) >= 8){
		my $ident = $fields[2];
		my $start = $fields[3];
		my $end = $fields[4];
		my $length = int($end - $start);
		my @name_mess = split /;/, $fields[8];
		if ($ident eq "five_prime_UTR"){
			my $name = substr($name_mess[1], 7);
			$five_list{$name} = $length;
		} elsif ($ident eq "three_prime_UTR"){
			my $name = substr($name_mess[1], 7);
			$three_list{$name} = $length;
		} elsif ($ident eq "CDS"){
			my $name = substr($name_mess[1], 7);
			if (defined $check_list{$name}){
				$trans_list{$name} = "Blood";
			}
		}
	}
}	
close($coding_trans);

print $output "TranscriptID\tRef_ID\tFiveUTRPB\tFiveUTRREF\tThreeUTRPB\tThreeUTRREF\n";
foreach my $trans (sort keys %check_list){
	my $ref = $check_list{$trans};
	my @data;
	if (defined $five_list{$trans}){
		push @data, $five_list{$trans};
	} else {
		push @data, "0";
	} if (defined $ref_five{$ref}){
		push @data, $ref_five{$ref};
	} else {
		push @data, "0";
	}
	if (defined $three_list{$trans}){
		push @data, $three_list{$trans};
	} else {
		push @data, "0";
	} if (defined $ref_three{$ref}){
		push @data, $ref_three{$ref};
	} else {
		push @data, "0";
	}
	print $output join("\t", $trans, $ref, @data), "\n";
}
close($output);