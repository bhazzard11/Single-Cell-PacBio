#!/usr/bin/perl -w                                                                                                                                                                                                                                                                                                                        

use strict;
use warnings;

open (my $blast, "<", "/Users/bhazzard/Documents/Blast_v54/Blood_vs_Spz.tab");
open (my $trans_list, "<", "/Users/bhazzard/Documents/PacBioPaperData/Spz_fulllength_seqs.pep");
open (my $output, ">", "/Users/bhazzard/Documents/Spz_comparison_blood.txt");
open (my $transdecoder, "<", "/Users/bhazzard/Documents/Spz_combined.coding.gtf");

my $counter = 0;
my %transcripts;
while (my $list_line = <$trans_list>){
	chomp $list_line;
	if ($list_line =~ ">"){
		my @fields = split / /, $list_line;
		my $name = substr($fields[0], 1);
		$transcripts{$name} = 1;
		$counter++;
	}
}
close($trans_list);
print "Total Transcripts = $counter\n";

my %ref_id;
my %ref_cov;
my %new_ref;
my %cds_len;
while (my $line = <$blast>){
	chomp $line;
	my @fields = split /\t/, $line;
	my $transcript = $fields[0];
	my $gene = $fields[1];
	my $per_id = $fields[2];
	my $align_len = $fields[3];
	my $ref_len = $fields[4];
	my $trans_len = $fields[6];
	my $per_cov = "NA";#($align_len/$ref_len)*100;
	if (defined $ref_id{$transcript}){
	} else {
		$ref_id{$transcript} = $per_id;
		$new_ref{$transcript} = $gene;
	}
	if (defined $ref_cov{$transcript}){
	} else {
		$ref_cov{$transcript} = $per_cov;
		$cds_len{$transcript} = $trans_len;
	}	
}

close($blast);

my %new_ref_cds;
while (my $new_ref_line = <$transdecoder>){
    chomp $new_ref_line;
    my @fields = split /\t/, $new_ref_line;
    if (scalar(@fields) > 3){
        my $ident = $fields[2];
        my $start = $fields[3];
        my $end = $fields[4];
        my $len = int($end-$start);
        my $info = $fields[8];
        if ($ident eq "CDS"){
            my @collapse_info = split /;/, $info;
            my $gene_id = substr($collapse_info[0], 7);
            if (defined $new_ref_cds{$gene_id}){
                my $old_len = $new_ref_cds{$gene_id};
                my $new_len = int($old_len + $len); #add CDS lengths                                                                                                             \
                                                                                                                                                                                  
                $new_ref_cds{$gene_id} = $new_len;
#               print "$gene_id\t$ref_cds{$gene_id}\n";                                                                                                                           
            } else {
                $new_ref_cds{$gene_id} = $len;
            }
        }
    }
}

my $great_ninty_id = 0;
my $ninty_to_thirty_id = 0;
my $less_thirty_id = 0;
foreach my $trans (sort keys %transcripts){
	if (defined $ref_id{$trans}){
		if ($ref_id{$trans} >= 90){
			$great_ninty_id++;
		} elsif ($ref_id{$trans} < 90){
			if ($ref_id{$trans} > 30){
				$ninty_to_thirty_id++;
			} elsif ($ref_id{$trans} < 30){
				$less_thirty_id++;
			}
		}
	} else {
		$less_thirty_id++;
	}	
}

print "Percent Identity\n>90 = $great_ninty_id\n90 to 30 = $ninty_to_thirty_id\nLess than 30 = $less_thirty_id\n";
print $output "Transcript_ID\tGene_ID\tGene_Name\tPercent_Identity\tCoding_Seq_Length\tPercent_Coverage\n";

my $great_ninty_cov = 0;
my $ninty_to_thirty_cov = 0;
my $less_thirty_cov = 0;
foreach my $trans_cov (sort keys %transcripts){
	if (defined $ref_cov{$trans_cov}){
		print $output "$trans_cov\t$transcripts{$trans_cov}\t$new_ref{$trans_cov}\t$ref_id{$trans_cov}\t$new_ref_cds{$trans_cov}\t$ref_cov{$trans_cov}\n";
		if ($ref_cov{$trans_cov} >= 90){	
			$great_ninty_cov++;
		} elsif ($ref_cov{$trans_cov} < 90){
			if ($ref_cov{$trans_cov} > 30){
				$ninty_to_thirty_cov++;
			} elsif ($ref_cov{$trans_cov} < 30){
				$less_thirty_cov++;
			}
		}	
	} else {
		print $output "$trans_cov\t$transcripts{$trans_cov}\tNA\t0\t$new_ref_cds{$trans_cov}\t0\n";
		$less_thirty_cov++;
	}
}
close($output);
print "Percent Coverage\n>90 = $great_ninty_cov\n90 to 30 = $ninty_to_thirty_cov\nLess than 30 = $less_thirty_cov\n";
