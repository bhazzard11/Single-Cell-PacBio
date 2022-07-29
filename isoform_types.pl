#!/usr/bin/perl -w                                                                                                                                                                                                                                                                                                                        

use strict;
use warnings;

###Determines the type of differences between multi-transcript genes: 5' and 3' UTR and Exon number#
open (my $coding_trans, "<", "/Users/bhazzard/Documents/Spz_combined.coding.gtf");
open (my $gene_list, "<", "/Users/bhazzard/Documents/PacBioPaperData/Spz_fulllength_geneisoformlist.txt");
#open (my $gene_list, "<", "/Users/bhazzard/Documents/PacBioPaperData/Bloodstage_isoform_list.txt");


my %iso_genes_list;
while (my $line = <$gene_list>){
	chomp $line;
	my @fields = split /\t/, $line;
	my $gene_id = $fields[0];
	my $trans = $fields[1];
	#$iso_genes_list{$gene_id}{$trans} = 1;
	$iso_genes_list{$trans} = $gene_id;
}
close($gene_list);

my %transcript_list;
my %trans_start;
my %trans_end;
my %exon_num;
my $temp_trans;
my %exon_len;
my %cds_start;
my %cds_len;
my $total;
while (my $all_trans = <$coding_trans>){
	chomp $all_trans;
	my @fields = split /\t/, $all_trans;
	my $descrip = $fields[2]; #mRNA, exon, or CDS
	if (scalar(@fields) >= 2){
	if ($descrip eq "mRNA"){
		my $start = $fields[3];
		my $end = $fields[4];
		my $len = int($end - $start);
		my $info = $fields[8];
		my @collapse_info = split /;/, $info;
		my $trans = substr ($collapse_info[0], 3);
		my $gene = $iso_genes_list{$trans};#substr ($collapse_info[1], 7);
		if (defined $iso_genes_list{$trans}){
			#$total++;
			$transcript_list{$gene}{$trans} = $len;
			$trans_start{$trans} = $start;
			$trans_end{$trans} = $end;
			$temp_trans = $trans;
		}
	} elsif ($descrip eq "exon"){
		my $start = $fields[3];
		my $end = $fields[4];
		my $len = int($end - $start);
		my $trans = substr ($fields[8], 7);
		$exon_num{$trans}++; #exon number
		if ($temp_trans = $trans){
			$exon_len{$trans}{$start}{$end} = 1; #exon lengths
		}
	} elsif ($descrip eq "CDS"){
		my $start = $fields[3];
		my $end = $fields[4];
		my $len = int($end - $start);
		my $trans = substr ($fields[8], 7);
		if ($temp_trans = $trans){
			if (defined $cds_start{$trans}){
				$cds_len{$trans} = int($cds_len{$trans} + $len); #cds length
			} else {
				$cds_start{$trans} = $start; #cds start can be used with length to determine end then UTRs
				$cds_len{$trans} = $len;
			}
		}
	} else {}
}
}
close($coding_trans);

my %alt_five_UTR;
my %alt_three_UTR;
my %dif_exon_num;
foreach my $gene (sort keys %transcript_list){
	#$total++;
	my $gene_exon_num = 0;
	my $cds_start_num = 0;
	my $cds_end_num = 0;
	foreach my $trans (sort keys $transcript_list{$gene}){
		$total++;
		###exon number difference
		my $trans_exon_num = $exon_num{$trans};
		if ($gene_exon_num >= 1){
			#print "$trans\t$gene_exon_num\t$trans_exon_num\n";
			if ($gene_exon_num == $trans_exon_num){
				#$dif_exon_num++;
			} else {
				$dif_exon_num{$gene}++;
			}
		} else {
			$gene_exon_num = $exon_num{$trans};
		}
		##UTR differences
		my $start_cds = $cds_start{$trans}; #cds start positon
		my $end_cds = int($cds_start{$trans} + $cds_len{$trans}); #cds end position
		my $five = int($start_cds - $trans_start{$trans}); #5' length
		my $three = int($trans_end{$trans} - $end_cds); #3' length
		###alt 5' UTR
		if ($cds_start_num >= 10){ ##alt 5' UTRs
			if ($cds_start_num != $five){
				$alt_five_UTR{$gene}{$trans} = 1;
			} else {}
		} else {
			$cds_start_num = $five;
		}
		###alt 3' UTR
		if ($cds_end_num >= 10){ ##alt 3' UTRs
			if ($cds_end_num != $three) {
				$alt_three_UTR{$gene}{$trans} = 1;
			} else {}
		} else {
			$cds_end_num = $three;
		}
	}
}

my $exon_count = keys %dif_exon_num;
my $five_count = keys %alt_five_UTR;
my $three_count = keys %alt_three_UTR;
print "Total Transcripts= $total\nExon number difference = $exon_count\nAlt 5' Length = $five_count\nAlt 3' Length = $three_count\n";

my %alt_either_utr;
my %exon_and_utr;
open (my $alt_fiveUTR_out, ">", "/Users/bhazzard/Documents/alt_fiveUTR_list.txt");
foreach my $start_fiveUTR_gene (sort keys %alt_five_UTR){
	foreach my $start_fiveUTR_iso (sort keys $alt_five_UTR{$start_fiveUTR_gene}){
		$alt_either_utr{$start_fiveUTR_gene} = 1;
		print $alt_fiveUTR_out "$start_fiveUTR_gene\t$start_fiveUTR_iso\n";
		if (defined $dif_exon_num{$start_fiveUTR_gene}){
			$exon_and_utr{$start_fiveUTR_gene} = 1;
		}
	}
}
close ($alt_fiveUTR_out);

open (my $alt_threeUTR_out, ">", "/Users/bhazzard/Documents/alt_threeUTR_list.txt");
foreach my $start_threeUTR_gene (sort keys %alt_three_UTR){
	foreach my $start_threeUTR_iso (sort keys $alt_three_UTR{$start_threeUTR_gene}){
		print $alt_threeUTR_out "$start_threeUTR_gene\t$start_threeUTR_iso\n";
		if (defined $alt_either_utr{$start_threeUTR_gene}){
		} else {
			$alt_either_utr{$start_threeUTR_gene} = 1;
		}
		if (defined $dif_exon_num{$start_threeUTR_gene}){
			$exon_and_utr{$start_threeUTR_gene} = 1;
		}
	}
}
close ($alt_threeUTR_out);
my $all_utr_isos = keys %alt_either_utr;
my $all_exon_inUTR = keys %exon_and_utr;
print "Alt UTR genes = $all_utr_isos\nExon in UTR = $all_exon_inUTR\n";
