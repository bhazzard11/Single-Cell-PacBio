#!/usr/bin/perl -w                                                                                                                                                                                                                                                                                                                        

use strict;
use warnings;

###Counts number of introns in coding region and compares to number of introns in UTR###


open (my $trans_list, "<", "/Users/bhazzard/Documents/PacBioPaperData/Bloodstage_fulllength_seqs.pep"); ##peptide list
my %transcripts;
while (my $list_line = <$trans_list>){
	chomp $list_line;
	if ($list_line =~ ">"){
		my @fields = split / /, $list_line;
		my $name = substr($fields[0], 1);
		$transcripts{$name} = 1;
	}
}
close($trans_list);

open (my $transdecoder, "<", "/Users/bhazzard/Documents/Bloodstage_combined.coding.gtf");
my $exon_count;
my $cds_count;
my %exon_num;
my %cds_num;
while (my $new_ref_line = <$transdecoder>){
	chomp $new_ref_line;
 	my @fields = split /\t/, $new_ref_line;
    if (scalar(@fields) > 3){
        my $ident = $fields[2];
        my $info = $fields[8];
        my $start = $fields[3];
        my $end = $fields[4];
        if ($ident eq "exon"){
        	my $trans_id = substr($info, 7);
        	if (defined $transcripts{$trans_id}){
        		if (defined $exon_num{$trans_id}){
        			$exon_num{$trans_id}++;
        		} else {
        			$exon_num{$trans_id} = 1;
        		}
        	}
        } elsif ($ident eq "CDS"){
        	my $trans_id = substr($info, 7);
        	if (defined $transcripts{$trans_id}){
        		if (defined $cds_num{$trans_id}){
        			$cds_num{$trans_id}++;
        		} else {
	        		$cds_num{$trans_id} = 1;
	        	}
	        }
	    }
	}
}
close($transdecoder);

open (my $output, ">", "/Users/bhazzard/Documents/Bloodstage_utr_exon_vs_coding_exon.txt");
print $output "TransID\tCodingExons\tUTRExons\n";
foreach my $trans (sort keys %transcripts){
	if (defined $exon_num{$trans}){
		if ($exon_num{$trans} > $cds_num{$trans}){
			my $utr_exons = int($exon_num{$trans} - $cds_num{$trans});
			print $output "$trans\t$cds_num{$trans}\t$utr_exons\n";
		} else {
			print $output "$trans\t$cds_num{$trans}\t0\n";
		}
	}
}
close($output);

my $coding_one_no_utr;
my $coding_one_yes_utr;
my $coding_intron_no_utr;
my $coding_intron_yes_utr;
open (my $output2, ">", "/Users/bhazzard/Documents/Bloodstage_utr_number_comp.txt");
foreach my $trans (sort keys %transcripts){
	if (defined $exon_num{$trans}){
		if ($cds_num{$trans} == 1){
			if ($exon_num{$trans} > $cds_num{$trans}){
				$coding_one_yes_utr++;
			} else {
				$coding_one_no_utr++;
			}
		} else {
			if ($exon_num{$trans} > $cds_num{$trans}){
				$coding_intron_yes_utr++;
			} else {
				$coding_intron_no_utr++;
			}
		}
	}
}
print $output2 "0UTR\t>1UTR\n1Coding\t$coding_one_no_utr\t$coding_one_yes_utr\n>1Coding\t$coding_intron_no_utr\t$coding_intron_yes_utr\n";
close($output2);