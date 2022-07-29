#!/usr/bin/perl -w                                                                                                                                                                                                                                                                                                                        

use strict;
use warnings;

###Looking for differences in coding peptide sequences for multi-transcript genes###
#open (my $gene_list, "<", "/Users/bhazzard/Documents/PacBioPaperData/Bloodstage_fulllength_geneisoformlist.txt");
open (my $gene_list, "<", "/Users/bhazzard/Documents/PacBioPaperData/Spz_isoform_list.txt");
open (my $pep_list, "<", "/Users/bhazzard/Documents/PacBioPaperData/Spz_fulllength_seqs.pep");

my %iso_genes_list;
while (my $line = <$gene_list>){
	chomp $line;
	my @fields = split /\t/, $line;
	my $gene_id = $fields[0];
	my $trans_id = $fields[1];
	$iso_genes_list{$trans_id} = $gene_id;
}
close($gene_list);

my $counter;
my %pep_gene_list;
my $transcript;
my $gene = "A";
my %gene_trans;
while (my $pep_seqs = <$pep_list>){
	chomp $pep_seqs;
	if ($pep_seqs =~ ">"){
		$gene = "A";
		my @fields = split / /, $pep_seqs;
		my $trans_id = substr ($fields[0], 1);
		#$trans_id =~ s/~.*//;
		#my @gene_info = split /\./, $trans_id;
		#my $gene_id = join (".", $gene_info[0], $gene_info[1]);
		#my $gene_id = substr($trans_id, 0, -2);
		#print "$gene_id\n";
		if (defined $iso_genes_list{$trans_id}){
			$transcript = substr ($fields[0], 1);
			$gene = $iso_genes_list{$trans_id};
			$gene_trans{$gene}{$transcript} = 1;
			#print "$gene\n";
		}
	} else {
		if (length($gene) >= 2) {
			if (defined $pep_gene_list{$gene}{$transcript}){
				my $add_to = $pep_gene_list{$gene}{$transcript};
				#print "$transcript\n$add_to\n";
				$pep_gene_list{$gene}{$transcript} = $add_to . $pep_seqs;
			} else {
				$pep_gene_list{$gene}{$transcript} = $pep_seqs;
			}
		}
	}
}
close ($pep_list);

my %diff_peps;
my %diff_lengths;
my %five_trunc;
my %three_trunc;
foreach my $genes (sort keys %pep_gene_list){
	#print "$genes\n";
	my $first = "A";
	my $trans_diff = 0;
	my $diff_len = 0;
	my $no_met = 0;
	my $no_stop = 0;
	foreach my $trans (sort keys $pep_gene_list{$genes}){
		my $M = substr ($pep_gene_list{$genes}{$trans}, 0, 1);
		my $stop = substr ($pep_gene_list{$genes}{$trans}, -1);
		if ($M eq "M"){
			if ($stop eq "*"){
				if (length($first) >= 2){
					if (length($first) == length($pep_gene_list{$genes}{$trans})){
						if ($pep_gene_list{$genes}{$trans} eq $first){
						} else {
							$trans_diff++;
						}
					} else {
						$diff_len++;
						my $start_first = substr($first, 0, 1);
						my $start_second = substr($pep_gene_list{$genes}{$trans}, 0, 1);
						my $first_piece = substr($first, 0, 10);
						my $second_piece = substr($pep_gene_list{$genes}{$trans}, 0, 10);
						my $first_bulk = substr($first, 0, 50);
						my $second_bulk = substr($pep_gene_list{$genes}{$trans}, 0, 50);
						if ($first_piece eq $second_piece){
							if ($first_bulk eq $second_bulk){
							} else {
								$trans_diff++;
							}
						} else {
							$trans_diff++;
						}
					}
				} else {
					$first = $pep_gene_list{$genes}{$trans};
					#print "$first\t$trans\n$pep_gene_list{$genes}{$trans}\n";
				}
			} else {
				$no_stop++;
			}
		} else {
			$no_met++;
		} 	
	}
	if ($trans_diff >= 1){
		$diff_peps{$genes} = 1;
	}
	if ($diff_len >= 1){
		$diff_lengths{$genes} = $diff_len;
	}	
	if ($no_met >= 1){
		$five_trunc{$genes} = 1;
	}
	if ($no_stop >= 1){
		$three_trunc{$genes} = 1;
	}
}

my $diff_genes = keys %diff_peps;
my $len_diff = keys %diff_lengths;
my $trunc = keys %five_trunc;
my $trunc3 = keys %three_trunc;
open (my $gene_isoforms, ">", "/Users/bhazzard/Documents/full_len_gene_blast_seqdif_isoformlist.txt");
foreach my $pep_genes (sort keys %diff_lengths){
	print $gene_isoforms "$pep_genes\n";
}
close ($gene_isoforms);
print "Genes with sequence length differences = $len_diff\nGenes with coding differences =  $diff_genes\nGenes with 5' Truncation = $trunc\nGenes with 3' Truncation = $trunc3\n";

###For genes with coding differences % with diff start, stop, or exon diff####
open (my $coding_trans, "<", "/Users/bhazzard/Documents/Spz_combined.coding.gtf");

my %transcript_list;
my %exon_num;
my $temp_trans;
my %exon_len;
my %cds_start;
my %cds_end;
my $total;
while (my $all_trans = <$coding_trans>){
	chomp $all_trans;
	if ($all_trans =~ "#"){
	} else {
	my @fields = split /\t/, $all_trans;
	my $descrip = $fields[2]; #mRNA, exon, or CDS
	if ($descrip eq "mRNA"){
		my $info = $fields[8];
		my @collapse_info = split /;/, $info;
		my $trans = substr ($collapse_info[0], 3);
		my $gene = $iso_genes_list{$trans};#substr ($collapse_info[1], 7);
		if (defined $diff_lengths{$gene}){
			#$total++;
			$transcript_list{$trans}= $gene;
			$temp_trans = $trans;
		}
	} elsif ($descrip eq "exon"){
		my $start = $fields[3];
		my $end = $fields[4];
		my $len = int($end - $start);
		my $trans = substr ($fields[8], 7);
		if (defined $transcript_list{$trans}){
			my $exon_gene = $transcript_list{$trans};
			$exon_num{$exon_gene}{$trans}++; #exon number
			$exon_len{$trans}{$start}{$end} = 1; #exon lengths
		}
	} elsif ($descrip eq "CDS"){
		my $start = $fields[3];
		my $end = $fields[4];
		my $len = int($end - $start);
		my $trans = substr ($fields[8], 7);
		if (defined $transcript_list{$trans}){
			my $cds_gene = $transcript_list{$trans};
			if (defined $cds_start{$cds_gene}{$trans}){
			} else {
				$cds_start{$cds_gene}{$trans} = $start; #cds start list
			}
			$cds_end{$cds_gene}{$trans} = $end;	#cds end list
		}
	} else {}
	}
}
close ($coding_trans);

##Alt start
my %diff_cds_start;
my %alt_start_list;
foreach my $start_gene (sort keys %cds_start){
	my $first = 0;
	my $diff_cds_start = 0;
	foreach my $start_trans (sort keys $cds_start{$start_gene}){
		if ($first >= 1){
			if ($first == $cds_start{$start_gene}{$start_trans}){
			} else {
				$diff_cds_start++;
				$alt_start_list{$start_gene}{$start_trans} = 1;
			}
		} else {
			$first = $cds_start{$start_gene}{$start_trans};
		}
	}
	if ($diff_cds_start >= 1){
		$diff_cds_start{$start_gene} = $diff_cds_start;
	}
}
my $cds_start_diff = keys %diff_cds_start;
print "Diff Coding Seq Start = $cds_start_diff\n";

##Alt end
my %diff_cds_end;
my %alt_end_list;
foreach my $end_gene (sort keys %cds_end){
	my $first = 0;
	my $diff_cds_end = 0;
	foreach my $end_trans (sort keys $cds_end{$end_gene}){
		if ($first >= 1){
			if ($first == $cds_end{$end_gene}{$end_trans}){
			} else {
				$diff_cds_end++;
				$alt_end_list{$end_gene}{$end_trans} = 1;
			}
		} else {
			$first = $cds_end{$end_gene}{$end_trans};
		}
	}
	if ($diff_cds_end >= 1){
		$diff_cds_end{$end_gene} = $diff_cds_end;
	}
}
my $cds_end_diff = keys %diff_cds_end;
print "Diff Coding Seq End = $cds_end_diff\n";

##Alt exon num
my %diff_exon_num;
my %alt_exon_list;
foreach my $exon_gene (sort keys %exon_num){
	my $first = 0;
	my $diff_exon = 0;
	foreach my $exon_trans (sort keys $exon_num{$exon_gene}){
		if ($first >= 1){
			if ($first == $exon_num{$exon_gene}{$exon_trans}){
			} else {
				$diff_exon++;
				$alt_exon_list{$exon_gene}{$exon_trans} = 1;
			}
		} else {
			$first = $exon_num{$exon_gene}{$exon_trans};
		}
	}
	if ($diff_exon >= 1){
		$diff_exon_num{$exon_gene} = $diff_exon;
	}
}
my $diff_exon_pat = keys %diff_exon_num;
print "Diff Exon Pattern = $diff_exon_pat\n";

#Print Chart
open (my $type_chart, ">", "/Users/bhazzard/Documents/full_len_blast_gene_isoform_type_chart.txt");

print $type_chart "Gene\tIsoformNumber\tAltCodingStart\tAltCodingEnd\tAltExonPat\n";
foreach my $new_gene_list (sort keys %pep_gene_list){
	my @gene_type;
	push @gene_type, $new_gene_list;
	my $isoform_num =  keys $pep_gene_list{$new_gene_list};
	push @gene_type, $isoform_num;
	if (defined $diff_cds_start{$new_gene_list}){
		if ($diff_cds_start{$new_gene_list} >= 1){
			push @gene_type, $diff_cds_start{$new_gene_list};
		} else {
			push @gene_type, "0";
		}
	} else {
		push @gene_type, "0";
	}
	if (defined $diff_cds_end{$new_gene_list}){
		if ($diff_cds_end{$new_gene_list} >= 1){
			push @gene_type, $diff_cds_end{$new_gene_list};
		} else {
			push @gene_type, "0";
		}
	} else {
		push @gene_type, "0";
	}
	if (defined $diff_exon_num{$new_gene_list}){
		if ($diff_exon_num{$new_gene_list} >= 1){
			push @gene_type, $diff_exon_num{$new_gene_list};
		} else {
			push @gene_type, "0";
		}
	} else {
		push @gene_type, "0";
	}
	my $joined = join ("\t", @gene_type);
	print $type_chart "$joined\n";
}
close ($type_chart);

open (my $alt_start_out, ">", "/Users/bhazzard/Documents/full_len_blast_alt_start_list.txt");
foreach my $start_out_gene (sort keys %alt_start_list){
	foreach my $start_out_iso (sort keys $alt_start_list{$start_out_gene}){
		print $alt_start_out "$start_out_gene\t$start_out_iso\n";
	}
}
close ($alt_start_out);

open (my $alt_end_out, ">", "/Users/bhazzard/Documents/full_len_blast_alt_end_list.txt");
foreach my $end_out_gene (sort keys %alt_end_list){
	foreach my $end_out_iso (sort keys $alt_end_list{$end_out_gene}){
		print $alt_end_out "$end_out_gene\t$end_out_iso\n";
	}
}
close ($alt_end_out);

open (my $alt_exon_out, ">", "/Users/bhazzard/Documents/full_len_blast_alt_exon_list.txt");
foreach my $exon_out_gene (sort keys %alt_exon_list){
	foreach my $exon_out_iso (sort keys $alt_exon_list{$exon_out_gene}){
		print $alt_exon_out "$exon_out_gene\t$exon_out_iso\n";
	}
}
close ($alt_exon_out);

