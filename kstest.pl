#!/usr/bin/perl -w                                                                                                                                                                                                                                                                                                                        

use strict;
use warnings;

####runs kstest on genes with isoforms to determine genes that are different over pseudotime####
####Outputs matrix of d and d-crit values####

#input 1: tab delimited file of list of genes with isoforms and the isoforms names
open (my $genelist, "<", "/path/to/isoform_list.txt");
#input 2: counts table of expression values
open (my $countmatrix, "<", "/path/to/transcript/expression/counts.txt");
#output 1: table of d and d-crit values
open (my $output, ">", "/path/to/iso_kstest_out.txt");

#make transcript list
my %genes;
my %transname;
while (my $list = <$genelist>){
	chomp $list;
	my @fields = split /\t/, $list;
	my $gene = $fields[0];
	my $trans = substr($fields[1],0, -3);
	$genes{$gene}{$trans} = 1;
	$transname{$trans} = 1;
}
close($genelist);
print "gene list made\n";

#make counts of transcripts
my $counter = 0;
my @pseudotime;
my %ks;
my %trans_n;
my %cells_per_iso;
my %cell_list;
while (my $matrix = <$countmatrix>){
	chomp $matrix;
	if ($counter == 0){
		@pseudotime = split /\t/, $matrix;
		$counter++;
	} else {
		my @cell = split /\t/, $matrix; #split line with transcript values into individual values
		my $transcript = $cell[0]; #capture transcript name
		if (defined $transname{$transcript}){
			my $value = 0;
			my $count = 1;
			my $cum = 0;
			my $cell_num = 0;
			while ($count <= 2718) {
				if ($cell[$count] >= 1){
					$cum = ($value + $cell[$count]);
					$value = $cell[$count];
					$ks{$transcript}{$pseudotime[$count - 1]} = $cell_num;
					$cell_num++;
					$cell_list{$pseudotime[$count - 1]}{$transcript}++;
				} else {
					$ks{$transcript}{$pseudotime[$count - 1]} = $cell_num;
				}
				$count++;
			}
			$trans_n{$transcript} = $cell_num;
			$cells_per_iso{$transcript} = $cell_num;
		}
	}
}
close($countmatrix);
print "matrix input\n";

#consider only genes where more than one transcript has it in 50 cells
my $iso_num = keys %cells_per_iso;
my $total_iso_cells = 0;
my %fifty;
my %cell_iso_num;
foreach my $iso_trans (sort keys %cells_per_iso){
	$total_iso_cells = ($total_iso_cells + $cells_per_iso{$iso_trans});
	if ($cells_per_iso{$iso_trans} >= 50){
		$fifty{$iso_trans} = 1;
	}
}
my %gene_fifty;
foreach my $gene_iso (sort keys %genes){
	foreach my $trans_iso (sort keys $genes{$gene_iso}){
		if (defined $fifty{$trans_iso}){
			$gene_fifty{$gene_iso}++;
		}
	}
	foreach my $ind_cell (sort keys %cell_list){
		my $iso_count = 0;
		foreach my $cell_trans (sort keys $genes{$gene_iso}){
			if (defined $cell_list{$ind_cell}{$cell_trans}){
				$iso_count++;
			}
		}
		$cell_iso_num{$ind_cell}{$gene_iso} = $iso_count;
	}
}
my $fifty_isos = 0;
my %gene_iso_list;
foreach my $gene_fifty_count (sort keys %gene_fifty){
	if ($gene_fifty{$gene_fifty_count} >= 2){
		$gene_iso_list{$gene_fifty_count} = 1;
		$fifty_isos++;
	}
}

#print number of genes that passed QC
my $total_genes = keys %genes;
my $avg_cells = ($total_iso_cells/$iso_num);
print "Total Genes = $total_genes\nTotal Isoforms = $iso_num\nAverage Cells Per Iso = $avg_cells\nGenes with at least 50 cells in >2 isoforms = $fifty_isos\n";

#per cell info
my $cells_with_isos;
my $cell_total;
my %avg_gene_iso;
foreach my $cell (sort keys %cell_iso_num){
	$cell_total++;
	my $isogene = 0;
	foreach my $cell_iso (sort keys $cell_iso_num{$cell}){
		if ($cell_iso_num{$cell}{$cell_iso} >= 2){
			$isogene++;
		}
	}
	if ($isogene >= 1){
		$cells_with_isos++;
		$avg_gene_iso{$cell} = $isogene;
	}
}
my $sum_cells = 0;
foreach my $cell_avg (sort keys %avg_gene_iso){
	$sum_cells = ($sum_cells + $avg_gene_iso{$cell_avg});
}

my $avg_iso_cell = ($sum_cells / $cell_total);
print "Total Cells = $cell_total\nCells with more than one isoform = $cells_with_isos\nAverage isoform presence per cell = $avg_iso_cell\n";

#calculate n and d values for each gene
my %gene_pseudo;
my %gene_n;
foreach my $gene (sort keys %gene_iso_list){
	my $counter = 0;
	my $isoforms = keys %{$genes{$gene}};
	my %sample_d;
	my $iso_num = 1;
	foreach my $trans (sort keys %{$genes{$gene}}){
		my $S = 0;
		foreach my $pseudo_value (sort keys %{$ks{$trans}}){
			if ($trans_n{$trans} == 0){
			} else {
				my $dvalue = ($ks{$trans}{$pseudo_value}/$trans_n{$trans});
				if ($dvalue >= $S){
					$S = $dvalue;
				}
			#}
				$sample_d{$iso_num}{$pseudo_value} = $S;
				my $gene_p = join ("_",$gene, $pseudo_value);
				my $iso_p = join("_", $trans, $pseudo_value);
				$gene_pseudo{$gene_p}{$iso_p} = $S;
			}
		}
		$iso_num++;
		$gene_n{$gene}{$trans} = $trans_n{$trans};
	}
}
print "d values found\n";

my %gene_D;
my %max_n_list;
my %min_n_list;
foreach my $gene_num (sort keys %gene_pseudo){
	my @trans_values;
	foreach my $trans_num (sort keys $gene_pseudo{$gene_num}){
		my $d_val = $gene_pseudo{$gene_num}{$trans_num};
		push @trans_values, $d_val;
	}
	my ($min, $max) = (sort {$a <=> $b} @trans_values)[0,-1];
	my $D = ($max - $min);
	my @split_name = split /\_/, $gene_num;
	my $gene_name = join("_", $split_name[0], $split_name[1]);
	$gene_D{$gene_name}{$gene_num} = $D;
}

my %sig_genes;
foreach my $gene_names (sort keys %gene_D){
	my $max = 0;
	foreach my $pseudo (sort keys $gene_D{$gene_names}){
		my $max_D = $gene_D{$gene_names}{$pseudo};
		if ($max_D > $max){
			$max = $max_D;
		}
	}
	$sig_genes{$gene_names} = $max;
}

#write to file
foreach my $sig_gene (sort keys %sig_genes){
	foreach my $sig_trans (sort keys $genes{$sig_gene}){
		if ($trans_n{$sig_trans} >= 50){
			print $output "$sig_gene\t$sig_trans\t$trans_n{$sig_trans}\t$sig_genes{$sig_gene}\n";
		}
	}
}
