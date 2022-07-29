#!/usr/bin/perl -w                                                                                                                                                                                                                                                                                                                        

use strict;
use warnings;

open (my $genelist, "<", "/Users/bhazzard/Documents/PacBioPaperData/Spz_isoform_list.txt");
my %isoforms;
my %transcripts;
while (my $iso_list = <$genelist>){
	chomp $iso_list;
	my @fields = split /\t/, $iso_list;
	my $name = $fields[0];
	my $trans = $fields[1];
	$isoforms{$name}{$trans} = 1;
	$transcripts{$trans} = 1;
}
close($genelist);

open (my $new_list, "<", "/Users/bhazzard/Documents/PacBioPaperData/Spz_fulllength_genes_with_utr_exons.txt");

my %exon;
while (my $utr_exon = <$new_list>){
	chomp $utr_exon;
	my @fields = split /\t/, $utr_exon;
	my $trans = $fields[0];
	if (defined $transcripts{$trans}){
		$exon{$trans} = 1;
	}
}
close($new_list);

my $change;
my %utr_gene_exon;
foreach my $gene (sort keys %isoforms){
	my $isos = 0;
	my $iso_exon = 0;
	foreach my $trans (sort keys $isoforms{$gene}){
		$isos++;
		if (defined $exon{$trans}){
			$iso_exon++;
			$utr_gene_exon{$gene} = 1;
		}
	}
	print "$gene\t$isos\t$iso_exon\n";
	if ($isos == $iso_exon){
	} elsif ($iso_exon == 0) {
	} else {
		$change++;
	}
}

print "UTR exon isoform = $change\n";


open (my $output, ">", "/Users/bhazzard/Documents/PacBioPaperData/Spz_genes_with_UTR_exon.txt");
foreach my $gene_utr (sort keys %utr_gene_exon){
	print $output "$gene_utr\n";
}
