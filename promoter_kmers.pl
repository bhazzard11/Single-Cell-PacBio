#!/usr/bin/perl -w                                                                                                                                                                                                                                                                                                                        

use strict;
use warnings;

###Looks at promoter region in genome up stream of UTRs and Kmer analysis###
open (my $iso_list, "<", "/Users/bhazzard/Documents/PacBioPaperData/Bloodstage_fulllength_geneisoformlist.txt");
my %isoforms;
while (my $iso_line = <$iso_list>){
	chomp $iso_line;
	my @names = split /\t/, $iso_line;
	my $isos = $names[1];
	$isoforms{$isos} = 1;
}
close($iso_list);
open (my $trans_list, "<", "/Users/bhazzard/Documents/PacBioPaperData/Bloodstage_fulllength_seqs.pep");
my %transcripts;
while (my $list_line = <$trans_list>){
	chomp $list_line;
	if ($list_line =~ ">"){
		my @fields = split / /, $list_line;
		my $name = substr($fields[0], 1);
		if (defined $isoforms{$name}){
		} else {
			$transcripts{$name} = 1;
		}
	}
}
close($trans_list);
print "finding transcripts...";
my %trans_promoter_start;
my %trans_promoter_end;
open (my $transdecoder, "<", "/Users/bhazzard/Documents/Bloodstage_combined.coding.gtf");
while (my $new_ref_line = <$transdecoder>){
	chomp $new_ref_line;
 	my @fields = split /\t/, $new_ref_line;
    if (scalar(@fields) > 3){
    	my $chr = $fields[0];
        my $ident = $fields[2];
        my $info = $fields[8];
        my $start = $fields[3];
        my $end = $fields[4];
        my $strand = $fields[6];
        if ($ident eq "mRNA"){
        	my @name = split /;/, $info;
        	my $trans_name = substr($name[0], 3);
        	if ($strand eq "+"){
	        	my $promoter_start = $start - 50;
    	    	$trans_promoter_start{$chr}{$trans_name} = $promoter_start;
        		$trans_promoter_end{$chr}{$trans_name} = $start;
        	} elsif ($strand eq "-"){
        		#my $promoter_start = $end + 100;
    	    	#$trans_promoter_start{$chr}{$trans_name} = $end;
        		#$trans_promoter_end{$chr}{$trans_name} = $promoter_start;
        	}
        }
    }
}

open (my $genome, "<", "/Users/bhazzard/Downloads/PlasmoDB-54_PvivaxP01_Genome.fasta");
my $genome_chr;
my %sequences;
print "complete\nreading genome...\n";
while (my $fasta = <$genome>){
	chomp $fasta;
	if ($fasta =~ ">"){
		my @chr_name = split / /, $fasta;
		$genome_chr = substr($chr_name[0], 1);
		print "$genome_chr\n";
	} else {
		if (defined $sequences{$genome_chr}){
			my $old_seq = $sequences{$genome_chr};
			my $new_seq = $old_seq.$fasta;
			$sequences{$genome_chr} = $new_seq;
		} else {
			$sequences{$genome_chr} = $fasta;
		}
	}
}
print "complete\n";
print "finding promoters...\n";
open (my $out_fasta, ">", "/Users/bhazzard/Documents/promoter_100bp.fasta");
my %promoter_seqs;
foreach my $chr (sort keys %trans_promoter_start){
	foreach my $trans (sort keys $trans_promoter_start{$chr}){
		my $promo_start = $trans_promoter_start{$chr}{$trans};
		my $promo_end = $trans_promoter_end{$chr}{$trans};
		my $promo_seq = substr($sequences{$chr}, $promo_start, 100);
		print $out_fasta ">$trans\n$promo_seq\n";
		$promoter_seqs{$trans} = $promo_seq;
	}
}	
close($genome);
close($out_fasta);

print "complete\nfinding kmers...";

my %kmers_list;
foreach my $trans_ID (sort keys %promoter_seqs){
	my $seq = $promoter_seqs{$trans_ID};
	my $len = length($seq);
	my $loop = 0;
	while ($loop < ($len - 5)){
		my $kmer = substr($seq, $loop, 5);
		if (defined $kmers_list{$kmer}){
			$kmers_list{$kmer}++;
		} else {
			$kmers_list{$kmer} = 1;
		}
		$loop++;
	}
}
print "complete\nprinting to file\n";
open (my $output, ">", "/Users/bhazzard/Documents/promoter_kmers_100bp.txt");
print $output "kmer\tnumber\n";
foreach my $kmer (sort keys %kmers_list){
	print $output "$kmer\t$kmers_list{$kmer}\n";
}
close($output);