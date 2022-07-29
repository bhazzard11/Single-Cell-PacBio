#!/usr/bin/perl -w                                                                                                                                                                                                                                                                                                                        

use strict;
use warnings;

###find and count transcripts with utr exons###
#input 1: peptide sequences list for names
open (my $trans_list, "<", "/path/to/peptide/fasta/transcripts.pep");
#output 1: list of transcripts with introns in the utrs
open (my $new_list, ">", "/path/to/output_fulllength_genes_with_utr_exons.txt");
#input 2: gtf of coding transcripts
open (my $transdecoder, "<", "/path/to/protein/list.coding.gtf");
#input 3: list of genes with isoforms
open (my $genelist, "<", "/path/to/genes/with/isoform/list/fulllength_geneisoformlist.txt");

#make transcript list, only uses names.
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

#Looks through gtf and finds exon designated transcripts and compares the number of exons per gene to the number of CDS exons
my %mRNAs_start;
my %mRNAs_end;
my %exons_start;
my %exons_end;
my $exon_num = 1;
while (my $new_ref_line = <$transdecoder>){
    chomp $new_ref_line;
    my @fields = split /\t/, $new_ref_line;
    if (scalar(@fields) > 3){
        my $ident = $fields[2]; #gtf designation, "mRNA", "exon", "CDS"
        my $info = $fields[8]; #gene name is in here
        my $start = $fields[3]; #location of start of designation
        my $end = $fields[4]; #end location of designation
        if ($ident eq "mRNA"){
        	my $exon_num = 1;
        } elsif ($ident eq "exon"){ 
        	my $trans_id = substr($info, 7);
        	my $exon_name;
        	if (defined $transcripts{$trans_id}){ #makes hashes of the exons for each transcript with there start and end
        		if (defined $exons_start{$trans_id}){
        			$exon_num++;
        			$exon_name = "exon_$exon_num";
        			$exons_start{$trans_id}{$exon_name} = $start;
        			$exons_end{$trans_id}{$exon_name} = $end;
        		} else {
        			$exon_name = "exon_$exon_num";
        			$exons_start{$trans_id}{$exon_name} = $start;
        			$exons_end{$trans_id}{$exon_name} = $end;
        		}
        	}
        } elsif ($ident eq "CDS"){
        	my $trans_id = substr($info, 7);
        	if (defined $transcripts{$trans_id}){
	            if (defined $mRNAs_start{$trans_id}){
        	    } else {
            	    $mRNAs_start{$trans_id} = $start; #hash of transcript start
            	}
            	if (defined $mRNAs_start{$trans_id}){
            		$mRNAs_end{$trans_id} = $end; #hash of transcript end
        	    } else {
            	}
            }
        }
    }
}
close($transdecoder);

#makes table of genes with exons, there start stop and type (utr or coding)#
my $gene_isoforms = 0;
my %UTR_exon;
my %UTR_five;
my %UTR_three;
print $new_list "Transcript\tExon_Type\tExon_Start\tExon_End\tCDS_Start\tCDS_End\n";
foreach my $trans (sort keys %transcripts){
	foreach my $three_exon (sort keys %{$exons_start{$trans}}){
		if ($exons_start{$trans}{$three_exon} > $mRNAs_end{$trans}){
			$UTR_exon{$trans}++;
			$UTR_three{$trans} = 1;
			print $new_list "$trans\t3Prime\t$exons_start{$trans}{$three_exon}\t$exons_end{$trans}{$three_exon}\t$mRNAs_start{$trans}\t$mRNAs_end{$trans}\n";
		}
	}
	foreach my $five_exon (sort keys %{$exons_end{$trans}}){
		if ($exons_end{$trans}{$five_exon} < $mRNAs_start{$trans}){
			$UTR_exon{$trans}++;
			$UTR_five{$trans} = 1;
			print $new_list "$trans\t5Prime\t$exons_start{$trans}{$five_exon}\t$exons_end{$trans}{$five_exon}\t$mRNAs_start{$trans}\t$mRNAs_end{$trans}\n";
		}
	}
}
close($new_list);

#makes hash of genes with isoforms
my %isoforms;
while (my $iso_list = <$genelist>){
	chomp $iso_list;
	my @fields = split /\t/, $iso_list;
	my $name = $fields[1];
	$isoforms{$name} = 1;
}
close($genelist);

#for genes with isoforms, is there a change in the exons creating that isoform
my $trans_exons;
my $greater_two;
my $both;
my $five;
my $isocount;
foreach my $mrna (sort keys %UTR_exon){
	$trans_exons++;
	if ($UTR_exon{$mrna} >= 2){
		$greater_two++;
	}
	if (defined $UTR_five{$mrna}){
		if (defined $UTR_three{$mrna}){
			$both++;
		} else {
			$five++;
		}
	}
	if (defined $isoforms{$mrna}){
		$isocount++;
	}
}
my $three = int($trans_exons - ($both + $five));
print "Transcripts with UTR exons: $trans_exons\nTranscripts with more than 2 exons: $greater_two\nFive : $five\nThree: $three\nBoth: $both\nIsoforms = $isocount\n";


