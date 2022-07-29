#!/usr/bin/perl -w                                                                                                                                                                                                                                                                                                                        

use strict;
use warnings;

####count genes sequences that alingned by blat to a transcript that is from a gene with an isoform####

#input: blat output
open (my $blat, "<", "/path/to/blat/output_blat_out.psl");
#input 2: tab file with pseudotime values for each cell
open (my $pseudotime, "<", "/path/to/pseudotime/cell/values/pseudotime.txt");
#output 1: count table of expression for each transcript in each cell
open (my $output, ">", "/path/to/expression/counts/from/blat/out_absolute_expression.txt");
#output 2: pseudotime, metadata table
open (my $output2, ">", "/path/to/metadata/output_info.txt");

#looks through blat output and determines the best match for each predicted transcript#
my %ref_id;
my %ref_len_cov;
my %trans_match;
my %q_len_cov;
my %q_dif;
my %new_ref;
my $read_count = 0;
my $line_count = 0;
while (my $line = <$blat>){
    chomp $line;
    $line_count++;
 	if ($line_count >= 12){
 		my @fields = split /\t/, $line;
    	my $strand = $fields[8];
    	my $name = $fields[9];
    	my $gene = $fields[13];
		my $align_len = $fields[0];
   		my $ref_len = $fields[14];
   		my $trans_len = $fields[10];
    	my $per_len = ($align_len / $trans_len)*100;
    	if ($per_len > 100){
			$per_len = ($trans_len / $align_len)*100;
    	}
    	my $len_dif = abs($ref_len - $align_len);
    	my @name_split = split /_/, $name;
    	my $cell = $name_split[3];
   		my $trans = $name_split[4];
    	if (defined $trans_match{$cell}{$trans}){
			my $len_comp = abs($len_dif - $q_dif{$cell}{$trans});
			if ($strand eq "+"){
	    		if ($per_len >= $q_len_cov{$cell}{$trans}){
					if ($len_dif < $q_dif{$cell}{$trans}){
		   				$q_len_cov{$cell}{$trans} = $per_len;
		    			$trans_match{$cell}{$trans} = $gene;
		    			$q_dif{$cell}{$trans} = $len_dif;
		    			$new_ref{$cell}{$trans} = $ref_len;
					}
	    		}
			}
    	} else {
			if ($per_len >= 90) {
	    		if ($strand eq "+"){
					$read_count++;
					$q_len_cov{$cell}{$trans} = $per_len;
					$trans_match{$cell}{$trans} = $gene;
					$q_dif{$cell}{$trans} = $len_dif;
					$new_ref{$cell}{$trans} = $ref_len;
	    		}
			}
		}
	} else {}
}
close($blat);
print "$read_count\n";

##makes hash for cell with pseudotime values
my %pseudo_num;
my $time = 0;
while (my $list = <$pseudotime>){
    chomp $list;
    my @fields = split /\t/, $list;
    my $name = $fields[0];
    $time++;
    my @split_name = split /_/, $name;
    my $sample = $split_name[0];
    my $cell = $split_name[1];
    if ($sample eq "sp5"){
	$pseudo_num{$cell} = $time;
    }
}
close($pseudotime);

#input 3: coding transcripts file
open (my $transdecoder, "<", "/local/projects-t3/SerreDLab-3/hazzardb/SingleCellData/Annotation/transdecoder_genome_v54/Bloodstage_combined.coding.gtf");
#input 4: name match file for 2 samples
open (my $alt_names, "<", "/local/projects-t3/SerreDLab-3/hazzardb/SingleCellData/Annotation/stringtie_genome_v54/Bloodstage_combined.annotated.gtf");

#makes hash of transcripts with the names for sample 1
my %transcript_list;
my %alt_list;
while (my $new_ref_line = <$transdecoder>){
    chomp $new_ref_line;
	my @fields = split /\t/, $new_ref_line;
	if (scalar(@fields) > 3){
   		my $ident = $fields[2];
   		my $start = $fields[3];
    	my $end = $fields[4];
    	my $len = int($end-$start);
    	my $info = $fields[8];
    	if ($ident eq "mRNA"){
			my @collapse_info = split /;/, $info;
			my $pre_trans_name = substr($collapse_info[0], 3);
			my $trans_name = substr($pre_trans_name, 0,-3);
       		$transcript_list{$trans_name} = 1;
    	}
	}
}
close($transdecoder);

#makes hash of transcripts with the names for sample 2
while (my $alt = <$alt_names>){
    chomp $alt;
	my @fields = split /\t/, $alt;
	my $ident = $fields[2];
	my $info = $fields[8];
	if ($ident eq "transcript"){
    	my @collapse_info = split /;/, $info;
    	my $trans_id = $collapse_info[0];
   		my $trans_alt = $collapse_info[4];
        my ($temp1, $trans_prop) = $trans_id =~ /([^\s"]+)/g;                                                                                                                                                                                             
        my ($temp2, $trans_old) = $trans_alt =~ /([^\s"]+)/g;
    	if (defined $transcript_list{$trans_prop}){
        	$alt_list{$trans_old} = $trans_prop;
    	}
	}
}
close($alt_names);

#transcripts by cell list
my %cell_trans;
foreach my $cell (sort keys %trans_match){
    foreach my $transcript (sort keys %{$trans_match{$cell}}){
		my $match = $trans_match{$cell}{$transcript};
		if (defined $transcript_list{$match}){
	   		$cell_trans{$cell}{$match}++;
		} elsif (defined $alt_list{$match}){
	    	my $switch = $alt_list{$match};
	    	$cell_trans{$cell}{$switch}++;
		}
    }
}


print $output "Cell";
print $output2 "Cell\tpseudotime1\tpseudotimeF\tcurve1\tcurve2\n";
foreach my $header (sort keys %transcript_list){
    print $output "\t$header";
}

foreach my $cell_name (sort keys %cell_trans){
    my @cell_data;
	if (defined $pseudo_num{$cell_name}){
   		print $output2 "$cell_name\t$pseudo_num{$cell_name}\t0\t1\t0\n";
    	foreach my $trans_name (sort keys %transcript_list){
			if (defined $cell_trans{$cell_name}{$trans_name}){
	    		push @cell_data, $cell_trans{$cell_name}{$trans_name};
			} else {
	    		push @cell_data, 0;
			}
    	}
    	my $cell_out = join("\t", @cell_data);
    	print $output "\n$pseudo_num{$cell_name}\t$cell_out";
	}
}
close($output);      
