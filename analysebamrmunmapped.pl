#!/usr/bin/perl -w                                                                                                                                 
                                                                                                                                                    
use strict;
use warnings;

#Removes PCR duplicates from sorted BAM file. Finds neccessary droplet and UMI information, prints count table to a file#    
                                  
#output 1: count table, large matrix                                                                                                                                                
open (my $output, ">", "/path/to/count/table/output_analysis.txt");
#input 1: mapped reads bam file
open (my $BAM, "-|", "/path/to/samtools", "view", "/path/to/hisat/output/mapped.bam");
#input 2: list of 10X barcodes
open (my $file1, "-|", "gzip", "-dc", "/path/to/barcode/list.txt.gz");
#input 3: fastq with readnames in case hisat doesn't propagate
open (my $readname, "-|", "gzip", "-dc", "/path/to/reads/merged.fastq.gz");

my %acceptedbarcode;
while (my $line = <$file1>){ #creates an array of barcodes                                                                                               
     chomp $line;
#     my @fields = split /\t/, $line;
#     my $barcode = $fields[0];
     my $barcode = $line;
     $acceptedbarcode{$barcode} = 1;
}

##for if for some reason hisat didn't propagate the actual barcode infor
my %read_names;
#my $counter = 0;
#my $sequencer = 0;
while (my $readline = <$readname>){
    chomp $readline;
	    my @fields = split ' ', $readline;
	    my $read = substr $fields[0], 1;
	    my $name = $fields[1];
	    $read_names{$read} = $name;
}

#parse bam file for counts. Uses 500 basepair windows based on bam chr and location#
my %seenbarcode;
my %results;
my %barcodes;
my %unique;
while (<$BAM>){
    chomp;
    my @fields = split /\t/; #splits lines into an array                                                                                             
    my $name = $fields[0];
    my $barcode;
    my $UMI;
    if (defined $read_names{$name}){
	my $info = $read_names{$name};
	my @cell_info = split /_/, $info;
#        my $barcode = substr $name, -29, 16;
	$barcode = $cell_info[2]; #my $barcode = substr $name, -27, 16; #10x droplet barcode                                                     
#	my $UMI = substr $name, -12, 12;
	$UMI = $cell_info[3];#my $UMI = substr $name, -10, 10; #UMI 10mer                              
	my $chromosome = $fields[2]; #chromosome read mapped to
	my $location = $fields[3]; #location read mapped to on chromosome                                                 
	my $flag = $fields[1]; #stand information of read map                                                                                           
	my $strand;
	if ($flag == 0){ #sam flags for forward and reverse reads that mapped to genome
	    $strand = "+";
	} elsif ($flag == 16){
	    $strand = "-";
	} else {
	}
	my $window =int ($location/500); #creates 500bp window
	if (defined $strand){
#	} else {
	    my @newlist = ($barcode, $UMI, $strand, $chromosome, $location, "\n"); #creates new array with just the fields I want                          
	    my @list2 = ($chromosome, $window, $strand);
	    my $identity = join("_", @newlist); #creates a flattened string with my fields seperated by a comma                                                                                        
	    my $G = join("_", @list2);
	    if (defined $acceptedbarcode{$barcode}) {
		if (defined $seenbarcode{$identity}){}
		else {
		    $seenbarcode{$identity} = 1;
		    $results{$barcode}{$G}++;
		    if ($unique{$G}){}
		    else {
			$unique{$G} = 1;
		    }
		    if (defined $barcodes{$barcode}) {
			$barcodes{$barcode}++;
		    } else {
			$barcodes{$barcode} = 1;
		    }
		}
	    } else {}
	}
    } else {}
}
print "Finding Barcodes\n";
my @genes = sort keys %unique;
my %parasites;
for my $count (sort keys %barcodes){
    if ($barcodes{$count} >= 5000){
        $parasites{$count} = 1;
    } else {}
}
print "Printing to file\n";
print $output join ("\t", 'Parasite', @genes), "\n";
for my $key (sort keys %barcodes){
    my @data;
    if (defined $parasites{$key}){
        for my $allgenes (sort keys %unique){
            if (defined $results{$key}{$allgenes}){
                push @data, ($results{$key}{$allgenes});
            } else {
                push @data, "0";
            }
        }
        print $output join("\t", $key, @data), "\n";
    } else {}
}

close ($BAM);
close ($output);
close ($file1);
