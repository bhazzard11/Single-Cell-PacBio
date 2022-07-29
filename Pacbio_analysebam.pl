#!/usr/bin/perl -w                                                                                                                                                                                                              

use strict;
use warnings;

###Removes PCR duplicates from sorted BAM file. Finds neccessary droplet and UMI information, prints to a file###                                                                                                                  

#input 1: bam file of mapped pacbio reads#
open (my $BAM, "-|", "/path/to/samtools", "view", "/path/to/bam_input.bam"); 
#output: text file of counts#
open (my $output, ">", "/path/to/new/barcode/list/out.txt");
#input 2: List of accepted barcodes#
open (my $file2, "<", "/path/to/illumina/barcode/list/in.txt");
#output 2: reduced sam of accepted reads#
open (my $out_reads, ">", "/path/to/new/rmdup_output.sam");

my %seenbarcode;
my %results;
my %barcodes;
my %unique;
my %acceptedbarcode;
my $read_count = 0;
my $dup_reads = 0;
while (my $line2 = <$file2>){
    chomp $line2;
    my $barcode2 = substr $line2, 0, 16;
    $acceptedbarcode{$barcode2} = 1;
}

while (my $line = <$BAM>){
    chomp $line;
    my @fields = split /\t/, $line; #splits lines into an array                                                                                                                                                                       
    my $name = $fields[0];
    my $barcode = substr $name, -29, 16; #10x droplet barcode                                                                                                                                                                 
    my $UMI = substr $name, -12, 12; #UMI 12mer                                                                                                                                                                               
    my $chromosome = $fields[2]; #chromosome read mapped to                                                                                                                                                                   
    my $location = $fields[3]; #location read mapped to on chromosome                                                                                                                                                         
    my $flag = $fields[1]; #stand information of read map                                                                                                                                                                    
    my $strand;
    if ($flag == 0){
        $strand = "+";
	$read_count++;
    } elsif ($flag == 16){
        $strand = "-";
	$read_count++;
    } else {
        $flag = "unknown";
	next;
    }
    my $window =int ($location/500);
    my @newlist = ($barcode, $UMI, $strand, $chromosome, $location, "\n"); #creates new array with just the fields I want                                                                                                     
    my @list2 = ($chromosome, $strand, $window);
    my $identity = join(",", @newlist); #creates a flattened string with my fields seperated by a comma                                                                                                                       
    my $G = join(",", @list2);
	 if (defined $acceptedbarcode{$barcode}){ #Checks if barcode is present in Illumina data
		if (defined $seenbarcode{$identity}){} #Checks if PCR duplicate with same location, barcode, and UMI
		else {
	    $dup_reads++;
	 	   print $out_reads "$line\n"; #prints accepted sequences to new sam file
	  	  $seenbarcode{$identity} = 1;
	 	   $results{$barcode}{$G}++;
	 	   if (defined $unique{$G}){}
	 	   else {
			$unique{$G} = 1;
	 	   }
	 	   if ($barcodes{$barcode}) {
			$barcodes{$barcode}++;
	 	   } else {
			$barcodes{$barcode} = 1;
	 	   }
		}
    } else {}
}
my @header = sort keys %barcodes;# $barcodes{$barcode} >=5000;                                                                                                                                                             
print "Finding Barcodes\n";
my @genes = sort keys %unique;
                                                                                                                                                                                                                        
print "Printing to file\n";
for my $key (sort keys %barcodes){
    print $output "$key\t$barcodes{$key}\n"; #creates new intermediate file with key barcode pairs                                                                                                                                                                    
}

print "Total Mapped Reads = $read_count\nTotal after Duplicate Removal = $dup_reads\n"; #prints to screen what is lost in this QC
close ($BAM);
close ($output);
close ($file2);
