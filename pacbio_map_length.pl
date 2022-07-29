#!/usr/bin/perl -w                                                                                                                                                        \                       

use strict;
use warnings;

####Minimap2 is a very liberal mapper and maintains reads that map lowly to the genome####
####Removes reads with less than 50% mapping to reference####

#input 1: fasta of used for mapping#
open (my $fasta, "<", "/path/to/ccs/reads.fasta"); #fasta
#input 2: bed of mapped reads from Pacbio_analysebam.pl output#
open (my $bed, "<", "/path/to/mapping/input_rmdup.bed"); #bed
#input 3: bam output of Pacbio_analysebam.pl#
open (my $BAM, "-|", "/path/to/samtools", "view", "/path/to/mapping/input_rmdup.bam"); #bam
#output: sam for output#
open (my $out_sam, ">", "/path/to/output.sam"); #out sam

my %fasta_len;
my $name;
my $counter;
while(my $line1 = <$fasta>){ #makes list of sequences from original fasta input and there lengths
    chomp $line1;
    if ($line1 =~ ">"){
	$counter =1;
	$name = substr($line1, 1);
    } elsif ($counter == 1) {
	my $len = length($line1);
	$fasta_len{$name} = $len;
	$counter++;
    }
}
close($fasta);

my %per_len;
my %seen_reads;
my %good_reads;
while (my $line2 = <$bed>){ #compares lengths of reads when fasta and reads after mapping, if 50% or more has been masked, remove. Stores read names to be kept
    chomp $line2;
    my @fields = split /\t/, $line2;
    my $start = $fields[1];
    my $end = $fields[2];
    my $read = $fields[3];
    my $len = ($end - $start);
    if (defined $seen_reads{$read}){
	my $real_len = $fasta_len{$read};
	my $old_len = $seen_reads{$read};
	my $new_len = ($old_len + $len);
	my $per = int(($new_len/$real_len)*100);
	$seen_reads{$read} = $new_len;
        $per_len{$read} = $per;
	if ($per >= 50){
	    $good_reads{$read} = 1;
	}
    } else {
	$seen_reads{$read} = $len;
	my $real_len = $fasta_len{$read};
	my $per = (($len/$real_len)*100);
	$per_len{$read} = $per;
	if ($per >= 50){
            $good_reads{$read} = 1;
	}
    }
}
close($bed);

while (my $line_bam = <$BAM>){ #removes sequences that did not pass previous test, prints accepted ones to new sam file
    chomp $line_bam;
    my @fields = split /\t/, $line_bam; #splits lines into an array                                                            
    my $ID = $fields[0]; #read ID 
    if (defined $good_reads{$ID}){
	print $out_sam "$line_bam\n";
    }
}
close($BAM);
close($out_sam);

