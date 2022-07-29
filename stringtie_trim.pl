#!/usr/bin/perl -w                                                                                                                                                                 

use strict;
use warnings;

####Takes stringtie gtf input and remove low abundance transcript predictions####

#input 1: concatenated stringtie output gtf#
open (my $abundance1,"<","/path/to/stringtie/output/merged.gtf");
#output 1: output.gtf#
open (my $output1, ">","/path/to/output_merged_cutoff.gtf");

my %isoforms;
my %duplications;
my $C01counts;
my $totaltrans;
my $total_cov;
while (my $line1 = <$abundance1>){ #goes through each transcript prediction and finds stringties calculated abundance. If >= 10 keeps, else throws out#
    chomp $line1;
    my @fields = split /\t/, $line1; #splits lines into an array                                                                                                                   
    my $ident = $fields[2];
    my $info = $fields[8];
    my $chr = $fields[0];
    my $start = $fields[3];
    my $end = $fields[4];
    my $strand = $fields[5];
    my $window = int($end/500);
    my @location = ($chr, $window, $strand);
    my $loc = join ("_", @location);
    if ($ident eq "transcript"){
	$totaltrans++;
        my @collapse_info = split /;/, $info;
        my $trans_id = $collapse_info[1];
        my $cov = $collapse_info[2];
#	my $FPKM = $collapse_info[3];
	my $altcov = $collapse_info[4];
        my ($temp1, $ID) = $trans_id =~ /([^\s"]+)/g;                                                                                                                     
        my ($temp2, $coverage) = $cov =~ /([^\s"]+)/g;
#	my ($temp3, $FPKM_val) = $FPKM =~ /([^\s"]+)/g;
        my ($temp4, $altcoverage) = $altcov =~ /([^\s"]+)/g;
        $total_cov = ($total_cov + $coverage);
	if ($cov =~ /cov/){
	    if ($coverage >= 10){
		$duplications{$loc} = 1;
		$isoforms{$ID} = 1;
		$C01counts++;
		print $output1 "$line1\n";
	    } else {}
	} else {
	    if ($altcoverage >= 10){
                $duplications{$loc} = 1;
                $isoforms{$ID} = 1;
                $C01counts++;
                print $output1 "$line1\n";
            } else {}
	}
    } elsif ($ident eq "exon"){
	my @collapse_info = split /;/, $info;
	my $trans_id = $collapse_info[1];
	my ($temp1, $ID) = $trans_id =~ /([^\s"]+)/g; 
	if (defined $isoforms{$ID}){
	    print $output1 "$line1\n";
	} else {}
    } else {}
}

close ($abundance1);
close ($output1);

print "C01 transcripts = $C01counts, total transcripts = $totaltrans, sample 1 complete\nTotal Coverage = $total_cov\n";
