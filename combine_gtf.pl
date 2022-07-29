#!/usr/bin/perl -w                                                                                                                                                                                                                                                               

use strict;
use warnings;

####Takes gffcompare output and removes transcripts not redundant in both sample files. Gffcompare outputs a "=" or "j" for transcripts where they exist in both files compared###

#input: gffcompare output gtf#
open (my $abundance1,"<","/path/to/gffcompare/input.annotated.gtf");
#output: output.gft#
open (my $output, ">", "/path/to/output.alltranscripts.gtf");

my %geneid;
while (my $line1 = <$abundance1>){
    chomp $line1;
    my @fields = split /\t/, $line1; #splits lines into an array                                                                                                           
    my $ident = $fields[2]; #transcript or exon                                                                                                                                   
    my $info = $fields[8]; #contains transcript_id and gene_id infor                                                                                                     
    my @collapse_info = split /;/, $info; 
    if ($ident eq "transcript"){
	my $size = @collapse_info;
        my $trans_id = $collapse_info[0];
	if ($size > 6){
	    my $class_code = $collapse_info[5];
	    my $altclass_code = $collapse_info[6];
	    my ($temp1, $transcript_id) = $trans_id =~ /([^\s"]+)/g;
	    my ($temp2, $class) = $class_code =~ /([^\s"]+)/g;
	    my ($temp2alt, $classalt) = $altclass_code =~ /([^\s"]+)/g;
	    if ($class eq "="){
		print $output "$line1\n";
		$geneid{$transcript_id} = 1;
	    } elsif ($classalt eq "="){
		print $output "$line1\n";
		$geneid{$transcript_id} = 1;
            } elsif ($classalt eq "j"){                                                                                                                                                     
                print $output "$line1\n";                           
                $geneid{$transcript_id} = 1; 
            } elsif ($class eq "j"){                                                                                                                                      
                print $output "$line1\n";                                                                                                  
                $geneid{$transcript_id} = 1; 
	    } else {}
	} else {
	    my $altclass_code = $collapse_info[3];
	    my ($temp1, $transcript_id) = $trans_id =~ /([^\s"]+)/g;
	    my ($temp2alt, $classalt) = $altclass_code =~ /([^\s"]+)/g;
            if ($classalt eq "="){
                print $output "$line1\n";
                $geneid{$transcript_id} = 1;
	    }
	}
    } elsif ($ident eq "exon"){
	my $trans_id = $collapse_info[0];
	my ($temp1, $transcript_id) = $trans_id =~ /([^\s"]+)/g;
        if (defined $geneid{$transcript_id}){
            print $output "$line1\n";
	} else {}
	}
}					  

close ($abundance1);
close ($output);						 
