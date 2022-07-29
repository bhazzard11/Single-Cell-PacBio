#!/usr/bin/perl -w                                                                                                                                                                                                    
use strict;
use warnings;

###Creates fasta file for mapping###
##Removes adapters on both ends, finds barcode/UMI and trims PolyA tail(trims from revcomplement)##

#input 1: ccs bam file#
open (my $BAM, "-|", "/path/to/samtools", "view", "/path/to/pacbio/reads.ccs.bam");

#output: name and path for output fasta# 
open (my $trim, ">", "/path/to/output.fasta");

#input 2: List of cell barcodes to consider, txt in this case a count table with the accepted barcodes#
open (my $cells, "<", "/path/to/accepted/barcodes/from/illumina.txt");

#Make hash of barcodes to check against#
my $counter = 0;
my %seenbarcode;
my $matches = 0;
my $recovered = 0;
my $trash = <$cells>;
my %illuminacells;
while (my $line1 = <$cells>){ #cell count list
    chomp $line1;
    my @fields = split / /, $line1; #splits lines into an array 
    my $barcode = $fields[0];
    if ($fields[1] >= 250){ #from count table, can change number to increase min reads to consider
	$illuminacells{$barcode} = 1;
	$seenbarcode{$barcode} = 1;
    }
}

#finds pacbio adapters and removes them (CCGATCT), reverse complements sequence if doesn't find adapter#
#looks for barcode which should be within 28bp of adapter and makes that and UMI read name, allows for one nucleotide mismatch#
#writes to a fasta file#
#keeps track of where sequences are lost, can't find adapter, can't find barcode#
my $has_firstadapter = 0;
my $has_barcode = 0;
my $has_secondadapter = 0;
while (<$BAM>){
    chomp;
    my @fields = split /\t/; #splits lines into an array                                                                                
    my $ID = $fields[0]; #read ID
    my $seq = $fields[9]; #Sequence
    my $barcodeUMI = 0;
    my $position;
    my $seqname;
    my $sequence;
    my $flag = 0;
	if ($seq =~ "CCGATCT"){ #find forward addapter
	    $has_firstadapter++;
	    $position = $+[0]; #position in sequence where adapter is found
	    $barcodeUMI = substr $seq, $+[0], 32; #store barcode and UMI that follows adapter
	    #	    $counter++;
	    $flag++;
	} else { #if don't find the adapter reverse transcribe the sequence and look again
	    my $revcomp = reverse($seq);
	    $revcomp =~ tr/ATCG/TAGC/;
	    $seq = $revcomp;
	    if ($revcomp =~ "CCGATCT"){
		$has_firstadapter++;
		$position = $+[0];
		$barcodeUMI = substr $revcomp, $+[0], 32;
		#		$counter++;
		$flag++;
	    } else {}
	}
    if ($flag == 1){
	my %rejects;
#	print "$ID\n";
#    for (my $i=0; $i<5; $i++){ #Look for barcode
	if ((length($barcodeUMI)) > 28){
#	    print "$barcodeUMI\n";
	    my $subseq = substr $barcodeUMI, 0, 16;#$i, 16;
	    if (defined $seenbarcode{$subseq}){
		$has_barcode++;
		my $UMI = substr $barcodeUMI, 16, 12;#($i + 16), 12;
		    my @name =  ($ID, $subseq, $UMI); #For unique sequences create name from ID, barcode, and UMI
#                    $barcodes{$subseq}++;                                                                                                                                                       
		$seqname = join("_", @name);
		    $sequence = substr $seq, ($position + 28);#($position + $i + 28);
#		print "$seqname\t$sequence\n";
		    if ($sequence =~ "CCCATGT"){ #"AAGCAGT"){ #find reverse adapter and remove it
      			$has_secondadapter++;
			$sequence = substr $sequence, 0, $-[0];
		    }
#			my $revseqtrim = reverse ($sequence);#($seqtrim); #reverse transcribe sequence so all polyA
#			$revseqtrim =~ tr/ATCG/TAGC/; 
		#			print $fasta ">$seqname\n$revseqtrim\n"; #print to file
#		print "$seqname\t$sequence\n";
		$sequence =~ s/^[^TTTTTTTT*]*TTTTTTTT*//; #trim sequence at 8T's
		if ((length($sequence)) > 10){
#		    print "$seqname\t$sequence\n";
		    my $revseqtrim = reverse ($sequence);#($seqtrim); #reverse transcribe sequence so all polyA  
		    $revseqtrim =~ tr/ATCG/TAGC/;
			    print $trim ">$seqname\n$revseqtrim\n";
			    $matches++; #count accepted sequences
			    #$i=6;
#			}
			}
#		    }                                                                                                                                                                   
#		$i=6;
	    } else {
		$rejects{$subseq} = 1;#$i;
	    }
	} else {
	}
#    }
    my $numrejects = keys %rejects;
#    if ($numrejects >= 5){
	foreach my $querybarcode (sort keys %rejects){
#	    my $value = $rejects{$querybarcode};
#	    my $UMI = substr $barcodeUMI, ($value + 16), 12;
#	    my @name =  ($ID, $querybarcode, $UMI);
#	    $seqname = join("_", @name);
#	    $sequence = substr $seq, ($position + $value + 28);
	    foreach my $testbarcode (sort keys %illuminacells){
		my $progress = 0;
		my $lettermatch = 0;
		my $badmatch = 0;
		if ($badmatch < 2){
		    foreach my $test (split //, $testbarcode){
		    my $query = substr $querybarcode, $progress, 1;
		    if ($test eq $query){
			$lettermatch++;
			$progress++;
		    } else {
			$progress++;
			$badmatch++;
		    }
		}
		    if ($lettermatch == 15){
			my $value = $rejects{$querybarcode};
			my $UMI = substr $barcodeUMI, ($value + 16), 12;
			my @name =  ($ID, $querybarcode, $UMI);
			$seqname = join("_", @name);
			$sequence = substr $seq, ($position + $value + 28);
			$has_barcode++;
			if ($sequence =~ "CCCATGT"){ #"AAGCAGT"){ #find reverse adapter and remove it                                                                     
			    $has_secondadapter++;
			    $sequence = substr $sequence, 0, $-[0];
			}
		#	my $revseqtrim = reverse ($sequence);#($seqtrim); #reverse transcribe sequence so all polyA                                                                        
		#	$revseqtrim =~ tr/ATCG/TAGC/;
			#	$revseqtrim =~ s/AAAAAAAA.*//; #trim sequence at 8A's
#			print "$seqname\t$sequence\n";
			$sequence =~ s/^[^TTTTTTTT*]*TTTTTTTT*//; #trim sequence at 8T's  
			if ((length($sequence)) > 10){
#			    print "$seqname\t$sequence\n";
			    my $revseqtrim = reverse ($sequence);#($seqtrim); #reverse transcribe sequence so all polyA                                                          
			    $revseqtrim =~ tr/ATCG/TAGC/;
			    print $trim ">$seqname\n$revseqtrim\n";
			    $recovered++;
			    $numrejects = 6;
			    $seenbarcode{$querybarcode} = 1;
			}
		    }
		}
	 #   }
	}
    }
	}
}

print "number with one adapter = $has_firstadapter\n number with both adapters= $has_secondadapter\nnumber with barcode = $has_barcode\nnumber recovered = $recovered\n";
print "accepted sequences = $matches\n";
close ($BAM);
close ($trim);
close ($cells);
