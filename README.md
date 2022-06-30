# Single-Cell-PacBio
Single cell PacBio analysis scripts

perl scripts for analyses of data from PacBio sequencing of single cell libraries

Pipeline:
1. Processing and Analysis of Illumina Reads:
  Identify barcode and UMI information in sequencing fasta files:
    Creates combined fasta for mapping for hisat2
  Make counts table from bam file: analysebamrmunmapped.pl

2. Processing of PacBio Reads:
  Find barcode and UMI from ccs file: pacbio_bamtofasta.pl
    Creates fasta for mapping with minimap2:
      minimap2 -a --secondary=no -x cdna -k14 -G 50000 Genome.fasta PacBio_out.fasta
  Removes duplicates from bam file: pacbio_analysebam.pl
  Removes lowly mapped reads: pacbio_map_length.pl
  
 3. Stringtie2 for transcript prediction:
   Split bame file into forward and reverse strands: 
     samtools view -h -F 16 in.bam > out.sam
   Add XStag to sam file
   Split bam file into mono vs multi exon reads: 
     samtools view -h in.forwardxs.bam | awk '{if($0 ~ /^@/ || $6 !~ /N/) {print $0}}' | samtools view -Sb - > out.forward.monoexon.bam
     samtools view -h in.forwardxs.bam | awk '{if($0 ~ /^@/ || $6 ~ /N/) {print $0}}' | samtools view -Sb - > out.forward.introns.bam
   Run stringtie:
    stringtie in.forward.introns.bam -o out_introns_f.gtf -R -L -l name
    (repeat for all 4 files)
   Remove low support transcript predictions (cov < 10): stringtie_trim.pl
   Combined filtered transcripts with gffcompare to collapse redundant: gffcompare -r in_cutoff.gtf -T -R -o out_combined out_cutoff.gtf
   Filter merged file to only include trancripts in both samples: combine_gtf.pl
   
 4. Transdecoder for coding transcript prediction:
    Run transdecoder on stingtie final output: 
      transdecoder-5.3.0/util/gtf_genome_to_cdna_fasta.pl in.alltranscripts.gtf Genome.fasta > out.alltranscripts.fasta
      transdecoder-5.3.0/TransDecoder.LongOrfs -t in.alltranscripts.fasta
      transdecoder-5.3.0/TransDecoder.Predict -t in.alltranscripts.fasta
      
 5. Analysis of predicted transcripts


  
