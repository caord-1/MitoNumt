# NUMT_get
This is a pepline for easy to get the true NUMT and MT sequence from the HIFI data.

Get the ref_numt sequence

   need the mt and genome ref

    minimap2 -x asm5 ./ref_mt.fa ./ref_genomic.fna > *.paf

   Then, selected the ref_numt sequence from to .paf

   <img width="1667" height="197" alt="image" src="https://github.com/user-attachments/assets/bb012dc2-2b8c-48c7-a393-0f418c9d0ff8" />

1. Extend 2KB from the starting position on each side.（Suitable for finding complete NUMT sequences）
   
such as:
    
	  samtools faidx ./ref_genomic.fna CM101310.1:68556902-685561280 > CM101310.1.68556902-68561280.fa

（1） get the NUMT HIFI reads

    minimap2 -t 16 --secondary=no -ax map-hifi ./CM101310.1.68556902-68561280.fa ./HIFI_subreads.fastq.gz | samtools view -@ 16 -b -F4 -F 0x800 -o reads.HiFiMapped.bam
   
    samtools sort -@ 16 -o reads.HiFiMapped.sorted.bam reads.HiFiMapped.bam
   
    samtools index -@ 16 reads.HiFiMapped.sorted.bam

（2） Next, IGV is used to visualize the bam file.


   <img width="765" height="314" alt="image" src="https://github.com/user-attachments/assets/e4aea820-d1d5-42f1-96ef-a4a5b969bd47" />


  Prioritize the HIFI reads that are aligned to the extended 2kb region. Secondly, choose the reads with a higher alignment rate.
   
 （3）performed the multiple sequence alignment with the selected HIFI reads and ref mt to get the numt sequence.


2. do not extend 2KB from the starting position on each side.（Suitable for finding quickly short NUMT sequences）

 such as: 
 
    samtools faidx ./ref_genomic.fna CM101310.1:68558902-68559280 > CM101310.1.68558902-68559280.fa

 （1） get the NUMT HIFI reads

    minimap2 -t 16 --secondary=no -ax map-hifi ./CM101310.1.68558902-68559280.fa ./HIFI_subreads.fastq.gz | samtools view -@ 16 -b -F4 -F 0x800 -o reads.HiFiMapped.bam
   
	 samtools sort -@ 16 -o reads.HiFiMapped.sorted.bam reads.HiFiMapped.bam
   
	 samtools index -@ 16 reads.HiFiMapped.sorted.bam

   <img width="1558" height="174" alt="image" src="https://github.com/user-attachments/assets/9fe82883-915d-4a06-ab37-96c802e7f160" />

   
	 python3 ./script/extract_high_identity_reads.py -b reads.HiFiMapped.sorted.bam -o CM101310.1.68558902-68559280 -d pos

   this step will generate two output files:CM101310.1.68558902-68559280.numt.high98.full.fa and CM101310.1.68558902-68559280.numt.high98.trim.fa

   (2) get consensus numt sequence
   
	 cat ref_mt.fa CM101310.1.68558902-68559280.numt.high98.trim.fa > CM101310.1.68558902-68559280.numt.high98.trim.mt.fa

    mafft --auto  CM101310.1.68558902-68559280.numt.high98.trim.mt.fa  CM101310.1.68558902-68559280.numt.high98.trim-mt.mafft.fa
   
    python3 ./script/numt_consensus_from_msa.py CM101310.1.68558902-68559280.numt.high98.trim-mt.mafft.fa ref_mt CM101310.1.68558902-68559280.numt


finally you will get the short consensus numt squence. CM101310.1.68558902-68559280.numt.fa


Get the mt sequence
  
  
  Just replace the reference sequence with "mt", and the remaining steps remain the same numt.

   

   

