# NUMT_get
Overview

This pipeline provides a practical and reproducible workflow to accurately get true NUMT (nuclear mitochondrial DNA segments) and mitochondrial (mtDNA) sequences from PacBio HiFi sequencing data.
requirements

**python3 minimap2 samtools mafft igv**

**Step 1. Identification of candidate NUMT regions in the nuclear genome**
Input requirements

Mitochondrial reference genome (ref_mt.fa)

Nuclear genome reference (ref_genomic.fna)

Candidate NUMT loci are identified by aligning the mitochondrial genome to the nuclear genome using minimap2:

    minimap2 -x asm5 ./ref_mt.fa ./ref_genomic.fna > *.paf

The resulting PAF file is then inspected to extract nuclear regions showing significant homology to the mitochondrial genome, which are considered putative NUMT loci.

   <img width="1667" height="197" alt="image" src="https://github.com/user-attachments/assets/bb012dc2-2b8c-48c7-a393-0f418c9d0ff8" />

**Step 2. Extraction of reference NUMT regions**

*Strategy 1*: Extended flanking regions (recommended for detecting complete or near-complete NUMTs)

To get potentially full-length NUMTs, the candidate NUMT region is extended by 2 kb upstream and downstream from the alignment coordinates:
   
	  samtools faidx ./ref_genomic.fna CM101310.1:68556902-685561280 > CM101310.1.68556902-68561280.fa

**Step 3. Retrieval of NUMT-associated HiFi reads**

    minimap2 -t 16 --secondary=no -ax map-hifi ./CM101310.1.68556902-68561280.fa ./HIFI_subreads.fastq.gz | samtools view -@ 16 -b -F4 -F 0x800 -o reads.HiFiMapped.bam
   
    samtools sort -@ 16 -o reads.HiFiMapped.sorted.bam reads.HiFiMapped.bam
   
    samtools index -@ 16 reads.HiFiMapped.sorted.bam

**Step 4. Read selection and visualization (optional but recommended)**

The sorted BAM file can be visualized in IGV to manually inspect read alignments.

Selection criteria:

Priority is given to HiFi reads spanning the entire NUMT region, including extended flanks when applicable.

Reads with higher alignment rate and continuity are preferred.

such as:


   <img width="765" height="314" alt="image" src="https://github.com/user-attachments/assets/e4aea820-d1d5-42f1-96ef-a4a5b969bd47" />


   
**Step 5. Get NUMT sequence**
    
	cat ref_mt.fa selected.HIFI.reads.fa > selected.numt-mt.fa

    mafft --auto  selected.numt-mt.fa > selected.numt-mt.mafft.fa  → MEGA



*Strategy 2*: No extension (recommended for rapid detection of more none-complete NUMTs)

For rapid screening of short NUMTs, only the aligned region is extracted without flanking extension:
 
    samtools faidx ./ref_genomic.fna CM101310.1:68558902-68559280 > CM101310.1.68558902-68559280.fa

**Step 3**. Retrieval of NUMT-associated HiFi reads

     minimap2 -t 16 --secondary=no -ax map-hifi ./CM101310.1.68558902-68559280.fa ./HIFI_subreads.fastq.gz | samtools view -@ 16 -b -F4 -F 0x800 -o reads.HiFiMapped.bam
   
	 samtools sort -@ 16 -o reads.HiFiMapped.sorted.bam reads.HiFiMapped.bam
   
	 samtools index -@ 16 reads.HiFiMapped.sorted.bam

   <img width="1558" height="174" alt="image" src="https://github.com/user-attachments/assets/9fe82883-915d-4a06-ab37-96c802e7f160" />

   
	 python3 ./script/extract_high_identity_reads.py -b reads.HiFiMapped.sorted.bam -o CM101310.1.68558902-68559280 -d pos

   This step generates two FASTA files:

*.numt.high98.full.fa — full-length reads

*.numt.high98.trim.fa — trimmed high-identity regions

**Step 4. Consensus NUMT sequence reconstruction**

The trimmed NUMT reads are combined with the mitochondrial reference and aligned
   
	 cat ref_mt.fa CM101310.1.68558902-68559280.numt.high98.trim.fa > CM101310.1.68558902-68559280.numt.high98.trim.mt.fa

     mafft --auto  CM101310.1.68558902-68559280.numt.high98.trim.mt.fa > CM101310.1.68558902-68559280.numt.high98.trim-mt.mafft.fa
   
     python3 ./script/numt_consensus_from_msa.py CM101310.1.68558902-68559280.numt.high98.trim-mt.mafft.fa ref_mt CM101310.1.68558902-68559280.numt


**Final output:**

**CM101310.1.68558902-68559280.numt.fa**


**Get the true mitochondrial sequence**

To get the true mitochondrial genome, the same workflow can be applied by replacing the nuclear NUMT reference with the mitochondrial reference, while keeping all other steps unchanged.

**Step 5. Check the Consensus NUMT sequence**

cat all numt sequences at  numts.fa

    python ./script/numt_diff_dedup_summary.py -n numts.fa -r ../PP646880.1.fa -o snow_leopard

we will get two file snow_leopard.diff_matrix.tsv snow_leopard.summary.tsv

**Step 6. Check the NUMT base in MT sequence**

need merge all mt sequence

python ./script/detect_mt_numt_contamination.py -m 157.ncbi.sequence.re_circularized_PP646880.2.fas -n snow_leopard.36.numts.fas -o snow_leopard




   

