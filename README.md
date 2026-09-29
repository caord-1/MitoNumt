# MitoNumt
Overview

This pipeline provides a practical and reproducible workflow for recovering mitochondrial DNA (mtDNA) sequences, identifying both complete and incomplete NUMTs from PacBio HiFi sequencing data, and detecting NUMT-derived contamination in the assembled mitochondrial genome.

  ![fig](numt_schematic_diagram.png)

*Requirements*

**python3 minimap2 samtools mafft igv**

Mitochondrial reference genome (ref_mt.fa)

Nuclear genome reference (ref_genomic.fna)

HIFI reads (HIFI_subreads.fastq.gz)

**Step 1. Get the true mitochondrial sequence**

1.1 Align the HIFI data to the mitochondrial genome

    minimap2 -t 16 --secondary=no -ax map-hifi ./ref_mt.fa ./HIFI_subreads.fastq.gz | samtools view -@ 16 -b -F4 -F 0x800 -o mt.reads.HiFiMapped.bam
   
    samtools sort -@ 16 -o mt.reads.HiFiMapped.sorted.bam mt.reads.HiFiMapped.bam
   
    samtools index -@ 16 mt.reads.HiFiMapped.sorted.bam

1.2 Reads selection and visualization 

The sorted BAM file can be visualized in IGV to manually inspect read alignments.

Selection criteria:

(1) Reads spanning the reference mtgenome were preferentially selected; 
(2) if no single read spanned the reference, overlapping reads were required to collectively bridge it; 
(3) sequences extending beyond either end of the reference mtgenome were required to match the corresponding sequence at the opposite end, as verified in MEGA; 
(4) among reads covering the same region, those with the lowest value of edit distance to the reference (NM) were preferentially selected to maximize sequence similarity to the reference mtgenome. 

Following these criteria, approximately 3–5 high-quality HiFi reads were sufficient to support a mitochondrial sequence.

Black arrows represent the recommended reads
  <img width="1200" height="690" alt="image" src="https://github.com/user-attachments/assets/21e99e2d-cc5c-4a98-bb74-0d4d44a2a72a" />

Sometimes the selected ref mt also contains contaminated regions. The criteria for selecting reads are the same, but there is a slight difference in the mapping rate part.such as SRR28532995.4886972, SRR28532995.296568 ...

  <img width="1129" height="563" alt="image" src="https://github.com/user-attachments/assets/02f844e1-36c8-4685-a4c6-ffd3f642767c" />

Next, we will compare the selected 3-4 reads with the reference sequence, and then perform visualization and correction in Mega. This will enable us to obtain the True mitochondrial sequence (True_ref_mt.fa).


**Step 2. Identification of candidate NUMT regions in the nuclear genome**


2.1 Candidate NUMT loci are identified by aligning the True_mitochondrial genome to the nuclear genome using minimap2:

    minimap2 -x asm5 ./True_ref_mt.fa ./ref_genomic.fna > *.paf

The resulting PAF file is then inspected to extract nuclear regions showing significant homology to the mitochondrial genome, which are considered putative NUMT loci.

   <img width="1034" height="673" alt="image" src="https://github.com/user-attachments/assets/40142038-4a1b-4d96-befb-4ce1051c1533" />

2.2 Extraction of reference NUMT regions


***Extension strategy***: recommended for detecting complete or near-complete NUMTs

To get potentially full-length NUMTs, the candidate NUMT region is extended by 2 kb upstream and downstream from the alignment coordinates:
   
	samtools faidx ./ref_genomic.fna CM109601.1:5456425-5717163 > CM109601.1:5456425-5717163.fa

2.2.1 Retrieval of NUMT-associated HiFi reads

    minimap2 -t 16 --secondary=no -ax map-hifi ./CM109601.1:5456425-5717163.fa ./HIFI_subreads.fastq.gz | samtools view -@ 16 -b -F4 -F 0x800 -o numt.reads.HiFiMapped.bam
   
    samtools sort -@ 16 -o numt.reads.HiFiMapped.sorted.bam numt.reads.HiFiMapped.bam
   
    samtools index -@ 16 numt.reads.HiFiMapped.sorted.bam

2.2.2 Read selection and visualization

The sorted BAM file can be visualized in IGV to manually inspect read alignments.

Selection criteria:

Priority is given to HiFi reads spanning the nuclear genomic region and the NUMT region.

Reads with higher alignment rate and continuity are preferred.

such as:
   
   <img width="1081" height="252" alt="image" src="https://github.com/user-attachments/assets/6eec23a2-55fa-4c03-baef-1b9787d286c0" />


 If the reads spanning the nuclear genomic region and the NUMT region do not form a complete circular NUMT, reads can still be searched for within the intervening region. However, these reads must be alignable to the previously identified non-circular NUMT sequence. This can be quickly verified in MEGA
   
   <img width="1081" height="252" alt="image" src="https://github.com/user-attachments/assets/46b9c226-9c0d-4611-a1a2-d9fa9955640c" />

2.2.3 Get NUMT sequence
    
	cat True_ref_mt.fa selected.HIFI.reads.fa > selected.numt-mt.fa

    mafft --auto  selected.numt-mt.fa > selected.numt-mt.mafft.fa  → MEGA

generate the complete.numt.fa


***Non-extension strategy***: recommended for rapid detection of incomplete NUMTs

2.2.1 For rapid screening of short NUMTs, only the each aligned region is extracted without flanking extension:

  <img width="1037" height="129" alt="image" src="https://github.com/user-attachments/assets/7aae8f7b-9c60-4a62-8fde-e8b0eb821c1d" />

    samtools faidx ./ref_genomic.fna CM109591.1:68539524-68540092 > CM109591.1:68539524-68540092.fa

2.2.2 Retrieval of NUMT-associated HiFi reads

     minimap2 -t 16 --secondary=no -ax map-hifi ./CM109591.1:68539524-68540092.fa ./HIFI_subreads.fastq.gz | samtools view -@ 16 -b -F4 -F 0x800 -o CM109591.1:68539524-68540092.reads.HiFiMapped.bam
   
	 samtools sort -@ 16 -o CM109591.1:68539524-68540092.reads.HiFiMapped.sorted.bam CM109591.1:68539524-68540092.reads.HiFiMapped.bam
   
	 samtools index -@ 16 CM109591.1:68539524-68540092.reads.HiFiMapped.sorted.bam

	 python3 ./script/extract_high_identity_reads.py -b reads.HiFiMapped.sorted.bam -o CM109591.1:68539524-68540092 -d pos

   This step generates two FASTA files:

*.numt.high98.full.fa — full-length HIFI reads

*.numt.high98.trim.fa — trimmed high-identity regions

**Step 3. Consensus NUMT sequence reconstruction**

The trimmed NUMT reads are combined with the mitochondrial reference and aligned
The long reads obtained using the extension strategy can also be processed with this script, provided that the portions aligned to the nuclear genome are trimmed.
   
	 cat True_ref_mt.fa CM109591.1:68539524-68540092.numt.high98.trim.fa > CM109591.1:68539524-68540092.numt.high98.trim.mt.fa

     mafft --auto  CM109591.1:68539524-68540092.numt.high98.trim.mt.fa > CM109591.1:68539524-68540092.numt.high98.trim-mt.mafft.fa
   
     python3 ./script/numt_consensus_from_msa.py CM109591.1:68539524-68540092.numt.high98.trim-mt.mafft.fa ref_mt prefix.numt


**Final output:**

**prefix.numt.2597-4082.fa** : This is a incomplete numt sequence. 2597-4082 is the start-end of NUMT at True_ref_mt.fa

Perform 2.2-3 operations on each small segment, and this will result in numerous numt sequences.

 Collect all the information of the numt sequences**

    cat *.numt.* complete.numt.fa > numts.fa

    python ./script/numt_ref_dif.py -n numts.fa -r ../Ture_ref.mt.fa -o prefix

we will get two file ***prefix.diff_matrix.tsv*** and ***prefix.summary.tsv***

Delete the numt sequences with a difference rate less than 0.01 and  fewer than 10 variant sites, and then remove the duplicate numt sequences.
	
**Step 4. Check the NUMT base in MT sequence**

Run the following command for each NUMT sequence.

    cat all_mt_sequence True_ref.mt.fa numt.fa > prefix.fa

    mafft --auto prefix.fa > prefix.mafft.fa

    python remove_ref_gap_columns.py prefix.mafft.fa ref_name prefix.mafft.no_gap.fa

    python ./script/numt_from_msa.py prefix.mafft.no_gap.fa True_ref.mt_name numt_name prefix 

The parameters used for NUMT site detection (numt_from_msa.py) were: flank = 150, min_ratio = 1, min_sites = 150, F_region = 1500, R_region = 1500, and min_del = 10

**flank**: Length of the flanking sequence used for comparison on each side of the candidate NUMT site; 

**min_ratio**: Minimum flanking sequence matching ratio required to classify a site as NUMT-derived; min_sites: 

**min_sites**: Minimum number of aligned bases required for flanking sequence comparison.
	
flank, min_ratio, and min_sites are the primary parameters used to determine whether a candidate site is supported as a NUMT-derived site based on the similarity of its flanking sequences. a flank length of 150 bp provides a practical default for identifying linked NUMT-derived sites in typical Illumina datasets. For datasets generated using longer reads, a larger flank value may improve specificity, whereas shorter reads may benefit from a smaller value to maintain detection sensitivity.
	
**F_region**: Length of the 5’ terminal region in which long contiguous deletion blocks are excluded to reduce false positives caused by terminal length variation (for example: the D-loop region). 
**R_region**: Length of the 3’ terminal region in which long contiguous deletion blocks are excluded to reduce false positives caused by terminal length variation. 
**min_del**: Minimum length of a contiguous deletion block to be removed from the defined terminal regions.
	
F_region, R_region, and min_del function together to reduce false-positive detections caused by long contiguous deletion blocks. Here, min_del specifies the minimum length of a contiguous deletion block, while F_region and R_region define the 5′ and 3′ terminal regions in which such deletion blocks are excluded from NUMT detection. By default, deletion blocks of at least 10 bp located within the first or last 1,500 bp of the mitochondrial alignment are not considered NUMT-derived sites. 
	

Finally, we will get two files.    

**prefix.numt.filtered.tsv**：The specific locations of NUMT contamination sites in each mitochondrion after filtration.

**prefix.numt.count.tsv**：The statistics of the number of NUMT contamination sites in each mitochondrion after filtration.

