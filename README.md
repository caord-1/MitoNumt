# NUMT_get
This is a pepline for easy to get the true NUMT and MT sequence from the HIFI data.

1. get the ref_numt sequence
   need the mt and genome ref
   minimap2 -x asm5 ./ref_mt.fa ./ref_genomic.fna > *.paf
   Then, selected the ref_numt sequence from to .paf
   
