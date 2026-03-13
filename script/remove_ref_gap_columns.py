#!/usr/bin/env python3
from Bio import SeqIO
import sys

aln = sys.argv[1]
ref_name = sys.argv[2]
out = sys.argv[3]

records = list(SeqIO.parse(aln, "fasta"))

ref = None
for r in records:
    if r.id == ref_name:
        ref = str(r.seq)
        break

if ref is None:
    sys.exit("Reference not found")

keep_cols = [i for i,b in enumerate(ref) if b != "-"]

with open(out, "w") as fw:
    for r in records:
        new_seq = "".join(r.seq[i] for i in keep_cols)
        fw.write(f">{r.id}\n{new_seq}\n")
