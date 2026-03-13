#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
NUMT consensus from MSA

Rules:
- ref sequence is used only for coordinates
- gap ('-') participates in consensus voting
- consensus threshold: >=80%
- unresolved sites are '-'
- output region = first and last non-gap base in consensus
  mapped to ref (gap excluded)

Output:
<output>.<start>-<end>.fa
>output.<start>-<end>
"""

from collections import Counter
from Bio import SeqIO
import sys


# ----------------------
# 参数
# ----------------------
if len(sys.argv) != 4:
    sys.exit("Usage: python numt_consensus_from_msa.py aln.fa ref_name output_prefix")

aln_fa, ref_name, output = sys.argv[1:4]

VALID = {"A", "C", "G", "T", "-"}

# ----------------------
# 读取 MSA
# ----------------------
seqs = {}
for rec in SeqIO.parse(aln_fa, "fasta"):
    seqs[rec.id] = str(rec.seq).upper()

if ref_name not in seqs:
    sys.exit("ERROR: reference sequence not found")

ref = seqs[ref_name]
numt_seqs = [s for k, s in seqs.items() if k != ref_name]

aln_len = len(ref)

# ----------------------
# 计算 ref 实际坐标（gap 不计）
# ----------------------
ref_coord = []
pos = 0
for b in ref:
    if b != "-":
        pos += 1
        ref_coord.append(pos)
    else:
        ref_coord.append(None)

# ----------------------
# 构建 NUMT 共识（alignment 坐标）
# ----------------------
consensus = []

for i in range(aln_len):
    r = ref[i]

    # ref gap → 不进入共识序列
    if r == "-":
        consensus.append(None)
        continue

    col = []
    for s in numt_seqs:
        b = s[i]
        if b in VALID:
            col.append(b)

    if not col:
        consensus.append("-")
        continue

    cnt = Counter(col)
    base, n = cnt.most_common(1)[0]

    if n / len(col) < 0.8:
        consensus.append("-")
    else:
        consensus.append(base)

# ----------------------
# 根据共识计算 start / end
# ----------------------
covered_coords = [
    ref_coord[i]
    for i, b in enumerate(consensus)
    if b not in (None, "-") and ref_coord[i] is not None
]

if not covered_coords:
    sys.exit("ERROR: consensus contains no non-gap bases")

start = min(covered_coords)
end = max(covered_coords)

region = f"{start}-{end}"
out_fa = f"{output}.{region}.fa"

# ----------------------
# 生成最终序列（去掉 ref gap）
# ----------------------
final_seq = [
    b for b in consensus
    if b is not None
]

final_seq = "".join(final_seq)

# ----------------------
# 输出 FASTA
# ----------------------
with open(out_fa, "w") as out:
    out.write(f">{output}.{region}\n")
    for i in range(0, len(final_seq), 60):
        out.write(final_seq[i:i+60] + "\n")

print(f"[OK] NUMT consensus written to: {out_fa}")

