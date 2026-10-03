#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import copy
import argparse
from Bio import SeqIO

# ----------------------
# ----------------------

parser = argparse.ArgumentParser(description="Detect NUMT sites from MSA")

parser.add_argument("aln_fa")
parser.add_argument("ref_name")
parser.add_argument("target_name")
parser.add_argument("prefix")

parser.add_argument("--flank", type=int, default=150)
parser.add_argument("--min_ratio", type=float, default=1)
parser.add_argument("--min_sites", type=int, default=150)

parser.add_argument("--F_region", type=int, default=1500)
parser.add_argument("--R_region", type=int, default=1500)

parser.add_argument("--min_del", type=int, default=10)

args = parser.parse_args()

aln_fa = args.aln_fa
ref_name = args.ref_name
target_name = args.target_name
prefix = args.prefix

flank = args.flank
min_ratio = args.min_ratio
min_sites = args.min_sites
F_region = args.F_region
R_region = args.R_region
min_del = args.min_del

out_raw = f"{prefix}.numt.tsv"
out_filtered = f"{prefix}.numt.filtered.tsv"
out_count = f"{prefix}.numt.count.tsv"

# ----------------------
# ----------------------

seqs = {}

for rec in SeqIO.parse(aln_fa, "fasta"):
    seqs[rec.id] = str(rec.seq).upper()

if ref_name not in seqs or target_name not in seqs:
    sys.exit("Error: ref or target not found in MSA")

ref = seqs[ref_name]
target = seqs[target_name]

sample_names = [s for s in seqs if s not in (ref_name, target_name)]

aln_len = len(ref)

# ----------------------
# ----------------------

ref_pos = []
pos = 0

for b in ref:

    if b != "-":
        pos += 1
        ref_pos.append(pos)
    else:
        ref_pos.append(None)

ref_length = pos

# ----------------------
# ----------------------

def ambiguous_match(a, b):
    return a == b


def flank_match(sample_seq, target_seq, idx):

    left_total = left_match = 0

    for i in range(max(0, idx - flank), idx):
        left_total += 1
        if ambiguous_match(sample_seq[i], target_seq[i]):
            left_match += 1

    left_ok = left_total >= min_sites and left_match / left_total >= min_ratio

    right_total = right_match = 0

    for i in range(idx + 1, min(aln_len, idx + flank + 1)):
        right_total += 1
        if ambiguous_match(sample_seq[i], target_seq[i]):
            right_match += 1

    right_ok = right_total >= min_sites and right_match / right_total >= min_ratio

    return left_ok or right_ok

# ----------------------
# ----------------------

with open(out_raw, "w") as out:

    header = ["pos", "ref", "target"] + sample_names
    out.write("\t".join(header) + "\n")

    for i in range(aln_len):

        if ref_pos[i] is None:
            continue

        r = ref[i]
        t = target[i]

        if r == t:
            continue

        if t == "-" and r != "-":
            t_out = "del"
        elif t != "-" and r == "-":
            t_out = "ins"
        else:
            t_out = t

        site_valid = False
        row = []

        for s in sample_names:

            sb = seqs[s][i]

            if sb == r:
                row.append("0")

            elif sb == t and flank_match(seqs[s], target, i):
                row.append("1")
                site_valid = True

            else:
                row.append("0")

        if site_valid:

            out.write(
                "\t".join([str(ref_pos[i]), r, t_out] + row) + "\n"
            )

# ----------------------
# ----------------------

with open(out_raw) as f:
    lines = [l.rstrip().split("\t") for l in f]

header = lines[0]
data = lines[1:]

sample_start = 3
filtered = copy.deepcopy(data)

# ----------------------
# ----------------------

for col in range(sample_start, len(header)):

    i = 0

    while i < len(filtered):

        if filtered[i][col] == "1" and filtered[i][2] == "del":

            block_start = i
            start_pos = int(filtered[i][0])

            while (
                i + 1 < len(filtered)
                and filtered[i + 1][col] == "1"
                and filtered[i + 1][2] == "del"
                and int(filtered[i + 1][0]) == int(filtered[i][0]) + 1
            ):
                i += 1

            block_end = i
            end_pos = int(filtered[i][0])

            block_len = block_end - block_start + 1

            if block_len >= min_del and (
                (start_pos <= F_region) or (end_pos >= ref_length - R_region)
            ):

                for k in range(block_start, block_end + 1):
                    filtered[k][col] = "0"

        i += 1

# ----------------------
# ----------------------

filtered = [
    row for row in filtered
    if any(row[col] == "1" for col in range(sample_start, len(header)))
]

# ----------------------
# ----------------------

for col in range(sample_start, len(header)):

    for i in range(len(filtered)):

        if filtered[i][col] != "1":
            continue

        prev_1 = i > 0 and filtered[i - 1][col] == "1"
        next_1 = i < len(filtered) - 1 and filtered[i + 1][col] == "1"

        if not (prev_1 or next_1):
            filtered[i][col] = "0"

# ----------------------
# ----------------------

for col in range(sample_start, len(header)):

    i = 0

    while i < len(filtered):

        if filtered[i][col] == "1":

            start = i

            while i + 1 < len(filtered) and filtered[i + 1][col] == "1":
                i += 1

            end = i
            block_len = end - start + 1

            if block_len == 2:

                pos1 = int(filtered[start][0])
                pos2 = int(filtered[end][0])

                if (pos2 - pos1) > 10:

                    filtered[start][col] = "0"
                    filtered[end][col] = "0"

        i += 1

# ----------------------
# ----------------------

with open(out_filtered, "w") as out:

    out.write("\t".join(header) + "\n")

    for r in filtered:
        out.write("\t".join(r) + "\n")

# ----------------------
# ----------------------

counts = {s: 0 for s in sample_names}

for row in filtered:

    for i, s in enumerate(sample_names):

        if row[sample_start + i] == "1":
            counts[s] += 1

with open(out_count, "w") as out:

    out.write("sample\tnumt_sites\n")

    for s in sample_names:
        out.write(f"{s}\t{counts[s]}\n")
