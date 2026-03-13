#!/usr/bin/env python3
import pysam
import argparse

# -----------------------------
# 反向互补
# -----------------------------
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

def revcomp(seq):
    return seq.translate(COMP)[::-1]


# -----------------------------
# 计算 aligned length
# -----------------------------
def aligned_length(aln):
    length = 0
    for op, l in aln.cigartuples:
        if op in (0, 1, 2, 7, 8):  # M I D = X
            length += l
    return length


# -----------------------------
# 主程序
# -----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Extract high-identity reads from BAM and trim soft clips"
    )
    parser.add_argument("-b", "--bam", required=True, help="sorted BAM file")
    parser.add_argument("-o", "--out", required=True, help="output prefix")
    parser.add_argument(
        "-d", "--direction",
        choices=["pos", "neg"],
        required=True,
        help="desired output direction (pos or neg)"
    )
    parser.add_argument(
        "--min_identity",
        type=float,
        default=0.98,
        help="minimum alignment identity (default: 0.98)"
    )

    args = parser.parse_args()

    bam = pysam.AlignmentFile(args.bam, "rb")

    fout_full = open(f"{args.out}.full.fa", "w")
    fout_trim = open(f"{args.out}.trim.fa", "w")

    for aln in bam.fetch(until_eof=True):

        if aln.is_unmapped:
            continue
        if aln.is_secondary or aln.is_supplementary:
            continue

        nm = aln.get_tag("NM") if aln.has_tag("NM") else None
        if nm is None:
            continue

        aln_len = aligned_length(aln)
        if aln_len == 0:
            continue

        identity = 1 - nm / aln_len
        if identity < args.min_identity:
            continue

        # read 序列
        seq = aln.query_sequence
        if seq is None:
            continue

        # soft clip
        left_clip = 0
        right_clip = 0
        cig = aln.cigartuples

        if cig[0][0] == 4:
            left_clip = cig[0][1]
        if cig[-1][0] == 4:
            right_clip = cig[-1][1]

        trimmed_seq = seq[left_clip: len(seq) - right_clip]

        # 方向处理
        if args.direction == "neg":
            seq = revcomp(seq)
            trimmed_seq = revcomp(trimmed_seq)

        # FASTA header
        header = (
            f">{aln.query_name} "
            f"identity={identity:.4f} "
            f"NM={nm} "
            f"aln_len={aln_len} "
            f"strand={'-' if aln.is_reverse else '+'}"
        )

        fout_full.write(header + "\n" + seq + "\n")
        fout_trim.write(header + "\n" + trimmed_seq + "\n")

    bam.close()
    fout_full.close()
    fout_trim.close()


if __name__ == "__main__":
    main()
