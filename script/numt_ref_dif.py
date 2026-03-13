#!/usr/bin/env python3
import argparse
import subprocess
import tempfile
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO


def run_mafft_pair(ref_rec, numt_rec):
    """
    Run MAFFT on ref + single numt
    Return aligned ref_seq, numt_seq
    """
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        in_fa = tmp / "in.fa"
        out_fa = tmp / "out.fa"

        SeqIO.write([ref_rec, numt_rec], in_fa, "fasta")

        subprocess.run(
            ["mafft", "--auto", str(in_fa)],
            stdout=open(out_fa, "w"),
            stderr=subprocess.PIPE,
            check=True
        )

        aln = list(SeqIO.parse(out_fa, "fasta"))
        if len(aln) != 2:
            raise RuntimeError("MAFFT did not return 2 sequences")

        if aln[0].id == ref_rec.id:
            return str(aln[0].seq), str(aln[1].seq)
        else:
            return str(aln[1].seq), str(aln[0].seq)


def ungapped_len(seq):
    return len(seq.replace("-", ""))


def main(numt_fa, ref_fa, prefix):
    ref_rec = next(SeqIO.parse(ref_fa, "fasta"))
    numt_records = list(SeqIO.parse(numt_fa, "fasta"))

    diff_matrix = defaultdict(dict)
    summary = []

    for numt in numt_records:
        ref_aln, numt_aln = run_mafft_pair(ref_rec, numt)

        ref_pos = 0
        diff_count = 0

        for r, n in zip(ref_aln, numt_aln):
            if r != "-":
                ref_pos += 1

            if r == "-" or n == "-":
                continue

            if r != n:
                diff_matrix[ref_pos][numt.id] = n
                diff_count += 1

        numt_len = ungapped_len(str(numt.seq))
        diff_rate = diff_count / numt_len if numt_len > 0 else 0

        summary.append(
            (numt.id, numt_len, diff_count, diff_rate)
        )

    # ---------- diff matrix ----------
    matrix_out = f"{prefix}.diff_matrix.tsv"
    numt_ids = [r.id for r in numt_records]

    with open(matrix_out, "w") as f:
        f.write("Ref_pos\tRef_base\t" + "\t".join(numt_ids) + "\n")

        for pos in sorted(diff_matrix):
            if not diff_matrix[pos]:
                continue

            ref_base = ref_rec.seq[pos - 1]
            row = [diff_matrix[pos].get(i, "") for i in numt_ids]
            f.write(
                f"{pos}\t{ref_base}\t" + "\t".join(row) + "\n"
            )

    # ---------- summary ----------
    summary_out = f"{prefix}.summary.tsv"
    with open(summary_out, "w") as f:
        f.write("numt_id\tlength\tdiff_count\tdiff_rate\n")
        for sid, l, c, r in summary:
            f.write(f"{sid}\t{l}\t{c}\t{r:.6f}\n")

    print(f"[OK] Summary     : {summary_out}")
    print(f"[OK] Diff matrix : {matrix_out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="NUMT vs mt pairwise diff statistics using MAFFT"
    )
    parser.add_argument("-n", "--numt", required=True, help="NUMT fasta (multi)")
    parser.add_argument("-r", "--ref", required=True, help="Reference mt fasta")
    parser.add_argument("-o", "--out", required=True, help="Output prefix")

    args = parser.parse_args()
    main(args.numt, args.ref, args.out)
