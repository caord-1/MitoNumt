# Recircularize all mtGenome sequences using forward 100bp of reference 

#python circularized.py --seed_len 100 --min_ratio 0.6 ref.fa input.fa output.fa 

#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(
        description="Recircularize sequences using reference seed"
    )
    parser.add_argument("ref", help="Reference fasta file")
    parser.add_argument("merge_file", help="Input fasta file to recircularize")
    parser.add_argument("out", help="Output fasta file")
    parser.add_argument("--seed_len", type=int, default=100,
                        help="Seed length from reference start (default: 100)")
    parser.add_argument("--min_ratio", type=float, default=0.6,
                        help="Minimum match ratio to accept position (default: 0.6)")
    return parser.parse_args()

def find_best_pos(seq, seed, min_ratio):
    best_pos = -1
    best_match = 0
    seed_len = len(seed)

    for i in range(len(seq) - seed_len + 1):
        match = sum(seq[i + j] == seed[j] for j in range(seed_len))
        if match > best_match:
            best_match = match
            best_pos = i

    if best_match < seed_len * min_ratio:
        return -1
    return best_pos

def main():
    args = parse_args()

    # 读取参考序列
    ref_record = next(SeqIO.parse(args.ref, "fasta"))
    ref_seq = str(ref_record.seq)

    if len(ref_seq) < args.seed_len:
        raise ValueError("Reference length shorter than seed_len")

    seed = ref_seq[:args.seed_len]

    with open(args.out, "w") as fw:
        for rec in SeqIO.parse(args.merge, "fasta"):
            seq = str(rec.seq)

            if len(seq) < args.seed_len:
                print(f"⚠️ {rec.id}: 序列长度小于 seed_len，跳过重排")
                new_seq = seq
            else:
                pos = find_best_pos(seq, seed, args.min_ratio)
                if pos == -1:
                    print(f"⚠️ {rec.id}: 无法定位 ref 起点，保持原序列")
                    new_seq = seq
                else:
                    new_seq = seq[pos:] + seq[:pos]

            fw.write(f">{rec.id}\n{new_seq}\n")

    print("✔️ 已完成重环化:", args.out)

if __name__ == "__main__":
    main()

