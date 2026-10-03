

#!/bin/bash

    ref_genome="ref_genomic.fna"
    ref_mt="ref_mt"
    SRR="HIFI_subreads.fastq.gz"
    species="species"

minimap2 -x asm5 ./${ref_mt}.fa ./${ref_genome} > ${ref_mt}-${ref_genome}.paf

awk 'BEGIN{OFS="\t"}
$11 > 150 {
    if ($5 == "+") $5 = "pos";
    else if ($5 == "-") $5 = "neg";

    key = $8 FS $9 FS $10;

    if (!(key in seen)) {
        seen[key] = 1;
        print $1, $3, $4, $5;
    }
}' ${ref_mt}-${ref_genome}.paf > numts.list

while read sample st ed dir; do

    prefix="${sample}.${st}-${ed}"
    bam="${prefix}.reads.HiFiMapped.bam"
    sbam="${prefix}.reads.HiFiMapped.sorted.bam"
    
    echo "[INFO] Processing ${prefix}"

    # ---------- 1. 比对（仅在 bam 不存在时） ----------
    if [ ! -s "$bam" ]; then
        echo "[INFO] BAM 不存在，开始比对"

        samtools faidx ./${ref_genome} \
            ${sample}:${st}-${ed} > ${prefix}.fa

        minimap2 -t 16 --secondary=no -ax map-hifi \
            ${prefix}.fa ./${SRR} \
        | samtools view -@ 16 -b -F4 -F 0x800 \
        > "$bam"
    else
        echo "[INFO] BAM 已存在，跳过比对"
    fi

    # ---------- 2. 排序 ----------
    if [ ! -s "$sbam" ]; then
        samtools sort -@ 16 -o "$sbam" "$bam"
    fi

    # ---------- 3. 建索引 ----------
    if [ ! -s "${sbam}.bai" ]; then
        samtools index -@ 16 "$sbam"
    fi

    # ---------- 4. 后续 NUMT 分析 ----------
    python ./script/extract_high_identity_reads.py \
        -b "$sbam" \
        -o ${prefix}.numt.high98 \
        -d "$dir"

    seqkit head -n 50 \
        ${prefix}.numt.high98.trim.fa \
        > ${prefix}.numt.high98.trim.50.fa

    cat ${prefix}.numt.high98.trim.50.fa ./${ref_mt}.fa \
        > ${prefix}.numt.high98.trim-numt-mt.mt.fa

    mafft --auto \
        ${prefix}.numt.high98.trim-numt-mt.mt.fa \
        > ${prefix}.numt.high98.trim-numt-mt.mt.mafft.fa

    python ./script/numt_consensus_from_msa.py \
        ${prefix}.numt.high98.trim-numt-mt.mt.mafft.fa \
        ${ref_mt} \
        ${species}.numt

done < numts.list

cat ${species}.numt* > numts.fa

python ./script/numt_ref_dif.py -n numts.fa -r ./${ref_mt}.fa -o ${species}


