#!/usr/bin/env python3
#' @description Extract Paired-Tag DNA insert from R2 reads and filter bad reads
#' @Author: Zhilin Long
#' @date 2026-9-7
"""
Extract the genomic insert of Paired-Tag DNA R2 reads and drop low-quality-structure reads.

R2 read structure (DNA, following Paired-Tag filter_merge_bam.sh):
    [~17bp cell barcode / UMI] GTGGCC [genomic insert + 7bp barcode] ATCCAC [barcode + polyA] CCTGCA...

Reads are filtered into classes:
    - Contamination : contains T[ATCG]{4}CCTGCA  (primer dimer / mis-ligation)
    - polyG         : starts with >=15 G
    - no-structure  : missing GTGGCC...ATCCAC  (Random)
    - short         : insert shorter than min_insert
    - kept          : valid structure, insert length >= min_insert

Usage:
    python extractPairedTagDNAinsert.py <in.R2.fastq.gz> <out.insert.fastq.gz> [min_insert=20] [barcode_3p=7]
"""
import sys, gzip, re

def main():
    if len(sys.argv) < 3:
        sys.stderr.write(
            "Usage: extractPairedTagDNAinsert.py <in.R2.fastq.gz> <out.insert.fastq.gz> "
            "[min_insert=20] [barcode_3p=7]\n")
        sys.exit(1)

    inpath = sys.argv[1]
    outpath = sys.argv[2]
    min_insert = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    barcode_3p = int(sys.argv[4]) if len(sys.argv) > 4 else 7

    ADAPTOR_5 = "GTGGCC"    # 5' adaptor right before the insert
    LINKER_3 = "ATCCAC"     # 3' linker right after the insert + barcode
    CONTAM = re.compile(r"T[ATCG]{4}CCTGCA")
    POLYG = re.compile(r"^G{15,}")

    opener = gzip.open if inpath.endswith(".gz") else open
    writer = gzip.open if outpath.endswith(".gz") else open

    n_total = n_keep = n_contam = n_polyg = n_nostruct = n_short = 0
    with opener(inpath, "rt") as fin, writer(outpath, "wt") as fout:
        while True:
            name = fin.readline().rstrip("\n")
            if not name:
                break
            seq = fin.readline().rstrip("\n")
            plus = fin.readline().rstrip("\n")
            qual = fin.readline().rstrip("\n")
            n_total += 1

            if CONTAM.search(seq):
                n_contam += 1
                continue
            if POLYG.match(seq):
                n_polyg += 1
                continue

            i5 = seq.find(ADAPTOR_5)
            if i5 < 0:
                n_nostruct += 1
                continue
            i3 = seq.find(LINKER_3, i5 + len(ADAPTOR_5))
            if i3 < 0:
                n_nostruct += 1
                continue

            start = i5 + len(ADAPTOR_5)
            end = i3 - barcode_3p   # drop barcode_3p bp of cell barcode before ATCCAC
            if end - start < min_insert:
                n_short += 1
                continue

            fout.write(name + "\n" + seq[start:end] + "\n" + plus + "\n" + qual[start:end] + "\n")
            n_keep += 1

    sys.stderr.write(
        f"[extractPairedTagDNAinsert] {inpath}\n"
        f"  total={n_total} kept={n_keep} contam={n_contam} "
        f"polyG={n_polyg} noStructure={n_nostruct} short={n_short}\n")

if __name__ == "__main__":
    main()
