#!/bin/bash
# https://github.com/7lagrange/scPT-Natural-Aging/blob/main/scripts/generate_dna_counts.sh
set -euo pipefail
# conda activate scDNA

tissues=("brainFC" "brainHip" "brainCB" "kidney" "lung" "heart" "BAT" "colon" "muscle" "liver" "testis" "spleen")
histones=("H3K27ac" "H3K27me3" "H3K9me3")
BASE_DIR="/storage/zhangyanxiaoLab/qihongjian/projects/paired_seq_tag"
DATA_DIR="$BASE_DIR/data/merge"
PIPELINE_DIR="/storage/zhangyanxiaoLab/qihongjian/github/Paired-Tag-Pipeline-modified/pipelines"
OUT_DIR="/storage/zhangyanxiaoLab/longzhilin/Projects/AgingMouse/PairedTag/3.histoneLevel"
THREADS=64

get_peak_file() {
  local histone="$1"
  case "$histone" in
    H3K27ac)
      echo "/storage/zhangyanxiaoLab/suzhuojie/projects/pairedtag_27ac_4me1/data/peak_merge/spm5_all_chr/all_celltype_summit_all_chr_merge.bed"
      ;;
    H3K27me3)
      echo "/storage/zhangyanxiaoLab/longzhilin/public_html/AgingMouse/PairedTag_new/MergeIGV/H3K27me3/H3K27me3_tissue_peak.bed"
      ;;
    H3K9me3)
      echo "/storage/zhangyanxiaoLab/longzhilin/public_html/AgingMouse/PairedTag_new/MergeIGV/H3K9me3/H3K9me3_tissue_peak.bed"
      ;;
    *)
      echo "[ERROR] Unsupported histone mark: $histone" >&2
      return 1
      ;;
  esac
}

for histone in "${histones[@]}"; do
  peak_file="$(get_peak_file "$histone")"
  if [[ ! -f "$peak_file" ]]; then
    echo "[WARN] BED file not found for $histone: $peak_file, skipping this histone."
    continue
  fi

  for tissue in "${tissues[@]}"; do
    echo "======================================="
    echo "[INFO] Processing $tissue - $histone"
    echo "======================================="

    if [[ ! -d "$OUT_DIR/$histone/peakMatrix" ]]; then
      mkdir -p "$OUT_DIR/$histone/peakMatrix"
    fi
    cd "$OUT_DIR/$histone/peakMatrix"

    merged_bam="${tissue}-${histone}-merged.bam"
    matrix_name="${tissue}-${histone}"

    echo "[INFO] Merging BAM..."
    files=("$DATA_DIR"/*"${tissue}"*"${histone}"*/dna_bam/histone_dna_bam/*"${histone}"*.bam)

    n=${#files[@]}
    if [ $n -eq 0 ]; then
        echo "[WARN] No BAM files found for $tissue $histone, skipping..."
        continue
    elif [ $n -eq 1 ]; then
        echo "[INFO] Only one BAM file found, copying..."
        cp "${files[@]}" "$merged_bam"
        samtools index -@ $THREADS "$merged_bam"
    else
        echo "[INFO] Merging $n BAM files..."
        samtools merge -f -@ $THREADS "$merged_bam" "${files[@]}"
        samtools index -@ $THREADS "$merged_bam"
    fi

    bash "$PIPELINE_DIR/bam2mtx.sh" \
      -n "$matrix_name" \
      -g "$peak_file" \
      -i "$merged_bam" \
      -l DNA \
      -o "$matrix_name" \
      --thread "$THREADS"

    rm -f "$merged_bam" "${merged_bam}.bai" "${matrix_name}_CountMatrix.txt"
    echo "[INFO] Done: $tissue - $histone"
  done
done

echo "[INFO] ALL DONE"
