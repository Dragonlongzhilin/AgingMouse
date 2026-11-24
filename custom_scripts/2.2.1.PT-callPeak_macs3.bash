#!/bin/bash
#' @description Call peak with MACS3
#' @author Zhilin Long
#' @date 2026-8-10
# conda activate SingleCellAnalysis
histone="H3K4me1" # "H3K27ac" "H3K4me3" "H3K4me1"
tissues=("brainCB" "brainFC" "brainHip" "heart" "liver" "kidney" "BAT" "colon" "muscle" "spleen" "lung" "skin" "testis")
#### ---- call peak for each tissue ---- ####
## Young + Old
for tissue in "${tissues[@]}"; do
  echo "Analysis ${tissue}..."
  bam_path="~/longzhilin/public_html/AgingMouse/PairedTag_new/${tissue}/${histone}/cellType"
  mkdir -p $bam_path/peak
  cd $bam_path/peak

  bam_files=$(find $bam_path/bam -name "*.bam" | grep -vE "^Merge")
  if [ -z "$bam_files" ]; then
    echo "Error: No ${histone} BAM file found in $bam_path/bam";
  else
    fileName="Merge";
    if [[ "${histone}" =~ (H3K27me3|H3K9me3|H3K36me3|H3K4me1) ]]; then
      macs3 callpeak --nomodel -f BAM -g mm -t $bam_files -n ${fileName}.q.0001 --broad --broad-cutoff 0.0001 -q 0.0001 --keep-dup all --outdir .;
    else
      macs3 callpeak --nomodel -f BAM -g mm -t $bam_files -n ${fileName}.q.0001 -q 0.0001 --keep-dup all --outdir .;
    fi
  fi
done

## Young
for tissue in "${tissues[@]}"; do
  echo "Analysis ${tissue}..."
  bam_path="~/longzhilin/public_html/AgingMouse/PairedTag_new/${tissue}/${histone}/cellType_age"
  mkdir -p $bam_path/peak
  cd $bam_path/peak

  bam_files=$(find $bam_path/bam -name "*.bam" | grep "_Young.bam$" | grep -vE "^Merge")
  if [ -z "$bam_files" ]; then
    echo "Error: No ${histone} BAM file found in $bam_path/bam";
  else
    fileName="Merge_Young";
    if [[ "${histone}" =~ (H3K27me3|H3K9me3|H3K36me3|H3K4me1) ]]; then
      macs3 callpeak --nomodel -f BAM -g mm -t $bam_files -n ${fileName}.q.0001 --broad --broad-cutoff 0.0001 -q 0.0001 --keep-dup all --outdir .;
    else
      macs3 callpeak --nomodel -f BAM -g mm -t $bam_files -n ${fileName}.q.0001 -q 0.0001 --keep-dup all --outdir .;
    fi
  fi
done

## Old
for tissue in "${tissues[@]}"; do
  echo "Analysis ${tissue}..."
  bam_path="~/longzhilin/public_html/AgingMouse/PairedTag_new/${tissue}/${histone}/cellType_age"
  mkdir -p $bam_path/peak
  cd $bam_path/peak

  bam_files=$(find $bam_path/bam -name "*.bam" | grep "_Old.bam$" | grep -vE "^Merge")
  if [ -z "$bam_files" ]; then
    echo "Error: No ${histone} BAM file found in $bam_path/bam";
  else
    fileName="Merge_Old";
    if [[ "${histone}" =~ (H3K27me3|H3K9me3|H3K36me3|H3K4me1) ]]; then
      macs3 callpeak --nomodel -f BAM -g mm -t $bam_files -n ${fileName}.q.0001 --broad --broad-cutoff 0.0001 -q 0.0001 --keep-dup all --outdir .;
    else
      macs3 callpeak --nomodel -f BAM -g mm -t $bam_files -n ${fileName}.q.0001 -q 0.0001 --keep-dup all --outdir .;
    fi
  fi
done

#### ---- call peak for each cell type ---- ####
## Young + Old
for tissue in "${tissues[@]}"; do
  echo "Analysis ${tissue}..."
  bam_path="~/longzhilin/public_html/AgingMouse/PairedTag_new/${tissue}/${histone}/cellType"
  if [ ! -d "$bam_path/peak/macs3" ]; then
    mkdir -p $bam_path/peak/macs3
  fi
  cd $bam_path/peak/macs3

  bam_files=$(find "$bam_path/bam" -name "*.bam" | grep -vE "^Merge")
  cellTypes=$(echo "$bam_files" | sed 's/.*cellType\/bam\///' | sed 's/.bam//' | sort -u)
  for cellType in $cellTypes; do
    bam_file=$(echo "$bam_files" | grep "${cellType}.bam$")
    echo "Running analysis for $bam_file..."
    if [ -z "$bam_file" ]; then
      echo "Error: No ${histone} BAM file found in $bam_path/bam";
    else
      if [[ "${histone}" =~ (H3K27me3|H3K9me3|H3K36me3|H3K4me1) ]]; then
        macs3 callpeak --nomodel -f BAM -g mm -t $bam_file -n ${cellType}.q.0001 --broad --broad-cutoff 0.0001 -q 0.0001 --keep-dup all --outdir .;
      else
        macs3 callpeak --nomodel -f BAM -g mm -t $bam_file -n ${cellType}.q.0001 -q 0.0001 --keep-dup all --outdir .;
      fi
    fi
  done
done

## Young
for tissue in "${tissues[@]}"; do
  echo "Analysis ${tissue}..."
  bam_path="~/longzhilin/public_html/AgingMouse/PairedTag_new/${tissue}/${histone}/cellType_age"
  if [ ! -d "$bam_path/peak/macs3_Young" ]; then
    mkdir -p $bam_path/peak/macs3_Young
  fi
  cd $bam_path/peak/macs3_Young

  bam_files=$(find "$bam_path/bam" -name "*.bam" | grep "_Young.bam$" | grep -vE "^Merge")
  cellTypes=$(echo "$bam_files" | sed 's/.*cellType_age\/bam\///' | sed 's/.bam//' | sort -u)
  for cellType in $cellTypes; do
    bam_file=$(echo "$bam_files" | grep "${cellType}.bam$")
    echo "Running analysis for $bam_file..."
    if [ -z "$bam_file" ]; then
      echo "Error: No ${histone} BAM file found in $bam_path/bam";
    else
      if [[ "${histone}" =~ (H3K27me3|H3K9me3|H3K36me3|H3K4me1) ]]; then
        macs3 callpeak --nomodel -f BAM -g mm -t $bam_file -n ${cellType}.q.0001 --broad --broad-cutoff 0.0001 -q 0.0001 --keep-dup all --outdir .;
      else
        macs3 callpeak --nomodel -f BAM -g mm -t $bam_file -n ${cellType}.q.0001 -q 0.0001 --keep-dup all --outdir .;
      fi
    fi
  done
done

## Old
for tissue in "${tissues[@]}"; do
  echo "Analysis ${tissue}..."
  bam_path="~/longzhilin/public_html/AgingMouse/PairedTag_new/${tissue}/${histone}/cellType_age"
  if [ ! -d "$bam_path/peak/macs3_Old" ]; then
    mkdir -p $bam_path/peak/macs3_Old
  fi
  cd $bam_path/peak/macs3_Old

  bam_files=$(find "$bam_path/bam" -name "*.bam" | grep "_Old.bam$" | grep -vE "^Merge")
  cellTypes=$(echo "$bam_files" | sed 's/.*cellType_age\/bam\///' | sed 's/.bam//' | sort -u)
  for cellType in $cellTypes; do
    bam_file=$(echo "$bam_files" | grep "${cellType}.bam$")
    echo "Running analysis for $bam_file..."
    if [ -z "$bam_file" ]; then
      echo "Error: No ${histone} BAM file found in $bam_path/bam";
    else
      if [[ "${histone}" =~ (H3K27me3|H3K9me3|H3K36me3|H3K4me1) ]]; then
        macs3 callpeak --nomodel -f BAM -g mm -t $bam_file -n ${cellType}.q.0001 --broad --broad-cutoff 0.0001 -q 0.0001 --keep-dup all --outdir .;
      else
        macs3 callpeak --nomodel -f BAM -g mm -t $bam_file -n ${cellType}.q.0001 -q 0.0001 --keep-dup all --outdir .;
      fi
    fi
  done
done
