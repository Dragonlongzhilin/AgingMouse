#!/bin/bash
#' @description Call peak with SICER2
#' @author Zhilin Long
#' @date 2026-8-10

# conda activate SingleCellAnalysis
histone="H3K9me3" # "H3K9me3" "H3K27me3" "H3K36me3"
tissues=("brainCB" "brainFC" "brainHip" "heart" "liver" "kidney" "BAT" "colon" "muscle" "spleen" "lung" "skin" "testis")
# H3K4me1: window_size=1000, gap_size=3000
#### ---- call peak for each tissue ---- ####
## Young + Old
for tissue in "${tissues[@]}"; do
  echo "Analysis ${tissue}..."
  bam_path="~/longzhilin/public_html/AgingMouse/PairedTag_new/${tissue}/${histone}/cellType"
  mkdir -p $bam_path/peak
  cd $bam_path/peak

  window_size=5000
  gap_size=10000
  bam_files=$(find $bam_path/bam -name "*.bam" | grep -vE "^Merge")
  if [ -z "$bam_files" ]; then
    echo "Error: No ${histone} BAM file found in $bam_path/bam";
  else
    fileName="Merge"
    samtools merge -@ 12 "${fileName}.bam" $bam_files
    echo "Call peak for ${fileName}:${bam_files} using sicer with -w ${window_size} -g ${gap_size}";
    ~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t "${fileName}.bam" -s mm10 -w $window_size -g $gap_size -cpu 20;
    scoreisland="${fileName}-W${window_size}-G${gap_size}.scoreisland";
    if [ ! -f "${scoreisland}" ]; then
      gap_size=5000
      echo "Re-Call peak for ${fileName} using sicer with -w ${window_size} -g ${gap_size}";
      ~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t "${fileName}.bam" -s mm10 -w $window_size -g $gap_size -cpu 20;
      scoreisland="${fileName}-W${window_size}-G${gap_size}.scoreisland";
    fi

    if [ ! -f "${scoreisland}" ]; then
      echo "Failed to Call peak for ${fileName} using sicer!!!";
      rm *.wig;
    else
      awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$4,"."}' $scoreisland > "${fileName}-W${window_size}-G${gap_size}.bed";
      mv $scoreisland "${fileName}-W${window_size}-G${gap_size}.scoreisland";
    fi
    rm "${fileName}.bam" "${fileName}.bed"
  fi
done

## Young
for tissue in "${tissues[@]}"; do
  echo "Analysis ${tissue}..."
  bam_path="~/longzhilin/public_html/AgingMouse/PairedTag_new/${tissue}/${histone}/cellType_age"
  mkdir -p $bam_path/peak
  cd $bam_path/peak

  window_size=5000
  gap_size=10000
  bam_files=$(find $bam_path/bam -name "*.bam" | grep "_Young.bam$" | grep -vE "^Merge")
  if [ -z "$bam_files" ]; then
    echo "Error: No ${histone} BAM file found in $bam_path/bam";
  else
    fileName="Merge_Young"
    samtools merge -@ 12 "${fileName}.bam" $bam_files
    echo "Call peak for ${fileName} using sicer with -w ${window_size} -g ${gap_size}";
    ~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t "${fileName}.bam" -s mm10 -w $window_size -g $gap_size -cpu 20;
    scoreisland="${fileName}-W${window_size}-G${gap_size}.scoreisland";
    if [ ! -f "${scoreisland}" ]; then
      gap_size=5000
      echo "Re-Call peak for ${fileName} using sicer with -w ${window_size} -g ${gap_size}";
      ~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t "${fileName}.bam" -s mm10 -w $window_size -g $gap_size -cpu 20;
      scoreisland="${fileName}-W${window_size}-G${gap_size}.scoreisland";
    fi

    if [ ! -f "${scoreisland}" ]; then
      echo "Failed to Call peak for ${fileName} using sicer!!!";
      rm *.wig;
    else
      awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$4,"."}' $scoreisland > "${fileName}-W${window_size}-G${gap_size}.bed";
      mv $scoreisland "${fileName}-W${window_size}-G${gap_size}.scoreisland";
    fi
    rm "${fileName}.bam" "${fileName}.bed"
  fi
done

## Old
for tissue in "${tissues[@]}"; do
  echo "Analysis ${tissue}..."
  bam_path="~/longzhilin/public_html/AgingMouse/PairedTag_new/${tissue}/${histone}/cellType_age"
  mkdir -p $bam_path/peak
  cd $bam_path/peak

  window_size=5000
  gap_size=10000
  bam_files=$(find $bam_path/bam -name "*.bam" | grep "_Old.bam$" | grep -vE "^Merge")
  if [ -z "$bam_files" ]; then
    echo "Error: No ${histone} BAM file found in $bam_path/bam";
  else
    fileName="Merge_Old"
    samtools merge -@ 12 "${fileName}.bam" $bam_files
    echo "Call peak for ${fileName} using sicer with -w ${window_size} -g ${gap_size}";
    ~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t "${fileName}.bam" -s mm10 -w $window_size -g $gap_size -cpu 20;
    scoreisland="${fileName}-W${window_size}-G${gap_size}.scoreisland";
    if [ ! -f "${scoreisland}" ]; then
      gap_size=5000
      echo "Re-Call peak for ${fileName} using sicer with -w ${window_size} -g ${gap_size}";
      ~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t "${fileName}.bam" -s mm10 -w $window_size -g $gap_size -cpu 20;
      scoreisland="${fileName}-W${window_size}-G${gap_size}.scoreisland";
    fi

    if [ ! -f "${scoreisland}" ]; then
      echo "Failed to Call peak for ${fileName} using sicer!!!";
      rm *.wig;
    else
      awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$4,"."}' $scoreisland > "${fileName}-W${window_size}-G${gap_size}.bed";
      mv $scoreisland "${fileName}-W${window_size}-G${gap_size}.scoreisland";
    fi
    rm "${fileName}.bam" "${fileName}.bed"
  fi
done

#### ---- call peak for each cell type ---- ####
## Young + Old
for tissue in "${tissues[@]}"; do
  echo "Analysis ${tissue}..."
  bam_path="~/longzhilin/public_html/AgingMouse/PairedTag_new/${tissue}/${histone}/cellType"
  if [ ! -d "$bam_path/peak/sicer" ]; then
    mkdir -p $bam_path/peak/sicer
  fi
  cd $bam_path/peak/sicer

  window_size=5000
  gap_size=10000
  bam_files=$(find "$bam_path/bam" -name "*.bam" | grep -vE "^Merge")
  cellTypes=$(echo "$bam_files" | sed 's/.*cellType\/bam\///' | sed 's/.bam//' | sort -u)
  for cellType in $cellTypes; do
    bam_file=$(echo "$bam_files" | grep "${cellType}.bam$")
    echo "Running analysis for $bam_file..."
    if [ -z "$bam_file" ]; then
      echo "Error: No ${histone} BAM file found in $bam_path/bam";
    else
      fileName="${cellType}"
      echo "Call peak for ${fileName} using sicer with -w ${window_size} -g ${gap_size}";
      ~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t $bam_file -s mm10 -w $window_size -g $gap_size -cpu 20;
      scoreisland="${fileName}-W${window_size}-G${gap_size}.scoreisland";
      if [ ! -f "${scoreisland}" ]; then
        gap_size=5000
        echo "Re-Call peak for ${fileName} using sicer with -w ${window_size} -g ${gap_size}";
        ~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t $bam_file -s mm10 -w $window_size -g $gap_size -cpu 20;
        scoreisland="${fileName}-W${window_size}-G${gap_size}.scoreisland";
      fi

      if [ ! -f "${scoreisland}" ]; then
        echo "Failed to Call peak for ${fileName} using sicer!!!";
        rm *.wig;
      else
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$4,"."}' $scoreisland > "${fileName}-W${window_size}-G${gap_size}.bed";
        mv $scoreisland "${fileName}-W${window_size}-G${gap_size}.scoreisland";
      fi
      rm "${bam_path}/bam/${fileName}.bed"
    fi
  done
done

## Young
for tissue in "${tissues[@]}"; do
  echo "Analysis ${tissue}..."
  bam_path="~/longzhilin/public_html/AgingMouse/PairedTag_new/${tissue}/${histone}/cellType_age"
  if [ ! -d "$bam_path/peak/sicer_Young" ]; then
    mkdir -p $bam_path/peak/sicer_Young
  fi
  cd $bam_path/peak/sicer_Young

  window_size=5000
  gap_size=10000
  bam_files=$(find "$bam_path/bam" -name "*.bam" | grep "_Young.bam$" | grep -vE "^Merge")
  cellTypes=$(echo "$bam_files" | sed 's/.*cellType_age\/bam\///' | sed 's/.bam//' | sort -u)
  for cellType in $cellTypes; do
    bam_file=$(echo "$bam_files" | grep "${cellType}.bam$")
    echo "Running analysis for $bam_file..."
    if [ -z "$bam_file" ]; then
      echo "Error: No ${histone} BAM file found in $bam_path/bam";
    else
      fileName="${cellType}"
      echo "Call peak for ${fileName} using sicer with -w ${window_size} -g ${gap_size}";
      ~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t $bam_file -s mm10 -w $window_size -g $gap_size -cpu 20;
      scoreisland="${fileName}-W${window_size}-G${gap_size}.scoreisland";
      if [ ! -f "${scoreisland}" ]; then
        gap_size=5000
        echo "Re-Call peak for ${fileName} using sicer with -w ${window_size} -g ${gap_size}";
        ~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t $bam_file -s mm10 -w $window_size -g $gap_size -cpu 20;
        scoreisland="${fileName}-W${window_size}-G${gap_size}.scoreisland";
      fi

      if [ ! -f "${scoreisland}" ]; then
        echo "Failed to Call peak for ${fileName} using sicer!!!";
        rm *.wig;
      else
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$4,"."}' $scoreisland > "${fileName}-W${window_size}-G${gap_size}.bed";
        mv $scoreisland "${fileName}-W${window_size}-G${gap_size}.scoreisland";
      fi
      rm "${bam_path}/bam/${fileName}.bed"
    fi
  done
done

## Old
for tissue in "${tissues[@]}"; do
  echo "Analysis ${tissue}..."
  bam_path="~/longzhilin/public_html/AgingMouse/PairedTag_new/${tissue}/${histone}/cellType_age"
  if [ ! -d "$bam_path/peak/sicer_Old" ]; then
    mkdir -p $bam_path/peak/sicer_Old
  fi
  cd $bam_path/peak/sicer_Old

  window_size=5000
  gap_size=10000
  bam_files=$(find "$bam_path/bam" -name "*.bam" | grep "_Old.bam$" | grep -vE "^Merge")
  cellTypes=$(echo "$bam_files" | sed 's/.*cellType_age\/bam\///' | sed 's/.bam//' | sort -u)
  for cellType in $cellTypes; do
    bam_file=$(echo "$bam_files" | grep "${cellType}.bam$")
    echo "Running analysis for $bam_file..."
    if [ -z "$bam_file" ]; then
      echo "Error: No ${histone} BAM file found in $bam_path/bam";
    else
      fileName="${cellType}"
      echo "Call peak for ${fileName} using sicer with -w ${window_size} -g ${gap_size}";
      ~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t $bam_file -s mm10 -w $window_size -g $gap_size -cpu 20;
      scoreisland="${fileName}-W${window_size}-G${gap_size}.scoreisland";
      if [ ! -f "${scoreisland}" ]; then
        gap_size=5000
        echo "Re-Call peak for ${fileName} using sicer with -w ${window_size} -g ${gap_size}";
        ~/longzhilin/Softwares/Epigenomics/SICER2/bin/sicer -t $bam_file -s mm10 -w $window_size -g $gap_size -cpu 20;
        scoreisland="${fileName}-W${window_size}-G${gap_size}.scoreisland";
      fi

      if [ ! -f "${scoreisland}" ]; then
        echo "Failed to Call peak for ${fileName} using sicer!!!";
        rm *.wig;
      else
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$4,"."}' $scoreisland > "${fileName}-W${window_size}-G${gap_size}.bed";
        mv $scoreisland "${fileName}-W${window_size}-G${gap_size}.scoreisland";
      fi
      rm "${bam_path}/bam/${fileName}.bed"
    fi
  done
done
