#' @description split DNA bam for other tissues
#' @Author: Zhilin Long
#' @date 2026-4-24

#### ---- construct the split files ---- ####
source("~/Functions/makeGroupForPseudoBulk.R")
tissue.seurat.final <- readRDS("tissue.seurat.final.rds")
split.files <- lapply(tissue.seurat.final, function(x){
  cat("Analysis", x$tissue[1], "\n")
  meta.data <- x@meta.data
  meta.data$batch <- gsub("-.*", "", meta.data$assay)
  tissue <- unique(meta.data$tissue)
  public_html <- file.path("/storage/zhangyanxiaoLab/longzhilin/public_html/AgingMouse/PairedTag_new", tissue)
  meta.data$cellType <- meta.data$cell_type
  res <- makeGroupForPseudoBulk(meta.data, public_html, "cellType", split.by = "histone")
  meta.data$cellType_age <- paste0(meta.data$cell_type, "_", meta.data$age_group)
  res <- makeGroupForPseudoBulk(meta.data, public_html, "cellType_age", split.by = "histone")
  meta.data$cellType_age_replication <- paste0(meta.data$cell_type, "_", meta.data$age_group, "_", meta.data$sample, "_", meta.data$batch)
  res <- makeGroupForPseudoBulk(meta.data, public_html, "cellType_age_replication", split.by = "histone")
})

#### count RNA & DNA depths
bam_path="~/longzhilin/public_html/AgingMouse/PairedTag_new"
output_txt="${bam_path}/bam_read_counts.txt"

NPROC=24
count_bam() {
  local n
  n=$(samtools idxstats -@ 12 "$1" | awk '{s += $3 + $4} END {print s + 0}')
  printf '%s\t%s\n' "$1" "$n"
}
export -f count_bam

printf 'bam_file\tread_count\n' > "$output_txt"
find "$bam_path" -type f -name "*.bam" -print0 \
  | xargs -0 -P "$NPROC" -n 1 bash -c 'count_bam "$1"' _ \
  >> "$output_txt"

echo "Done. Read counts written to: $output_txt"

#### ---- split DNA ---- ####
merged_DNA_path="~/projects/paired_seq_tag/data/merge"
tissueNames=("BAT" "brainCB" "brainFC" "brainHip" "colon" "heart" "kidney" "liver" "lung" "muscle" "spleen" "skin" "testis") 
histones=("H3K27ac" "H3K27me3" "H3K9me3" "H3K4me1" "H3K4me3" "H3K36me3")
for tissueName in "${tissueNames[@]}"; do
  echo "Processing $tissueName..."
  workPath="~/longzhilin/public_html/AgingMouse/PairedTag_new/${tissueName}"
  cd $workPath
  # histones=$(ls -d H3*/ | sed 's/\///')
  for histone in "${histones[@]}"; do
    cd "$workPath/$histone"
    echo "Running analysis for $workPath/$histone..."

    bam_file=$(find $merged_DNA_path -type f -name "dna.metacell_${histone}.bam" | grep "${tissueName}" | grep -vE "tmp|dna_cutoff_test")
    file_count=$(echo "$bam_file" | wc -l)
    if [ "$file_count" -ge 2 ]; then
      echo "Multiple bam files found for $histone in $merged_DNA_path: $bam_file"
      Merge_DNA_file="$workPath/$histone/${tissueName}.${histone}.merged_DNA.bam"
      samtools merge -@ 12 $Merge_DNA_file $bam_file
      samtools index -@ 12 $Merge_DNA_file
      merged_DNA="$Merge_DNA_file"
    elif [ "$file_count" -eq 0 ]; then
      echo "Error: no bam file found for $histone in $merged_DNA_path"
      exit 1
    else
      echo "single bam file for $histone in $merged_DNA_path: $bam_file"
      merged_DNA="${bam_file}"
    fi

    if [[ "${histone}" =~ (H3K27me3|H3K9me3|H3K36me3) ]]; then
      smoothBinSize=5000
      smoothLength=50000
    else
      smoothBinSize=50
      smoothLength=200
    fi

    split_files=("cellType_age" "cellType_age_replication" "cellType")
    for group in "${split_files[@]}"; do
      file="$workPath/$histone/$group/${group}.txt"
      if [ -e $file ]; then
        outDir=$(dirname $file);
        cd $outDir

        echo "Running with $file ...";
        ## split bam file
        if [ -d "$outDir/bam" ]; then
          rm -rf $outDir/bam;
          rm -rf $outDir/bw;
        fi
        bash ~/Functions/splitBam.sh $merged_DNA $file $histone $smoothBinSize $smoothLength;
      else
        echo "$file not exist!"
      fi
    done

    if [ -f "$Merge_DNA_file" ]; then
      rm "$Merge_DNA_file" "${Merge_DNA_file}.bai"
    fi
  done
done
