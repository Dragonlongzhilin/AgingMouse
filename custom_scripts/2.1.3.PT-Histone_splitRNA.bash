#' @description split RNA bam for other tissues
#' @Author: Zhilin Long
#' @date 2026-4-24

#### ---- split RNA ---- ####
merged_RNA_path="~/projects/paired_seq_tag/data/merge"
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

    bam_file=$(find $merged_RNA_path -type f -name "RNA_final_sorted.bam" | grep "${histone}" | grep "${tissueName}" | grep -vE "tmp|dna_cutoff_test")
    file_count=$(echo "$bam_file" | wc -l)
    if [ "$file_count" -ge 2 ]; then
      echo "Multiple bam files found for $histone in $merged_RNA_path: $bam_file"
      Merge_RNA_file="$workPath/$histone/${tissueName}.${histone}.merged_RNA.bam"
      samtools merge -@ 12 $Merge_RNA_file $bam_file
      samtools index -@ 12 $Merge_RNA_file
      merged_RNA="$Merge_RNA_file"
    elif [ "$file_count" -eq 0 ]; then
      echo "Error: no bam file found for $histone in $merged_RNA_path"
      exit 1
    else
      echo "single bam file for $histone in $merged_RNA_path: $bam_file"
      merged_RNA="${bam_file}"
    fi

    smoothBinSize=50
    smoothLength=200

    split_files=("cellType_age")
    for group in "${split_files[@]}"; do
      file="$workPath/$histone/$group/${group}.txt"
      if [ -e $file ]; then
        outDir=$(dirname $file);
        mkdir -p $outDir/RNA
        cd $outDir/RNA

        echo "Running with $file ...";
        ## split bam file
        if [ -d "$outDir/RNA/bam" ]; then
          rm -rf $outDir/RNA/bam;
          rm -rf $outDir/RNA/bw;
        fi
        bash ~/Functions/splitBam.sh $merged_RNA $file $histone $smoothBinSize $smoothLength;
      else
        echo "$file not exist!"
      fi
    done

    if [ -f "$Merge_RNA_file" ]; then
      rm "$Merge_RNA_file" "${Merge_RNA_file}.bai"
    fi
  done
done

#### ---- re-smooth bw ---- ####
bam_files=$(find ~longzhilin/public_html/AgingMouse/PairedTag_new -type f -name "*.bam" | grep -E "cellType_age/RNA")
for bam_file in $bam_files; do
  outDir=$(dirname $(dirname "$bam_file"))
  outDir="${outDir}/bw"
  cd $outDir
  fileName=$(basename "$bam_file" | sed 's/.bam/.bw/')
  bamCoverage -p 16 --binSize 50 --smoothLength 200 -b "$bam_file" -o $fileName --outFileFormat bigwig --normalizeUsing RPKM
done
