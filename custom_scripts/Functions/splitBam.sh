#!/bin/bash
#' @description make igv session
#' @author Zhilin Long
#' @date 2025-7-3

# annother version for pseudo bulk

# sample_barcode.txt
# '''
# barcode sample
# 01      H3K27ac
# 02      H3K27ac
# 03      H3K27ac
# 04      H3K27ac
# 05      H3K27ac
# 06      H3K27ac
# 07      H3K9me3
# 08      H3K9me3
# 09      H3K9me3
# 10      H3K9me3
# 11      H3K9me3
# 12      H3K9me3
# '''

input_bam=$1  # merge_DNA.bam file
sample_barcode=$2  # sample_barcode.txt
outPrefix=$3  # the prefix of outfile
binSize=$4 # Default: 50
smoothLength=$5
# if the --binSize is set to 20 and the --smoothLength is set to 60, then, for each bin, the average of the bin and its left and right neighbors is considered. default: None

echo "work dir: $(pwd)"

if [ -d "./bam" ]; then
	rm -rf ./bam
fi
if [ -d "./bw" ]; then
	rm -rf ./bw
fi
mkdir -p ./bam
mkdir -p ./bw

## extract sample barcode from bam file to pseudo bulk
python /storage/zhangyanxiaoLab/longzhilin/ProjectCode/AgingProject/extract_sample_barcode_from_bam.py \
	--bam $input_bam \
	--statH $sample_barcode \
	--outPrefix $outPrefix

mv *.bam ./bam
for file in ./bam/*.bam; do
  mv -- "$file" "${file// /_}"
done

for BAM_FILE in ./bam/*.bam; do
	echo "Processing $BAM_FILE"
	samtools index "$BAM_FILE"

	bw_file=$(basename "$BAM_FILE" | sed 's/.bam/.bw/g')
	# bigWig
	bamCoverage -p 16 --binSize $binSize --smoothLength $smoothLength -b "$BAM_FILE" -o "./bw/$bw_file" --outFileFormat bigwig --normalizeUsing RPKM
done
