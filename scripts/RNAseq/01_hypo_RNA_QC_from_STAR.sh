#!/bin/bash
# 01_hypo_RNA_QC_from_STAR.sh
# Parse STAR Log.final.out files to extract uniquely mapped read counts
# and percentages for hypothalamus RNA-seq samples.

BASE_DIR="/stg3/data1/yanmiao/yanmiao2/addition_model_Mouse_Mor/hippo/rna/rna_p2/RNA_BAM"
output_file="$BASE_DIR/rna_qc.tsv"

echo -e "Sample\tUniqueMappedReads\tUniqueMappedReads%" > "$output_file"

find "$BASE_DIR" -type f -name "*_Log.final.out" | while read log; do
    sample=$(basename "$log" "_Log.final.out")
    unique_mapped_num=$(grep "Uniquely mapped reads number" "$log" | awk '{print $NF}')
    unique_mapped_rate=$(grep "Uniquely mapped reads %" "$log" | awk '{print $NF}')
    echo -e "$sample\t$unique_mapped_num\t$unique_mapped_rate" >> "$output_file"
done

echo "RNA QC metrics saved in $output_file"
