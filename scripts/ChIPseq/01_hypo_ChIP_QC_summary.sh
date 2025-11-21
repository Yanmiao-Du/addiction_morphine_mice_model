#!/bin/bash

# Define base directory
BASE_DIR="/stg3/data1/yanmiao/yanmiao2/addition_model_Mouse_Mor/hippo/chip"
PEAK_DIR="/stg3/data1/yanmiao/yanmiao2/addition_model_Mouse_Mor/hippo/chip/narrowpeak_file"

# Define output files
output_file="$BASE_DIR/combined_qc.tsv"
frip_output="$BASE_DIR/frip_scores.tsv"

echo -e "Sample\tUniqueMappedReads\tOverallAlignment%\tNumPeaks" > "$output_file"

# Step 1: Extract Alignment Metrics from Bowtie2 logs
declare -A unique_reads_map
declare -A overall_alignment_map


# Use a loop without piping `find` into `while read`
while IFS= read -r log; do
    sample=$(basename "$log" .log)

    # Extract metrics
    unique_reads=$(grep "aligned concordantly exactly 1 time" "$log" | awk '{print $1}')
    total_reads=$(grep "reads; of these:" "$log" | awk '{print $1}')
    overall_alignment=$(grep "overall alignment rate" "$log" | awk '{print $1}' | tr -d '%')

    # Assign values to associative arrays
    unique_reads_map["$sample"]="$unique_reads"
    overall_alignment_map["$sample"]="$overall_alignment"
done < <(find "$BASE_DIR" -type f -path "*/Aligned_SAM/*.log")

# Print results
for sample in "${!unique_reads_map[@]}"; do
    echo "Sample: $sample | Unique Reads: ${unique_reads_map[$sample]} | Alignment Rate: ${overall_alignment_map[$sample]}%"
done

# Step 2: Extract Peak Counts from narrowPeak files
declare -A num_peaks_map

while IFS= read -r peak_file; do
    sample=$(basename "$peak_file" _peaks.narrowPeak)
    num_peaks=$(wc -l < "$peak_file")
    
    num_peaks_map["$sample"]="$num_peaks"
done < <(find "$PEAK_DIR" -type f -name "*_peaks.narrowPeak")


# Step 3: Merge Data into Combined File
for sample in "${!unique_reads_map[@]}"; do
    unique_reads="${unique_reads_map[$sample]:-NA}"
    overall_alignment="${overall_alignment_map[$sample]:-NA}"
    num_peaks="${num_peaks_map[$sample]:-NA}"
    
    echo -e "$sample\t$unique_reads\t$overall_alignment\t$num_peaks" >> "$output_file"
done

echo "Combined QC metrics saved in $output_file"
