#!/bin/bash
#SBATCH --job-name=make_RNA_bigwig
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -o make_bigwig_%j.out
#SBATCH -e make_bigwig_%j.err
#SBATCH --mail-user=yad008@ucsd.edu
#SBATCH --mail-type=END,FAIL

# Create CPM-normalized bigWig tracks for RNA-seq BAMs

ROOT="/stg3/data1/yanmiao/yanmiao2/addition_model_Mouse_Mor/hypothamus/Hypothalamus"
RNA_DIR="${ROOT}/RNA_BAM"
OUT_DIR="${ROOT}/rna/bw"

mkdir -p "${OUT_DIR}"

BIN_SIZE=25
SMOOTH_LEN=75
THREADS=8

echo ">>> ROOT     : ${ROOT}"
echo ">>> RNA_DIR  : ${RNA_DIR}"
echo ">>> OUT_DIR  : ${OUT_DIR}"
echo

# Expect sorted BAMs in rna/, e.g.:
# rna/SAL-SAL_10-1_F.bam
# rna/MOR-SAL_13-1_M.bam
#
# Adjust pattern if needed.
find "${RNA_DIR}" -type f -name "*.bam" ! -path "${OUT_DIR}/*" | sort | while read -r BAM; do
    bn=$(basename "${BAM}")
    sample=${bn%.bam}

    # output file like RNA_MOR-SAL_13-1_M.bw
    bw_out="${OUT_DIR}/RNA_${sample}.bw"

    echo ">>> Processing BAM: ${BAM}"
    echo "    → ${bw_out}"

    # index if missing
    if [ ! -f "${BAM}.bai" ]; then
        echo "    indexing BAM..."
        samtools index -@ "${THREADS}" "${BAM}"
    fi

    # CPM-normalized RNA coverage
    bamCoverage \
        --bam "${BAM}" \
        --outFileName "${bw_out}" \
        --outFileFormat bigwig \
        --binSize "${BIN_SIZE}" \
        --smoothLength "${SMOOTH_LEN}" \
        --normalizeUsing CPM \
        --extendReads \
        --centerReads \
        --ignoreDuplicates \
        --minMappingQuality 30 \
        --numberOfProcessors "${THREADS}"

    echo "    done."
    echo
done

echo "✅ Finished generating RNA CPM bigWigs in: ${OUT_DIR}"