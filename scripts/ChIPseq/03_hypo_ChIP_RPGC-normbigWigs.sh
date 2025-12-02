#!/bin/bash
#SBATCH --job-name=make_Chip_bigwig
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -o make_bigwig_%j.out
#SBATCH -e make_bigwig_%j.err
#SBATCH --mail-user=yad008@ucsd.edu
#SBATCH --mail-type=END,FAIL

# Create RPGC-normalized bigWig tracks for H3K27ac ChIP-seq BAMs

# ---- paths ----
ROOT="/stg3/data1/yanmiao/yanmiao2/addition_model_Mouse_Mor/hypothamus/Hypothalamus"
CHIP_DIR="${ROOT}/chip"
OUT_DIR="${CHIP_DIR}/bw_norm"

mkdir -p "${OUT_DIR}"

# ---- genome / normalization settings (mm10) ----
# Effective genome size from deepTools docs for mm10 (~2.15e9)
GENOME_SIZE=2150570000
BIN_SIZE=25
SMOOTH_LEN=150
THREADS=8

echo ">>> ROOT      : ${ROOT}"
echo ">>> CHIP_DIR  : ${CHIP_DIR}"
echo ">>> OUT_DIR   : ${OUT_DIR}"
echo ">>> Genome    : mm10 (effective size ${GENOME_SIZE})"
echo

# ---- loop over BAMs ----
# Expect sorted BAMs anywhere under chip/, e.g.:
# chip/MOR-LPS-F/H3K27ac_MOR-LPS-F_10-1.bam
# chip/SAL-LPS-M/H3K27ac_SAL-LPS-M_13-1.bam
#
# Adjust the pattern if your BAM names differ.
find "${CHIP_DIR}" -type f -name "*unique.bam" ! -path "${OUT_DIR}/*" | sort | while read -r BAM; do
    bn=$(basename "${BAM}")
    sample=${bn%.bam}

    # output file name like H3K27ac_MOR-SAL_13-1_M.rpgc.bw
    bw_out="${OUT_DIR}/${sample}.rpgc.bw"

    echo ">>> Processing BAM: ${BAM}"
    echo "    → ${bw_out}"

    # index if missing
    if [ ! -f "${BAM}.bai" ]; then
        echo "    indexing BAM..."
        samtools index -@ "${THREADS}" "${BAM}"
    fi

    # RPGC-normalized bigWig
    bamCoverage \
        --bam "${BAM}" \
        --outFileName "${bw_out}" \
        --outFileFormat bigwig \
        --binSize "${BIN_SIZE}" \
        --smoothLength "${SMOOTH_LEN}" \
        --normalizeUsing RPGC \
        --effectiveGenomeSize "${GENOME_SIZE}" \
        --extendReads 200\
        --centerReads \
        --ignoreDuplicates \
        --minMappingQuality 30 \
        --numberOfProcessors "${THREADS}"

    echo "    done."
    echo
done

echo "✅ Finished generating H3K27ac RPGC bigWigs in: ${OUT_DIR}"
