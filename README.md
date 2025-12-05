# addiction_morphine_mice_model

Analysis code for a three-trimester mouse model of perinatal morphine exposure.  
The repository contains scripts used to process and integrate H3K27ac ChIP-seq, bulk RNA-seq, and related phenotyping data, with a current focus on the hypothalamus analyses for the manuscript.

---

## Repository contents

- `scripts/`
  - `RNAseq/` – RNA-seq QC, DESeq2, PCA, and plotting scripts.
  - `ChIPseq/` – ChIP-seq QC and summary scripts.
  - `utilities/` – small helper scripts sourced by the main analyses.
- `figs/` – example figure panels / exported plots.
- `docs/methods/` – methods notes and text fragments used in the manuscript.
- (Local, not committed) `data/metadata/` and `results/` – sample sheets, count matrices, and intermediate outputs expected by the scripts.

---

## How to run the main analyses

1. **Clone this repo** and make sure the paths in the script headers match your local project layout.
2. **Prepare metadata** in `data/metadata/`  
   (sample IDs, treatment group: `SAL-SAL`, `SAL-LPS`, `MOR-SAL`, `MOR-LPS`, sex, batch, modality).
3. **RNA-seq pipeline (hypothalamus):**
   - QC + STAR summary: `scripts/RNAseq/01_hypo_RNA_QC_from_STAR.sh`
   - Differential expression (e.g. MOR-SAL vs SAL-SAL):  
     `scripts/RNAseq/02_hypo_RNA_DESeq2_MOR-SAL_vs_SAL-SAL_Fig2A.R`
   - PCA and integrated RNA/ChIP plots (Fig. 1d–e):  
     `scripts/RNAseq/03_hypo_PCA_RNA_and_ChIP_Fig1d_e.R`
4. **ChIP-seq QC:**
   - `scripts/ChIPseq/01_hypo_ChIP_QC_summary.sh`

Each script has a header describing required inputs, software versions, and example command lines. Edit those paths if you are not using the original UCSD cluster layout.

---

## Requirements

- R (version as noted in script headers) with tidyverse, DESeq2, DiffBind, etc.
- UNIX environment with `bash`, `awk`, `sed`
- External tools: STAR, samtools, MACS2 for alignment and peak calling.

---

## Contact / citation

Please cite the associated manuscript (when available) if you use this code, and feel free to reach out with questions or bug reports.

Maintainer: **Yanmiao Du**  
Repository: `addiction_morphine_mice_model`
