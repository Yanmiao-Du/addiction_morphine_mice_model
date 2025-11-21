# addiction_morphine_mice_model

Integrative multi-omic analysis of a three-trimester mouse model of perinatal morphine exposure. This repository contains the analysis code and supporting files for combining H3K27ac ChIP-seq and bulk RNA-seq (plus cytokine and physiological data) across multiple brain regions. Current scripts concentrate on the hypothalamus analyses used in the manuscript.

---

## Table of contents

- Project overview
- Repository layout
- Scripts (what to run)
- Metadata (required files)
- Reproducing analyses
- Contact & citation

---

## Project overview

This project investigates long-term immune and metabolic dysregulation after perinatal morphine exposure using complementary genomic and phenotypic assays. Key data modalities:

- H3K27ac ChIP-seq (epigenomic regulation)
- Bulk RNA-seq (transcriptional changes)
- Cytokine profiling and physiological phenotypes

Analyses include QC, PCA, differential expression (DESeq2), ChIP QC, and integrative figure generation for manuscript figures (e.g., PCA and DE results for hypothalamus samples).

## Repository layout

Top-level layout (relevant folders):

```
scripts/
  RNAseq/        # RNA-seq analysis scripts (STAR QC, DESeq2, PCA, plotting)
  ChIPseq/       # ChIP-seq QC and summary scripts
  utilities/     # small helper scripts used by analyses
data/
  metadata/      # sample manifests and design matrices (required)
docs/
  methods/       # method writeups and notes
figs/            # figure outputs used for manuscript
results/         # intermediate analysis outputs and summary tables
README.md        # this file
```

Example files you will find:

- `scripts/RNAseq/01_hypo_RNA_QC_from_STAR.sh` — collect STAR QC metrics and summarize sample-level QC
- `scripts/RNAseq/02_hypo_RNA_DESeq2_MOR-SAL_vs_SAL-SAL_Fig2A.R` — DE analysis and figure generation for Fig.2A
- `scripts/RNAseq/03_hypo_PCA_RNA_and_ChIP_Fig1d_e.R` — PCA plots combining RNA and ChIP
- `scripts/ChIPseq/01_hypo_ChIP_QC_summary.sh` — ChIP-seq QC aggregation

## Metadata

All structured sample metadata for the hypothalamus analyses lives under `data/metadata/`. The primary metadata file contains:

- Sample names and sequencing IDs (RNA & ChIP)
- Treatment groups (`SAL-SAL`, `SAL-LPS`, `MOR-SAL`, `MOR-LPS`)
- Sex, batch, and modality annotations
- Mappings between raw sequencing files and the analysis sample labels

This metadata is required for:

- Building design matrices for DESeq2
- PCA and multi-omic sample grouping (Fig.1d–e)
- Matching STAR logs, MACS2 outputs, and count matrices to sample labels

## Reproducing analyses

Prerequisites:

- R (version used for manuscript; see script headers)
- STAR, samtools, MACS2 for alignment and peak calling
- Standard Unix tools (bash, awk, sed)

Typical workflow (high level):

1. Prepare metadata in `data/metadata/` (ensure sample names and IDs match raw files).
2. Run RNA-seq QC and alignment summary: `scripts/RNAseq/01_hypo_RNA_QC_from_STAR.sh`.
3. Perform differential expression: `scripts/RNAseq/02_hypo_RNA_DESeq2_MOR-SAL_vs_SAL-SAL_Fig2A.R`.
4. Generate PCA figures and combined RNA/ChIP plots: `scripts/RNAseq/03_hypo_PCA_RNA_and_ChIP_Fig1d_e.R`.
5. Aggregate ChIP QC: `scripts/ChIPseq/01_hypo_ChIP_QC_summary.sh`.

Each script contains header comments describing inputs, expected file paths, and required R or shell packages. Inspect those headers for exact runtime commands and version notes.

## Outputs

- Figures used in the manuscript are stored under `results/figures/` and `figs/` (PDFs and summary plots).
- Intermediate tables and count matrices appear under `results/` and may be referenced by the R scripts.

## Contact & citation

If you use this repository or the associated data, please cite the manuscript (when available) and contact the maintainer for data access or questions.

Maintainer: Yanmiao Du
Repository: `addiction_morphine_mice_model`
