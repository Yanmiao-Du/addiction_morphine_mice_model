# 02_hypo_RNA_DESeq2_MOR-SAL_vs_SAL-SAL_Fig2A.R
# - Filter hypothalamus RNA-seq samples by QC
# - Build DESeq2 object with Group_noSex design
# - Run all pairwise contrasts
# - Save DESeq2 results tables
# - Generate Fig.2A bar summaries for Enhancer / Promoter / DEG counts

library(openxlsx)
library(dplyr)
library(readr)
library(stringr)
library(DESeq2)
library(tidyverse)
library(ggtext)
library(patchwork)

# Set this to your local hypo analysis directory
ROOT <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo"

# --- Load RNA QC table (from Fig.1 code) ---
rna_qc <- read.table("/Users/yanmiaodu/Downloads/MOR_SAL/hypo/QC/rna-QC/merged_RNA_qc.tsv",
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# unify column names if needed
rna_qc <- rna_qc %>%
  mutate(
    AlignmentRate   = dplyr::coalesce(`Uniquely.Mapped....`),
    exon_percentage = dplyr::coalesce(`Exon_Percentage`),
    TotalReads      = dplyr::coalesce(`Total.Reads`)
  ) %>%
  mutate(
    AlignmentRate   = as.numeric(AlignmentRate),
    exon_percentage = as.numeric(exon_percentage),
    TotalReads      = as.numeric(TotalReads),
    RNA_pass = !is.na(AlignmentRate) & AlignmentRate >= 70 &
               !is.na(exon_percentage) & exon_percentage >= 20 &
               !is.na(TotalReads) & TotalReads >= 5e6
  )

# --- Identify pass/fail samples ---
qc_pass_samples <- rna_qc %>%
  filter(RNA_pass) %>%
  mutate(Sample_ID = stringr::str_extract(Sample, "\\d+-\\d+")) %>%
  pull(Sample_ID) %>% unique()

cat("✅ RNA samples passing QC:", length(qc_pass_samples), "\n")
print(qc_pass_samples)

# Keep only QC-passing samples in metadata and expression matrix
meta_path = "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/RNA_meta.xlsx"
meta = read.xlsx(meta_path)
meta_filtered <- meta %>% filter(RNA %in% qc_pass_samples)
expression_file_path = "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/combined_by_gene_filtered.tsv"
expression_data <- read.delim(expression_file_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(expression_data) <- gsub("^X", "", colnames(expression_data))
colnames(expression_data) <- gsub("\\.", "-", colnames(expression_data))

expr_filtered <- expression_data[, c("Gene", meta_filtered$RNA)]

# keep only QC-pass samples in both expression and metadata
valid_rna <- intersect(qc_pass, meta$RNA)

meta_filtered <- meta %>% filter(RNA %in% valid_rna)
expr_filtered <- expression_data[, c("Gene", valid_rna)]

# confirm that all columns match
stopifnot(all(valid_rna %in% colnames(expr_filtered)))

# rename expression columns to sample names in meta
matched_idx <- match(colnames(expr_filtered)[-1], meta_filtered$RNA)
new_colnames <- meta_filtered$name[matched_idx]
colnames(expr_filtered)[-1] <- new_colnames

# prepare DESeq2 input
expr_fixed <- expr_filtered
rownames(expr_fixed) <- expr_fixed$Gene
expr_fixed <- expr_fixed[, -1]

meta_filtered <- meta_filtered[match(colnames(expr_fixed), meta_filtered$name), ]
rownames(meta_filtered) <- meta_filtered$name
meta_filtered <- meta_filtered %>%
  mutate(
    Group_noSex = str_remove(Group, "-[FM]$"),   # remove trailing "-F" or "-M"
    Group_noSex = factor(Group_noSex,
                         levels = c("SAL-SAL", "SAL-LPS", "MOR-SAL", "MOR-LPS"))  # optional fixed order
  )

levels(meta_filtered$Group_noSex) <- gsub("-", "_", levels(meta_filtered$Group_noSex))

# confirm
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = expr_fixed,
  colData   = meta_filtered,
  design    = ~ Group_noSex
)

# filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# run DESeq2
dds <- DESeq(dds)

saveRDS(dds, "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/hypo_deg_qcpass_dds.rds")
dds <- readRDS("/Users/yanmiaodu/Downloads/MOR_SAL/hypo/hypo_deg_qcpass_dds.rds")
all_groups <- levels(dds$Group_noSex)
pairwise_comparisons <- combn(all_groups, 2, simplify = FALSE)

output_dir <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/DEGs_treatment_only"
dir.create(output_dir, showWarnings = FALSE)


for (pair in pairwise_comparisons) {
  g1 <- pair[1]; g2 <- pair[2]
  comp <- paste0(g1, "_vs_", g2)
  cat("Running:", comp, "\n")
  
  # Run DESeq2 contrast
  res <- results(dds, contrast = c("Group_noSex", g1, g2))
  
  # Preserve gene names from rownames
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column(var = "gene_name") %>%
    arrange(padj)
  
  # Write full and significant tables
  readr::write_tsv(res_df, file.path(output_dir, paste0(comp, "_ALL.tsv")))
  readr::write_tsv(subset(res_df, padj < 0.05), file.path(output_dir, paste0(comp, "_sig.tsv")))
}


