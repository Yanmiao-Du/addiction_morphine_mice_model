# 03_hypo_PCA_RNA_and_ChIP_Fig1d_e.R
# - Apply QC filters to RNA-seq and ChIP-seq samples
# - Perform PCA on log2-normalized RNA counts and ChIP signal matrix
# - Save separate PCA plots and combined Fig.1d–e panel
# - Optionally export interactive 3D PCA HTML files

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
  library(plotly)
  library(htmlwidgets)
  library(matrixStats)
  library(patchwork)
})

ROOT <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo"
expr_file        <- file.path(ROOT, "combined_by_gene_filtered.tsv")
meta_file        <- file.path(ROOT, "meta_samples.tsv")
chip_signal_file <- file.path(ROOT, "chip_signal_matrix.tsv")
chip_qc_file     <- file.path(ROOT, "QC/chip-QC/ChIP_QC_combined_hypo.tsv")
rna_qc_file      <- file.path(ROOT, "QC/rna-QC/merged_RNA_qc.tsv")

treat_order <- c("SAL-SAL","SAL-LPS","MOR-SAL","MOR-LPS")
okabe_ito <- c("MOR-LPS"="#E69F00","MOR-SAL"="#0072B2",
               "SAL-LPS"="#009E73","SAL-SAL"="#D55E00")

theme_journal <- theme_classic(base_size = 10) +
  theme(
    panel.grid.major.y = element_line(color="grey90", linewidth=.3),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(color="black"),
    axis.text.y = element_text(color="black"),
    strip.background = element_blank(),
    strip.text = element_text(face="bold", size=10),
    plot.title = element_text(face="bold"),
    legend.position = "right",
    plot.margin = margin(6,8,6,8)
  )

POINT_SZ <- 2

# ==============================================================
# ---------------------- 1. RNA section -------------------------
# ==============================================================

# QC thresholds
thr_rna_align <- 70
thr_rna_exon  <- 20
thr_rna_reads <- 5e6

# --- Load expression matrix and metadata ---
expr <- read_tsv(expr_file)
gene_col <- colnames(expr)[1]
expr_mat <- expr %>% column_to_rownames(gene_col) %>% as.matrix()

meta <- read_tsv(meta_file) %>%
  mutate(
    Sample_ID = str_extract(Sample, "\\d+-\\d+"),
    Treatment = str_replace_all(Treatment, "_","-"),
    Treatment = factor(Treatment, levels = treat_order)
  )

# --- RNA QC filter ---
rna_qc_df <- read_tsv(rna_qc_file, show_col_types = FALSE)
rna_qc_df <- rna_qc_df %>%
  mutate(
    AlignmentRate   = dplyr::coalesce(`Uniquely.Mapped....`),
    exon_percentage = dplyr::coalesce(`Exon_Percentage`, NA_real_),
    TotalReads      = dplyr::coalesce(`Total.Reads`),
    UniqueMappedReads = dplyr::coalesce(`Uniquely.Mapped.Reads`)
  )%>%
  mutate(
    AlignmentRate   = as.numeric(AlignmentRate),
    exon_percentage = as.numeric(exon_percentage),
    TotalReads      = as.numeric(TotalReads),
    RNA_pass = !is.na(AlignmentRate) &
               AlignmentRate >= thr_rna_align &
               !is.na(TotalReads) & TotalReads >= thr_rna_reads
  )

qc_pass_rna <- rna_qc_df %>%
  filter(RNA_pass) %>%
  mutate(Sample_ID = str_extract(Sample, "\\d+-\\d+")) %>%
  pull(Sample_ID) %>%
  unique()

cat("✅ RNA QC-passed samples:", length(qc_pass_rna), "\n")

# --- Subset RNA data ---
meta_rna <- meta %>% filter(Modality == "RNA")
colnames(expr_mat) <- gsub("^X", "", colnames(expr_mat))
colnames(expr_mat) <- gsub("\\.", "-", colnames(expr_mat))
expr_mat_rna <- expr_mat[, colnames(expr_mat) %in% qc_pass_rna, drop = FALSE]
meta_rna_sub <- meta_rna %>% filter(Sample_ID %in% qc_pass_rna)

# --- Log2 + PCA ---
X <- log2(expr_mat_rna + 1)
X <- X[rowSums(is.finite(X)) == ncol(X), ]
gene_var <- rowVars(X)
X <- X[gene_var > 0 & is.finite(gene_var), ]

top_n <- min(5000, nrow(X))
if (nrow(X) > top_n) {
  X <- X[order(rowVars(X), decreasing = TRUE)[1:top_n], ]
}

pca_rna <- prcomp(t(X), center = TRUE, scale. = TRUE)
ve_rna <- (pca_rna$sdev^2 / sum(pca_rna$sdev^2)) * 100

pca_df_rna <- as.data.frame(pca_rna$x[, 1:3]) %>%
  rownames_to_column("Sample_ID") %>%
  left_join(meta_rna_sub %>% select(Sample_ID, Sample, Treatment), by = "Sample_ID") %>%
  drop_na(Treatment)

# --- Plot RNA PCA ---
P1_rna <- ggplot(pca_df_rna, aes(PC1, PC2, fill = Treatment)) +
  geom_point(shape = 21, size = POINT_SZ, color = "black", stroke = 0.2) +
  scale_fill_manual(values = okabe_ito, drop = FALSE) +
  labs(title = "RNA-seq PCA",
       x = sprintf("PC1 (%.1f%%)", ve_rna[1]),
       y = sprintf("PC2 (%.1f%%)", ve_rna[2])) +
  theme_journal

ggsave(file.path(ROOT, "Fig_PCA_RNA_PC1_PC2.pdf"), P1_rna,
       width = 4.5, height = 3.8, device = cairo_pdf)

# ==============================================================
# ---------------------- 2. ChIP section ------------------------
# ==============================================================

# ChIP QC thresholds
thr_chip_align <- 70
thr_chip_frip  <- 5
thr_chip_peaks <- 2000

# --- Load signal matrix & QC ---
chip_signal <- read_tsv(chip_signal_file)
chip_qc_df <- read.table(chip_qc_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  rename(`AlignmentRate%` = `OverallAlignment.`, `FRIP%` = `FRiP.`) %>%
  mutate(
    AlignmentRate = as.numeric(`AlignmentRate%`),
    FRIP = as.numeric(`FRIP%`),
    NumPeaks = as.numeric(NumPeaks),
    CHIP_pass = !is.na(AlignmentRate) & !is.na(FRIP) &
                AlignmentRate >= thr_chip_align & FRIP >= thr_chip_frip &
                !is.na(NumPeaks) & NumPeaks >= thr_chip_peaks
  )

qc_pass_chip <- chip_qc_df %>%
  filter(CHIP_pass) %>%
  mutate(Sample_ID = str_extract(Sample, "\\d+-\\d+")) %>%
  pull(Sample_ID) %>%
  unique()

cat("✅ ChIP QC-passed samples:", length(qc_pass_chip), "\n")

# --- Clean ChIP matrix ---
chip_signal <- chip_signal %>%
  drop_na() %>%
  filter(rowSums(across(-Peak), na.rm = TRUE) > 0)

chip_mat <- chip_signal %>%
  column_to_rownames("Peak") %>%
  as.matrix() %>%
  { log2(. + 1) } %>%
  t()

rownames(chip_mat) <- gsub("^X", "", rownames(chip_mat))
rownames(chip_mat) <- gsub("\\.", "-", rownames(chip_mat))
chip_mat <- chip_mat[, apply(chip_mat, 2, sd, na.rm = TRUE) > 0]
chip_mat <- chip_mat[rownames(chip_mat) %in% qc_pass_chip, , drop = FALSE]

# --- Filter metadata ---
chip_meta <- meta %>%
  filter(Modality == "ChIP") %>%
  mutate(Sample_ID = gsub("\\.", "-", Sample_ID)) %>%
  filter(Sample_ID %in% qc_pass_chip)

# --- PCA ---
pca_chip <- prcomp(chip_mat, scale. = TRUE, center = TRUE)
ve_chip <- (pca_chip$sdev^2 / sum(pca_chip$sdev^2)) * 100

pca_df_chip <- as.data.frame(pca_chip$x[, 1:3]) %>%
  tibble::rownames_to_column("Sample_ID") %>%
  left_join(chip_meta, by = "Sample_ID") %>%
  drop_na(Treatment)

# --- Plot ChIP PCA ---
P1_chip <- ggplot(pca_df_chip, aes(PC1, PC2, fill = Treatment)) +
  geom_point(shape = 21, size = POINT_SZ, color = "black", stroke = 0.2) +
  scale_fill_manual(values = okabe_ito, drop = FALSE) +
  labs(title = "ChIP-seq PCA",
       x = sprintf("PC1 (%.1f%%)", ve_chip[1]),
       y = sprintf("PC2 (%.1f%%)", ve_chip[2])) +
  theme_journal

ggsave(file.path(ROOT, "Fig_PCA_CHIP_PC1_PC2.pdf"), P1_chip,
       width = 4.5, height = 3.8, device = cairo_pdf)

# ==============================================================
# --------------- Optional: interactive 3D PCAs -----------------
# ==============================================================

# RNA 3D
p3d_rna <- plot_ly(pca_df_rna,
  x = ~PC1, y = ~PC2, z = ~PC3,
  color = ~Treatment, colors = unname(okabe_ito[levels(pca_df_rna$Treatment)]),
  type = "scatter3d", mode = "markers",
  marker = list(size = 5, line = list(width = 0.5, color = "black")),
  text = ~paste("Sample:", Sample, "<br>Treatment:", Treatment)
)
htmlwidgets::saveWidget(p3d_rna, file.path(ROOT, "RNA_PCA_3D.html"), selfcontained = FALSE)

# ChIP 3D
p3d_chip <- plot_ly(pca_df_chip,
  x = ~PC1, y = ~PC2, z = ~PC3,
  color = ~Treatment, colors = unname(okabe_ito[levels(pca_df_chip$Treatment)]),
  type = "scatter3d", mode = "markers",
  marker = list(size = 5, line = list(width = 0.5, color = "black")),
  text = ~paste("Sample:", Sample, "<br>Treatment:", Treatment)
)
htmlwidgets::saveWidget(p3d_chip, file.path(ROOT, "ChIP_PCA_3D.html"), selfcontained = FALSE)

cat("\n✅ PCA analysis complete. RNA and ChIP figures written to:\n", ROOT, "\n")

library(patchwork)

# Combine RNA (Fig.1e) and ChIP (Fig.1d) panels
composite_pca <- (P1_rna |P1_chip ) +
  plot_annotation(tag_levels = "a")  # labels "a", "b"

# Save combined figure (Fig.1d–e style)
ggsave(file.path(ROOT, "Fig1d_e_PCA_RNA_ChIP.pdf"),
       composite_pca, width = 9, height = 3.8, device = cairo_pdf)

cat("✅ Combined PCA figure saved as Fig1d_e_PCA_RNA_ChIP.pdf\n")
