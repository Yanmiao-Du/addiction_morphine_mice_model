#!/usr/bin/env Rscript

# Hypothalamus MOR–SAL vs SAL–SAL
# Figure 2C (volcano), Figure 2E (Venn),
# Figure 2F + Supplementary GO panels, and immune GO.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(ggtext)
  library(patchwork)
  library(VennDiagram)
  library(grid)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggalluvial)
  library(ggplot2)
  library(stringr)
})

out_dir <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo"

## ---------- Input ----------
deg_file  <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/DEGs_treatment_only/SAL_SAL_vs_MOR_SAL_ALL.tsv"
chip_file <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/Diff_analysis/Annotated_MORSAL_vs_SALSAL.tsv"
diffbind_rda   <- file.path(ROOT, "hypo_diffbind_DESeq2_allContrasts.RData")

## ---------- Parameters ----------
fdr_cut <- 0.05
lfc_cut <- 0.25

col_deg       <- "#009E73"
col_enhancer  <- "#E69F00"
col_promoter  <- "#56B4E9"
col_ns        <- "grey80"

theme_pub <- theme_classic(base_size = 11, base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 9,  color = "black"),
    axis.title  = element_text(face = "bold"),
    plot.title  = element_text(size = 13, face = "bold", hjust = 0.5),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    plot.margin        = margin(6, 6, 6, 6)
  )

# ===============================================================
# 1. Load data (shared across all Fig2 panels)
# ===============================================================

## --- DEGs: SAL_SAL_vs_MOR_SAL_ALL.tsv ---
deg_raw <- read_tsv(deg_file, show_col_types = FALSE)

deg_tbl <- deg_raw %>%
  dplyr::rename(FDR = padj) %>%
  mutate(
    ## flip so + = MOR–SAL / SAL–SAL
    log2FC = -log2FoldChange,
    log2FC = ifelse(is.na(log2FC), 0, log2FC),
    sig = if_else(FDR < fdr_cut & abs(log2FC) > lfc_cut, "Sig", "NS")
  )

## --- ChIP: annotated differential peaks (MORSAL vs SALSAL) ---
load(diffbind_rda)   # loads dbObj_clean and others

## contrast 2 was MOR-LPS vs SAL-LPS in your earlier code
mor_sal_vs_sal_sal_all <- DiffBind::dba.report(
  dbObj_clean,
  contrast = 5,
  th = 1  # th=1 -> all peaks with FDR
)

chip_all <- as.data.frame(mor_sal_vs_sal_sal_all) %>%
  as_tibble() %>%
  mutate(
    seqnames = as.character(seqnames),
    start = start,
    end = end
  )

## --- Annotation file (only FDR<0.05) for region + SYMBOL ---
chip_anno <- read_tsv(chip_file, show_col_types = FALSE)

chip_anno_min <- chip_anno %>%
  dplyr::select(seqnames, start, end, annotation, SYMBOL)

## join annotation onto full DiffBind table
chip_tbl <- chip_all %>%
  left_join(chip_anno_min,
            by = c("seqnames", "start", "end")) %>%
  mutate(
    region = case_when(
      !is.na(annotation) &
        str_detect(annotation, regex("promoter", ignore_case = TRUE)) ~ "Promoter",
      !is.na(annotation)                                             ~ "Enhancer",
      TRUE                                                           ~ "Unannotated"
    ),
    sig = case_when(
      FDR < fdr_cut & abs(Fold) > lfc_cut & region == "Promoter" ~ "Promoter",
      FDR < fdr_cut & abs(Fold) > lfc_cut & region == "Enhancer" ~ "Enhancer",
      TRUE                                                       ~ "NS"
    )
  )

## Keep a convenience object of FDR-significant, annotated peaks (for bars/GO)
chip_sig <- chip_tbl %>%
  filter(FDR < fdr_cut & region %in% c("Enhancer", "Promoter"))

# ===============================================================
# Figure 2C — Volcano plots for differential peaks & DEGs
# ===============================================================

## ---- choose top DEGs to label ----
top_deg <- deg_tbl %>%
  filter(sig == "Sig") %>%
  arrange(FDR, desc(abs(log2FC))) %>%
  slice_head(n = 10)

p_deg <- ggplot(deg_tbl, aes(x = log2FC, y = -log10(FDR))) +
  geom_point(aes(color = sig), size = 1.6, alpha = 0.8) +
  scale_color_manual(values = c("Sig" = col_deg, "NS" = col_ns)) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(fdr_cut),
             linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = top_deg,
    aes(label = gene_name),
    size = 3, color = "black", max.overlaps = 30
  ) +
  labs(
    title = "RNA-seq (DEGs)",
    x = expression(log[2]~fold~change~"(MOR-SAL / SAL-SAL)"),
    y = expression(-log[10]~adjusted~italic(P))
  ) +
  theme_pub +
  theme(legend.position = "none")

## ---- choose top ChIP peaks to label ----
top_chip <- chip_tbl %>%
  filter(sig != "NS") %>%
  arrange(FDR, desc(abs(Fold))) %>%
  slice_head(n = 10)

chip_plot_tbl <- chip_tbl %>%
  filter(!is.na(FDR), FDR < 0.999)

p_chip <- ggplot(chip_plot_tbl, aes(x = Fold, y = -log10(FDR))) +
  geom_point(aes(color = sig), size = 1.8, alpha = 0.8) +
  scale_color_manual(values = c(
    "Enhancer" = col_enhancer,
    "Promoter" = col_promoter,
    "NS"       = col_ns
  )) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut),
             linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(fdr_cut),
             linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = top_chip,
    aes(label = SYMBOL),
    size = 3,
    max.overlaps = 15
  ) +
  labs(
    title = "H3K27ac differential peaks",
    x = expression(log[2]~fold~change~"(MOR-SAL / SAL-SAL)"),
    y = expression(-log[10]~adjusted~italic(P))
  ) +
  theme_pub +
  theme(
    legend.position = "top",
    legend.title    = element_blank()
  )

p_combined <- p_chip | p_deg +
  plot_annotation(
    title = "MOR–SAL vs SAL–SAL",
    subtitle = "Differential peaks (Enhancer/Promoter) and DEGs",
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, color = "grey25", hjust = 0.5)
    )
  )

ggsave(file.path(out_dir, "Fig2C_volcano_MOR-SAL_vs_SAL-SAL.pdf"),
       p_combined, width = 7, height = 3.5, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig2C_volcano_MOR-SAL_vs_SAL-SAL.png"),
       p_combined, width = 7, height = 3.5, dpi = 600)

# ===============================================================
# Figure 2E — Numeric-only Venn diagram (FDR < 0.05 only)
# ===============================================================

deg_genes <- deg_tbl %>%
  filter(FDR < fdr_cut) %>%
  pull(gene_name) %>%
  unique()

enhancer_genes <- chip_sig %>%
  filter(region == "Enhancer", !is.na(SYMBOL), SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique()

promoter_genes <- chip_sig %>%
  filter(region == "Promoter", !is.na(SYMBOL), SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique()

venn_list <- list(
  DEGs     = deg_genes,
  Enhancer = enhancer_genes,
  Promoter = promoter_genes
)

venn_plot <- draw.triple.venn(
  area1 = length(venn_list$DEGs),
  area2 = length(venn_list$Enhancer),
  area3 = length(venn_list$Promoter),
  n12   = length(intersect(venn_list$DEGs, venn_list$Enhancer)),
  n23   = length(intersect(venn_list$Enhancer, venn_list$Promoter)),
  n13   = length(intersect(venn_list$DEGs, venn_list$Promoter)),
  n123  = length(Reduce(intersect, venn_list)),
  category = c("DEGs", "Enhancer", "Promoter"),
  fill = c("#00A087", "#E69F00", "#56B4E9"),
  alpha = rep(0.7, 3),
  cat.cex = 1.3,
  cex = 1.3,
  cat.fontface = "bold",
  cat.col = "black",
  fontfamily = "sans",
  cat.fontfamily = "sans",
  cat.pos  = c(-10, 10, 0),
  cat.dist = c(0.05, 0.05, 0.05),
  margin   = 0.08,
  scaled   = FALSE,
  euler.d  = FALSE
)

pdf(file.path(out_dir, "Fig2E_Venn_numeric_MOR-SAL_vs_SAL-SAL_centered.pdf"),
    width = 4.5, height = 4.5, family = "Helvetica")
grid.draw(venn_plot)
dev.off()

png(file.path(out_dir, "Fig2E_Venn_numeric_MOR-SAL_vs_SAL-SAL_centered.png"),
    width = 1200, height = 1200, res = 300, type = "cairo")
grid.draw(venn_plot)
dev.off()

# ===============================================================
# GO enrichment (Fig2F + Suppl S2B/S2C)
# ===============================================================

## --- define gene sets using SAME objects as above ----
deg_all_sig <- deg_tbl %>%
  filter(FDR < fdr_cut & abs(log2FC) > lfc_cut) %>%
  pull(gene_name) %>%
  unique()

deg_up <- deg_tbl %>%
  filter(FDR < fdr_cut & log2FC >  lfc_cut) %>%
  pull(gene_name) %>%
  unique()

deg_down <- deg_tbl %>%
  filter(FDR < fdr_cut & log2FC < -lfc_cut) %>%
  pull(gene_name) %>%
  unique()

enh_genes <- chip_tbl %>%
  filter(FDR < fdr_cut & region == "Enhancer",
         !is.na(SYMBOL), SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique()

prom_genes <- chip_tbl %>%
  filter(FDR < fdr_cut & region == "Promoter",
         !is.na(SYMBOL), SYMBOL != "") %>%
  pull(SYMBOL) %>%
  unique()

diff_genes <- unique(c(enh_genes, prom_genes))

# --- enrichGO (same as before, just fed with updated gene sets) ---
ego_deg <- enrichGO(
  gene          = deg_all_sig,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
write_tsv(as.data.frame(ego_deg),
          "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_DEG_all_GO_BP_enrichment.tsv")

ego_deg_up <- enrichGO(
  gene          = deg_up,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
write_tsv(as.data.frame(ego_deg_up),
          "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_DEG_up_GO_BP_enrichment.tsv")

ego_deg_down <- enrichGO(
  gene          = deg_down,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
write_tsv(as.data.frame(ego_deg_down),
          "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_DEG_down_GO_BP_enrichment.tsv")

ego_enh <- enrichGO(
  gene = enh_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05
)
write_tsv(as.data.frame(ego_enh),
          "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_enhancer_peakall_GO_BP_enrichment.tsv")

ego_pro <- enrichGO(
  gene = prom_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05
)
write_tsv(as.data.frame(ego_pro),
          "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_promoter_peakall_GO_BP_enrichment.tsv")

ego_diff <- enrichGO(
  gene = diff_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05
)
write_tsv(as.data.frame(ego_diff),
          "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_diff_peakall_GO_BP_enrichment.tsv")

## --- updated GO dotplot helper (same as Fig3F) ---
plot_go_clean <- function(ego, title, top_n = 10, wrap_width = 40) {

  df <- as.data.frame(ego) %>%
    dplyr::mutate(
      p.adjust  = as.numeric(p.adjust),
      Count     = as.numeric(Count),
      GeneRatio = as.character(GeneRatio)
    ) %>%
    dplyr::filter(!is.na(p.adjust), p.adjust > 0) %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      Description_wrapped = stringr::str_wrap(Description, wrap_width),
      GO = factor(Description_wrapped, levels = rev(unique(Description_wrapped))),
      GeneRatio_num = sapply(strsplit(GeneRatio, "/"), function(x) {
        as.numeric(x[1]) / as.numeric(x[2])
      }),
      negLog10FDR = -log10(p.adjust)
    )

  # x-axis limits: start a bit smaller than the min
  xmin <- min(df$GeneRatio_num) * 0.9
  xmax <- max(df$GeneRatio_num) * 1.05

  ggplot(df, aes(x = GeneRatio_num, y = GO,
                 size = Count, color = negLog10FDR)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#56B4E9", high = "#E64B35") +
    scale_x_continuous(
      name   = "Gene Ratio",
      limits = c(xmin, xmax),
      breaks = scales::pretty_breaks(4),
      expand = c(0, 0)
    ) +
    labs(
      title = title,
      y     = NULL,
      color = expression(-log[10]~FDR),
      size  = "Count"
    ) +
    theme_bw(base_size = 11, base_family = "Helvetica") +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text.y  = element_text(size = 9,  color = "black"),
      axis.text.x  = element_text(size = 8.5, color = "black"),
      axis.title   = element_text(face = "bold"),
      plot.title   = element_text(face = "bold", size = 12, hjust = 0.5),
      legend.position    = "right",
      legend.key.height  = unit(0.5, "lines"),
      plot.margin        = margin(5.5, 10, 5.5, 5.5)
    )
}


p_deg  <- plot_go_clean(ego_deg,  "RNA-seq (DEGs)")
p_enh  <- plot_go_clean(ego_enh,  "H3K27ac Enhancer Peaks")
p_pro  <- plot_go_clean(ego_pro,  "H3K27ac Promoter Peaks")
p_diffpeak <- plot_go_clean(ego_diff, "H3K27ac Differential Peaks")

p_go <- (p_deg | p_enh | p_pro) +
  plot_annotation(
    title = "Functional enrichment of genes altered by perinatal morphine exposure",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0)
    )
  )

ggsave(
  "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/Fig2F_GOterms_MOR-SAL_vs_SAL-SAL_clean.pdf",
  p_go, width = 18, height = 4.5, device = cairo_pdf
)
ggsave(
  "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/Fig2F_GOterms_MOR-SAL_vs_SAL-SAL_clean.png",
  p_go, width = 18, height = 4.5, dpi = 600
)

p_go_all <- (p_deg | p_diffpeak) +
  plot_annotation(
    title = "Functional enrichment of genes altered by perinatal morphine exposure",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0)
    )
  )

ggsave(
  "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/Fig2F_GOterms_MOR-SAL_vs_SAL-SAL_diffpeaks_clean.pdf",
  p_go_all, width = 12, height = 4.5, device = cairo_pdf)
ggsave(
  "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/Fig2F_GOterms_MOR-SAL_vs_SAL-SAL_diffpeaks_clean.png",
  p_go_all, width = 12, height = 4.5, dpi = 600)
# ===============================================================
# Suppl S2B/S2C and immune panels
# (same logic as before, but using ego_deg_up/ego_deg_down with
# the corrected MOR–SAL / SAL–SAL orientation)
# ===============================================================

combine_go_for_bar <- function(ego, label, top_n = 10) {
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  ego@result %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(group = label)
}

df_combined <- bind_rows(
  combine_go_for_bar(ego_deg,  "DEG"),
  combine_go_for_bar(ego_enh,  "Enhancer"),
  combine_go_for_bar(ego_pro,  "Promoter")
)

pcom <- ggplot(df_combined,
               aes(x = reorder(Description, -log10(p.adjust)),
                   y = -log10(p.adjust), fill = group)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c(
    "DEG"      = col_deg,
    "Enhancer" = col_enhancer,
    "Promoter" = col_promoter
  )) +
  labs(
    title = "Top GO terms (MOR–SAL vs SAL–SAL)",
    x = NULL, y = expression(-log[10]~adjusted~italic(P)),
    fill = NULL
  ) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    plot.title  = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )

ggsave(
  "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/FigS2B_GOterms_MOR-SAL_vs_SAL-SAL_bar_com.pdf",
  pcom, width = 10, height = 5, device = cairo_pdf
)

## ---- S2C up/down bars (semantics fixed: up = MOR–SAL up) ----

combine_go_for_bar2 <- function(ego, label, top_n = 10) {
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  ego@result %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(group = label)
}

df_deg_combined <- bind_rows(
  combine_go_for_bar2(ego_deg_up,   "Upregulated in MOR–SAL"),
  combine_go_for_bar2(ego_deg_down, "Downregulated in MOR–SAL")
) %>%
  mutate(
    Description = str_wrap(Description, 50),
    GO = factor(Description, levels = rev(unique(Description)))
  )

p_deg_bar <- ggplot(df_deg_combined,
                    aes(x = reorder(Description, -log10(p.adjust)),
                        y = -log10(p.adjust),
                        fill = group)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Upregulated in MOR–SAL"   = "#56B4E9",
    "Downregulated in MOR–SAL" = "#E64B35"
  )) +
  labs(
    title = "Top GO terms for DEGs (MOR–SAL vs SAL–SAL)",
    x = NULL,
    y = expression(-log[10]~adjusted~italic(P)),
    fill = NULL
  ) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title  = element_text(face = "bold"),
    plot.title  = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )

ggsave(
  "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/FigS2C_GOterms_DEG_updown_bar.pdf",
  p_deg_bar, width = 8, height = 6, device = cairo_pdf
)

ggsave(
  "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/FigS2C_GOterms_DEG_updown_bar.png",
  p_deg_bar, width = 8, height = 6, dpi = 600
)

## ---- Immune GO: use DOWN genes (MOR–SAL < SAL–SAL) ----

immune_down <- ego_deg_down@result %>%
  filter(
    str_detect(Description,
               regex("immune|inflamm|cytokine|defense|leukocyte|macrophage|microglia|chemokine|interferon",
                     ignore_case = TRUE))
  ) %>%
  arrange(p.adjust)

write_tsv(immune_down,
          "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_downregulated_immune_GO.tsv")

# (the rest of your Sankey + dot immune plotting code can stay
# as-is, since it just reads that TSV; no path or orientation change)

overlap_peaks <- chip_tbl %>%
  dplyr::filter(SYMBOL %in% overlap_all3_genes) %>%
  dplyr::select(SYMBOL, region, seqnames, start, end, Fold, FDR)

cat("\nPeaks for 3-way overlap genes:\n")
print(overlap_peaks)

# 3) make an IGV-ready BED file for those peaks
#    (0-based start, 1-based end; standard BED convention)
bed_igv <- overlap_peaks %>%
  dplyr::mutate(
    chrom = as.character(seqnames),
    chromStart = as.integer(start) - 1L,
    chromEnd   = as.integer(end),
    name = paste(SYMBOL, region, sep = "_")
  ) %>%
  dplyr::select(chrom, chromStart, chromEnd, name)

igv_bed_file <- file.path(out_dir, "Fig3_overlap2_genes_forIGV.bed")
readr::write_tsv(bed_igv, igv_bed_file, col_names = FALSE)

cat("\n IGV BED written to:\n  ", igv_bed_file, "\n")

# (optional) if you also want DEG table rows for those genes:
overlap_deg_rows <- deg_tbl %>%
  dplyr::filter(gene_name %in% overlap_all3_genes) %>%
  dplyr::select(gene_name, log2FC, FDR, baseMean)

cat("\nDEG stats for 3-way overlap genes:\n")
print(overlap_deg_rows)

## ===============================================================
## Find strongest promoter diff peaks overlapping DEGs
##  - region == "Promoter"
##  - gene in DEGs ∩ Promoter
##  - rank by FDR first, then |Fold|
## ===============================================================

# genes in the DEG–Promoter overlap (the 53 in the Venn)
overlap_deg_enh_genes <- intersect(deg_genes, enhancer_genes)

length(overlap_deg_enh_genes)
print(overlap_deg_enh_genes)  # preview

enh_overlap_peaks <- chip_tbl %>%
  dplyr::filter(
    region == "Enhancer",
    SYMBOL %in% overlap_deg_enh_genes,
    !is.na(FDR),
    FDR < fdr_cut
  )

cat("Number of enhancer peaks linked to DEG genes (FDR <", fdr_cut, "): ",
    nrow(enh_overlap_peaks), "\n")

## ---- top overall (strongest + most significant) ----
top_enhancer_peaks <- enh_overlap_peaks %>%
  dplyr::arrange(FDR, dplyr::desc(abs(Fold))) %>%
  dplyr::select(SYMBOL, seqnames, start, end, Fold, FDR) %>%
  dplyr::slice_head(n = 5)

cat("\nTop enhancer diff peaks overlapping DEGs (ranked by FDR then |Fold|):\n")
print(top_enhancer_peaks)

## OPTIONAL: IGV BED for those top peaks
bed_top_enh <- top_enhancer_peaks %>%
  dplyr::mutate(
    chrom      = as.character(seqnames),
    chromStart = as.integer(start) - 1L,
    chromEnd   = as.integer(end),
    name       = paste(SYMBOL, sprintf("Fold=%.2f", Fold), sep = "_")
  ) %>%
  dplyr::select(chrom, chromStart, chromEnd, name)

enh_igv_bed <- file.path(out_dir, "Fig2_enhancer_DEG_overlap_topPeaks_forIGV.bed")
readr::write_tsv(bed_top_enh, enh_igv_bed, col_names = FALSE)
