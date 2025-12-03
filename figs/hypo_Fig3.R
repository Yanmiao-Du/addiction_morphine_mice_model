## ===============================================================
## Hypothalamus – Figure 3 (MOR–LPS vs SAL–LPS)
## 3B summary bars, 3C/3D volcano, 3E Venn (numeric + no-text),
## 3F GO terms
## ===============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(ggtext)
  library(patchwork)
  library(VennDiagram)
  library(grid)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

## ---------- Paths ----------
ROOT      <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo"
chip_anno_file <- file.path(ROOT, "Diff_analysis/Annotated_MORLPS_vs_SALLPS.tsv")
deg_file       <- file.path(ROOT, "DEGs_treatment_only/SAL_LPS_vs_MOR_LPS_ALL.tsv")
diffbind_rda   <- file.path(ROOT, "hypo_diffbind_DESeq2_allContrasts.RData")
out_dir        <- ROOT

## ---------- Parameters ----------
fdr_cut <- 0.05
lfc_cut <- 0.25

col_deg       <- "#009E73"  # DEG bar / points
col_enhancer  <- "#E69F00"
col_promoter  <- "#56B4E9"
col_up        <- "#00A087"  # increased in MOR-LPS (for volcano)
col_down      <- "#E64B35"  # decreased in MOR-LPS
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

## ===============================================================
## 1. Load data
## ===============================================================

## --- DEGs (SAL_LPS_vs_MOR_LPS_ALL.tsv) ---
deg_raw <- read_tsv(deg_file, show_col_types = FALSE)

deg_tbl <- deg_raw %>%
  dplyr::rename(FDR = padj) %>%
  mutate(
    ## file is SAL_LPS_vs_MOR_LPS -> flip so + = MOR-LPS / SAL-LPS
    log2FC = -log2FoldChange,
    log2FC = ifelse(is.na(log2FC), 0, log2FC),
    sig = if_else(FDR < fdr_cut & abs(log2FC) > lfc_cut, "Sig", "NS")
  )

## --- DiffBind full results (all peaks, not only significant) ---
load(diffbind_rda)   # loads dbObj_clean and others

## contrast 2 was MOR-LPS vs SAL-LPS in your earlier code
mor_lps_vs_sal_lps_all <- DiffBind::dba.report(
  dbObj_clean,
  contrast = 2,
  th = 1  # th=1 -> all peaks with FDR
)

chip_all <- as.data.frame(mor_lps_vs_sal_lps_all) %>%
  as_tibble() %>%
  mutate(
    seqnames = as.character(seqnames),
    start = start,
    end = end
  )

## --- Annotation file (only FDR<0.05) for region + SYMBOL ---
chip_anno <- read_tsv(chip_anno_file, show_col_types = FALSE)

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

## ===============================================================
## 2. Figure 3B — unified summary barplot (DEG / Enhancer / Promoter)
##      - counts by FDR only (no LFC filter)
## ===============================================================

deg_summary <- deg_tbl %>%
  filter(FDR < fdr_cut) %>%
  mutate(region = "DEG") %>%
  count(region, name = "n")

chip_summary <- chip_sig %>%
  count(region, name = "n")

summary_df <- bind_rows(deg_summary, chip_summary) %>%
  mutate(region = factor(region, levels = c("DEG", "Enhancer", "Promoter")))

cat_colors <- c("DEG" = col_deg,
                "Enhancer" = col_enhancer,
                "Promoter" = col_promoter)

p_fig3B <- ggplot(summary_df, aes(x = region, y = n, fill = region)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_text(aes(label = n),
            vjust = -0.3, size = 3.8, fontface = "bold", color = "black") +
  scale_fill_manual(values = cat_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = NULL,
    y = "Number of significant features",
    title = "MOR–LPS vs SAL–LPS"
  ) +
  theme_pub +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "Fig3B_MOR-LPS_vs_SAL-LPS_unified_pubready.pdf"),
       p_fig3B, width = 4, height = 3.2, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig3B_MOR-LPS_vs_SAL-LPS_unified_pubready.png"),
       p_fig3B, width = 4, height = 3.2, dpi = 600)

## ===============================================================
## 3C & 3D – Volcano plots (RNA-seq DEGs + H3K27ac peaks)
## ===============================================================

## ---- choose top DEGs to label (strong FC + low FDR) ----
top_deg <- deg_tbl %>%
  filter(sig == "Sig") %>%
  arrange(FDR, desc(abs(log2FC))) %>%   # strong FC and low FDR
  slice_head(n = 10)                    # change n if you want more/fewer labels

## --- DEG volcano (3C) ---
p_deg_volcano <- ggplot(deg_tbl, aes(x = log2FC, y = -log10(FDR))) +
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
    x = expression(log[2]~fold~change~"(MOR-LPS / SAL-LPS)"),
    y = expression(-log[10]~adjusted~italic(P))
  ) +
  theme_pub +
  theme(legend.position = "none")

## ---- choose top ChIP peaks to label ----
top_chip <- chip_tbl %>%
  filter(sig != "NS") %>%
  arrange(FDR, desc(abs(Fold))) %>%
  slice_head(n = 10)

## ---- ChIP volcano data: drop only the FDR=1 floor ----
chip_plot_tbl <- chip_tbl %>%
  filter(
    !is.na(FDR),
    FDR < 0.999 
  )

## --- H3K27ac volcano (3D) ---
p_chip_volcano <- ggplot(chip_plot_tbl, aes(x = Fold, y = -log10(FDR))) +
  geom_point(aes(color = sig), size = 1.8, alpha = 0.8) +
  scale_color_manual(values = c(
    "Enhancer" = "#E69F00",
    "Promoter" = "#56B4E9",
    "NS"       = "grey80"
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
    x = expression(log[2]~fold~change~"(MOR-LPS / SAL-LPS)"),
    y = expression(-log[10]~adjusted~italic(P))
  ) +
  theme_pub +
  theme(
    legend.position = "top",
    legend.title    = element_blank()
  )

## --- combine C+D panel ---
p_fig3CD <- p_chip_volcano | p_deg_volcano +
  plot_annotation(
    title = "MOR–LPS vs SAL–LPS",
    subtitle = "Differential peaks (Enhancer/Promoter) and DEGs",
    theme = theme(
      plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, color = "grey25", hjust = 0.5)
    )
  )

ggsave(file.path(out_dir, "Fig3C_D_volcano_MOR-LPS_vs_SAL-LPS.pdf"),
       p_fig3CD, width = 7, height = 3.5, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig3C_D_volcano_MOR-LPS_vs_SAL-LPS.png"),
       p_fig3CD, width = 7, height = 3.5, dpi = 600)


## ===============================================================
## 4. Figure 3E — Venn diagram (numeric + no-text)
##      sets defined by FDR<0.05 only (no LFC filter)
## ===============================================================

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

## --- numeric Venn (keep this) ---
venn_numeric <- draw.triple.venn(
  area1 = length(venn_list$DEGs),
  area2 = length(venn_list$Enhancer),
  area3 = length(venn_list$Promoter),
  n12   = length(intersect(venn_list$DEGs, venn_list$Enhancer)),
  n23   = length(intersect(venn_list$Enhancer, venn_list$Promoter)),
  n13   = length(intersect(venn_list$DEGs, venn_list$Promoter)),
  n123  = length(Reduce(intersect, venn_list)),
  category = c("DEGs", "Enhancer", "Promoter"),
  fill = c(col_deg, col_enhancer, col_promoter),
  alpha = rep(0.7, 3),
  cat.cex = 1.3,
  cex = 1.3,
  cat.fontface = "bold",
  cat.col = "black",
  fontfamily = "sans",
  cat.fontfamily = "sans",
  margin = 0.08,
  scaled = FALSE,
  euler.d = FALSE
)

pdf(file.path(out_dir, "Fig3E_Venn_MOR-LPS_vs_SAL-LPS_numeric.pdf"),
    width = 4.5, height = 4.5, family = "Helvetica")
grid.draw(venn_numeric)
dev.off()

png(file.path(out_dir, "Fig3E_Venn_MOR-LPS_vs_SAL-LPS_numeric.png"),
    width = 1200, height = 1200, res = 300, type = "cairo")
grid.draw(venn_numeric)
dev.off()

## --- no-number Venn (only circles + category labels) ---
venn_notext <- draw.triple.venn(
  area1 = length(venn_list$DEGs),
  area2 = length(venn_list$Enhancer),
  area3 = length(venn_list$Promoter),
  n12   = length(intersect(venn_list$DEGs, venn_list$Enhancer)),
  n23   = length(intersect(venn_list$Enhancer, venn_list$Promoter)),
  n13   = length(intersect(venn_list$DEGs, venn_list$Promoter)),
  n123  = length(Reduce(intersect, venn_list)),
  fill = c(col_deg, col_enhancer, col_promoter),
  alpha = rep(0.7, 3),
  cat.cex = 1.3,
  cex = 1.3,
  cat.fontface = "bold",
  cat.col = "black",
  fontfamily = "sans",
  cat.fontfamily = "sans",
  margin = 0.08,
  scaled = FALSE,
  euler.d = FALSE
)

pdf(file.path(out_dir, "Fig3E_Venn_MOR-LPS_vs_SAL-LPS_notext.pdf"),
    width = 4.5, height = 4.5, family = "Helvetica")
grid.draw(venn_notext)
dev.off()

png(file.path(out_dir, "Fig3E_Venn_MOR-LPS_vs_SAL-LPS_notext.png"),
    width = 1200, height = 1200, res = 300, type = "cairo")
grid.draw(venn_notext)
dev.off()

## ===============================================================
## 5. Figure 3F — GO BP enrichment
##      - FDR < 0.05 sets, NO LFC filter
## ===============================================================

deg_sig <- deg_tbl %>%
  filter(FDR < fdr_cut) %>%
  pull(gene_name) %>%
  unique()

enh_genes_sig <- enhancer_genes
prom_genes_sig <- promoter_genes

ego_deg <- enrichGO(
  gene          = deg_sig,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ego_enh <- enrichGO(
  gene          = enh_genes_sig,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ego_pro <- enrichGO(
  gene          = prom_genes_sig,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

diff_genes <- Reduce(intersect,
                     list(enh_genes_sig, prom_genes_sig))
ego_diffpeak <- enrichGO(
  gene          = diff_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

ego_all <- enrichGO(
  gene          = Reduce(union, list(deg_sig, enh_genes_sig, prom_genes_sig)),
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

write_tsv(as.data.frame(ego_all),
          file.path(out_dir, "Fig3F_GO_BP_AllFeatures_MOR-LPS_vs_SAL-LPS.tsv"))

write_tsv(as.data.frame(ego_diffpeak),
          file.path(out_dir, "Fig3F_GO_BP_DiffPeaks_MOR-LPS_vs_SAL-LPS.tsv"))

write_tsv(as.data.frame(ego_deg),
          file.path(out_dir, "Fig3F_GO_BP_DEG_MOR-LPS_vs_SAL-LPS.tsv"))
write_tsv(as.data.frame(ego_enh),
          file.path(out_dir, "Fig3F_GO_BP_Enhancer_MOR-LPS_vs_SAL-LPS.tsv"))
write_tsv(as.data.frame(ego_pro),
          file.path(out_dir, "Fig3F_GO_BP_Promoter_MOR-LPS_vs_SAL-LPS.tsv"))


## ===============================================================
## 5. Figure 3F — GO BP enrichment (fixed axes / no overlap / no log10 error)
## ===============================================================

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



p_deg_go <- plot_go_clean(ego_deg, "RNA-seq (DEGs)")
p_enh_go <- plot_go_clean(ego_enh, "H3K27ac Enhancer Peaks")
p_pro_go <- plot_go_clean(ego_pro, "H3K27ac Promoter Peaks")

p_fig3F <- (p_deg_go | p_enh_go | p_pro_go) +
  plot_annotation(
    title = "Functional enrichment (MOR–LPS vs SAL–LPS)",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0)
    )
  )

ggsave(file.path(out_dir, "Fig3F_GOterms_MOR-LPS_vs-SAL-LPS.pdf"),
       p_fig3F, width = 16, height = 4, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig3F_GOterms_MOR-LPS_vs-SAL-LPS.png"),
       p_fig3F, width = 16, height = 4, dpi = 600)



## ===============================================================
## Identify the 3-way overlap genes and export regions for IGV
## ===============================================================

# 1) which genes are in all three sets? (DEG ∩ Enhancer ∩ Promoter)
overlap_all3_genes <- Reduce(intersect,
                             list(deg_genes, enhancer_genes, promoter_genes))

cat("Genes in 3-way overlap (should be length 2):\n")
print(overlap_all3_genes)

# 2) grab the ChIP peaks for those genes (both enhancer + promoter)
#    Annotated_MORLPS_vs_SALLPS.tsv from ChIPseeker has seqnames/start/end/SYMBOL
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
overlap_deg_prom_genes <- intersect(deg_genes, promoter_genes)

length(overlap_deg_prom_genes)
print(overlap_deg_prom_genes[1:5])  # preview

promoter_overlap_peaks <- chip_tbl %>%
  dplyr::filter(
    region == "Promoter",
    SYMBOL %in% overlap_deg_prom_genes,
    !is.na(FDR),
    FDR < fdr_cut
  )

cat("Number of promoter peaks linked to DEG genes (FDR <", fdr_cut, "): ",
    nrow(promoter_overlap_peaks), "\n")

## ---- top overall (strongest + most significant) ----
top_promoter_peaks <- promoter_overlap_peaks %>%
  dplyr::arrange(FDR, dplyr::desc(abs(Fold))) %>%
  dplyr::select(SYMBOL, seqnames, start, end, Fold, FDR) %>%
  dplyr::slice_head(n = 5)

cat("\nTop promoter diff peaks overlapping DEGs (ranked by FDR then |Fold|):\n")
print(top_promoter_peaks)

## OPTIONAL: split into strongest up and strongest down
top_promoter_up <- promoter_overlap_peaks %>%
  dplyr::filter(Fold > 0) %>%
  dplyr::arrange(FDR, dplyr::desc(Fold)) %>%
  dplyr::select(SYMBOL, seqnames, start, end, Fold, FDR, baseMean) %>%
  dplyr::slice_head(n = 5)

top_promoter_down <- promoter_overlap_peaks %>%
  dplyr::filter(Fold < 0) %>%
  dplyr::arrange(FDR, Fold) %>%     # more negative first
  dplyr::select(SYMBOL, seqnames, start, end, Fold, FDR, baseMean) %>%
  dplyr::slice_head(n = 5)

cat("\nTop 5 UP promoter peaks overlapping DEGs:\n")
print(top_promoter_up)

cat("\nTop 5 DOWN promoter peaks overlapping DEGs:\n")
print(top_promoter_down)

## OPTIONAL: IGV BED for those top peaks
bed_top_promoter <- top_promoter_peaks %>%
  dplyr::mutate(
    chrom      = as.character(seqnames),
    chromStart = as.integer(start) - 1L,
    chromEnd   = as.integer(end),
    name       = paste(SYMBOL, sprintf("Fold=%.2f", Fold), sep = "_")
  ) %>%
  dplyr::select(chrom, chromStart, chromEnd, name)

promoter_igv_bed <- file.path(out_dir, "Fig3_promoter_DEG_overlap_topPeaks_forIGV.bed")
readr::write_tsv(bed_top_promoter, promoter_igv_bed, col_names = FALSE)

cat("\nIGV BED for top promoter–DEG overlap peaks written to:\n  ",
    promoter_igv_bed, "\n")
