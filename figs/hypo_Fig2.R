#!/usr/bin/env Rscript

# Hypothalamus MOR–SAL vs SAL–SAL
# Figure 2C (volcano), Figure 2E (Venn),
# Figure 2F + Supplementary GO panels, and immune GO.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
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

# ===============================================================
# Figure 2C — Volcano plots for differential peaks & DEGs
# MOR–SAL vs SAL–SAL
# ===============================================================

theme_pub <- theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# ---------- Input ----------
deg_file  <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/DEGs_treatment_only/SAL_SAL_vs_MOR_SAL_ALL.tsv"
chip_file <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/Diff_analysis/Annotated_MORSAL_vs_SALSAL.tsv"
sig_anno_chip <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/Diff_analysis/Annotated_MORSAL_vs_SALSAL.tsv"
sig_chip <- read_tsv(sig_anno_chip, show_col_types = FALSE)

# ---------- Parameters ----------
fdr_cut <- 0.05
lfc_cut <- 0.25
col_promoter <- "#56B4E9"
col_enhancer <- "#E69F00"
col_ns <- "grey80"

# ---------- DEG volcano ----------
deg <- read_tsv(deg_file, show_col_types = FALSE) %>%
  rename(FDR = padj) %>%  # rename padj -> FDR
  mutate(
    log2FoldChange = ifelse(is.na(log2FoldChange), 0, log2FoldChange),
    sig = ifelse(FDR < fdr_cut & abs(log2FoldChange) > lfc_cut, "Sig", "NS")
  )

top_deg <- deg %>%
  filter(sig == "Sig") %>%
  slice_max(abs(log2FoldChange), n = 10)

p_deg <- ggplot(deg, aes(x = log2FoldChange, y = -log10(FDR))) +
  geom_point(aes(color = sig), size = 1.8, alpha = 0.8) +
  scale_color_manual(values = c("Sig" = "#009E73", "NS" = col_ns)) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "grey50") +
  geom_text_repel(data = top_deg, aes(label = gene_name),
                  size = 3, color = "black", max.overlaps = 12) +
  labs(title = "RNA-seq (DEGs)",
       x = expression(log[2]~fold~change),
       y = expression(-log[10]~adjusted~italic(P))) +
  theme_pub +
  theme(legend.position = "none")

# ---------- DiffPeak volcano ----------
chip <- read_tsv(chip_file, show_col_types = FALSE) %>%
  mutate(
    Fold = -Fold,  # reverse direction for consistency (MOR–SAL > SAL–SAL)
    region = if_else(grepl("promoter", annotation, ignore.case = TRUE),
                     "Promoter", "Enhancer"),
    sig = ifelse(FDR < fdr_cut & abs(Fold) > 0, region, "NS")
  )

top_chip <- chip %>%
  filter(sig != "NS") %>%
  slice_max(abs(Fold), n = 10)

p_chip <- ggplot(chip, aes(x = Fold, y = -log10(FDR))) +
  geom_point(aes(color = sig), size = 1.8, alpha = 0.8) +
  scale_color_manual(values = c("Enhancer" = col_enhancer,
                                "Promoter" = col_promoter,
                                "NS" = col_ns)) +
  geom_hline(yintercept = -log10(fdr_cut), linetype = "dashed", color = "grey50") +
  geom_text_repel(data = top_chip, aes(label = SYMBOL),
                  size = 3, color = "black", max.overlaps = 12) +
  labs(title = "H3K27ac differential peaks",
       x = expression(log[2]~fold~change),
       y = expression(-log[10]~adjusted~italic(P))) +
  theme_pub +
  theme(legend.position = "top",
        legend.title = element_blank())

# ---------- Combine ----------
p_combined <- p_chip | p_deg +
  plot_annotation(
    title = "MOR–SAL vs SAL–SAL",
    subtitle = "Differential peaks (Enhancer/Promoter) and DEGs",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, color = "grey25", hjust = 0.5)
    )
  )

out_dir <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo"

ggsave(file.path(out_dir, "Fig2C_volcano_MOR-SAL_vs_SAL-SAL.pdf"),
       p_combined, width = 7, height = 3.5, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig2C_volcano_MOR-SAL_vs_SAL-SAL.png"),
       p_combined, width = 7, height = 3.5, dpi = 600)

p_combined

# ===============================================================
# Figure 2E — Numeric-only Venn diagram
# ===============================================================

library(VennDiagram)
library(grid)

# ---------- Input ----------
deg_file  <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/DEGs_treatment_only/SAL_SAL_vs_MOR_SAL_ALL.tsv"
chip_file <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/Diff_analysis/Annotated_MORSAL_vs_SALSAL.tsv"

# ---------- Parameters ----------
fdr_cut <- 0.05
lfc_cut <- 0.25

# ---------- Load DEGs ----------
deg <- read_tsv(deg_file, show_col_types = FALSE) %>%
  rename(FDR = padj, log2FC = log2FoldChange) %>%
  filter(!is.na(FDR) & FDR < fdr_cut & abs(log2FC) > lfc_cut)

deg_genes <- unique(deg$gene_name)

# ---------- Load differential peaks ----------
chip <- read_tsv(chip_file, show_col_types = FALSE) %>%
  mutate(
    region = if_else(grepl("promoter", annotation, ignore.case = TRUE),
                     "Promoter", "Enhancer")
  ) %>%
  filter(FDR < fdr_cut) %>%
  select(region, SYMBOL) %>%
  filter(!is.na(SYMBOL) & SYMBOL != "")

enhancer_genes <- unique(chip$SYMBOL[chip$region == "Enhancer"])
promoter_genes <- unique(chip$SYMBOL[chip$region == "Promoter"])

# ---------- Prepare gene sets ----------
venn_list <- list(
  DEGs = deg_genes,
  Enhancer = enhancer_genes,
  Promoter = promoter_genes
)

# ===============================================================
# Figure 2E — Numeric-only Venn (centered labels)
# ===============================================================

venn_plot <- draw.triple.venn(
  area1 = length(venn_list$DEGs),
  area2 = length(venn_list$Enhancer),
  area3 = length(venn_list$Promoter),
  n12 = length(intersect(venn_list$DEGs, venn_list$Enhancer)),
  n23 = length(intersect(venn_list$Enhancer, venn_list$Promoter)),
  n13 = length(intersect(venn_list$DEGs, venn_list$Promoter)),
  n123 = length(Reduce(intersect, venn_list)),

  category = c("DEGs", "Enhancer", "Promoter"),
  fill = c("#00A087", "#E69F00", "#56B4E9"),
  alpha = rep(0.7, 3),

  # ---- Label & number styling ----
  cat.cex = 1.3,
  cex = 1.3,
  cat.fontface = "bold",
  cat.col = "black",
  fontfamily = "sans",
  cat.fontfamily = "sans",

  # ---- Position adjustments ----
  cat.pos = c(-10, 10, 0),   # shifts labels inward
  cat.dist = c(0.05, 0.05, 0.05), # distance from circles
  margin = 0.08,              # adds padding so labels fit inside
  scaled = FALSE,
  euler.d = FALSE
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

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)

deg_file  <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/DEGs_treatment_only/SAL_SAL_vs_MOR_SAL_ALL.tsv"
chip_file <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/Diff_analysis/Annotated_MORSAL_vs_SALSAL.tsv"

fdr_cut <- 0.05
lfc_cut <- 0.25

deg <- read_tsv(deg_file, show_col_types = FALSE)

# --- DEGs ---
deg <- read_tsv(deg_file, show_col_types = FALSE) %>%
  filter(padj < fdr_cut & abs(log2FoldChange) > lfc_cut) %>%
  pull(gene_name) %>%
  unique()

deg_up <- read_tsv(deg_file, show_col_types = FALSE) %>%
  filter(padj < fdr_cut & log2FoldChange > lfc_cut) %>%
  pull(gene_name) %>%
  unique()

deg_down <- read_tsv(deg_file, show_col_types = FALSE) %>%
  filter(padj < fdr_cut & log2FoldChange < -lfc_cut) %>%
  pull(gene_name) %>%
  unique()

# --- Differential peaks (Enhancer & Promoter genes) ---
chip <- read_tsv(chip_file, show_col_types = FALSE)
enh_genes <- chip %>%
  filter(FDR < fdr_cut & !grepl("promoter", annotation, ignore.case = TRUE)) %>%
  pull(SYMBOL) %>%
  na.omit() %>%
  unique()

prom_genes <- chip %>%
  filter(FDR < fdr_cut, grepl("promoter", annotation, ignore.case = TRUE)) %>%
  pull(SYMBOL) %>% na.omit() %>% unique()

diff_genes <- unique(c(enh_genes, prom_genes))

# DEG enrichment
ego_deg <- enrichGO(
  gene          = deg,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
write_tsv(as.data.frame(ego_deg), "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_DEG_all_GO_BP_enrichment.tsv")

ego_deg_up <- enrichGO(
  gene          = deg_up,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
write_tsv(as.data.frame(ego_deg_up), "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_DEG_up_GO_BP_enrichment.tsv")

ego_deg_down <- enrichGO(
  gene          = deg_down,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)
write_tsv(as.data.frame(ego_deg_down), "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_DEG_down_GO_BP_enrichment.tsv")

# Enhancer GO
ego_enh <- enrichGO(
  gene = enh_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05
)
write_tsv(as.data.frame(ego_enh), "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_enhancer_peakall_GO_BP_enrichment.tsv")

# Promoter GO
ego_pro <- enrichGO(
  gene = prom_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05
)
write_tsv(as.data.frame(ego_pro), "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_promoter_peakall_GO_BP_enrichment.tsv")

ego_diff <- enrichGO(
  gene = diff_genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05
)
write_tsv(as.data.frame(ego_diff), "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_diff_peakall_GO_BP_enrichment.tsv")

plot_go_clean <- function(ego, title, top_n = 10) {
  df <- ego@result %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(
      Description = str_wrap(Description, 45),
      GO = factor(Description, levels = rev(unique(Description)))
    )

  ggplot(df, aes(x = GeneRatio, y = GO, size = Count, color = -log10(p.adjust))) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#56B4E9", high = "#E64B35") +
    labs(
      title = title,
      x = "Gene Ratio",
      y = NULL,
      color = expression(-log[10]~FDR)
    ) +
    theme_classic(base_size = 12, base_family = "Helvetica") +
    theme(
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      legend.position = "right",
      panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
    )
}

library(patchwork)

p_deg  <- plot_go_clean(ego_deg,  "RNA-seq (DEGs)")
p_enh  <- plot_go_clean(ego_enh,  "H3K27ac Enhancer Peaks")
p_pro  <- plot_go_clean(ego_pro,  "H3K27ac Promoter Peaks")

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

# --- combined GO bar (S2B) ---
combine_go_for_bar <- function(ego, label, top_n = 10) {
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

pcom = ggplot(df_combined, aes(x = reorder(Description, -log10(p.adjust)),
                        y = -log10(p.adjust), fill = group)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c("DEG" = "#009E73", "Enhancer" = "#E69F00", "Promoter" = "#56B4E9")) +
  labs(
    title = "Top GO terms (MOR–SAL vs SAL–SAL)",
    x = NULL, y = expression(-log[10]~adjusted~italic(P)),
    fill = NULL
  ) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )

ggsave(
  "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/FigS2B_GOterms_MOR-SAL_vs_SAL-SAL_bar_com.pdf",
  pcom, width = 10, height = 5, device = cairo_pdf
)

# ===============================================================
# Supplementary: GO enrichment up/down DEGs (S2C)
# ===============================================================

library(tidyverse)
library(ggplot2)

combine_go_for_bar <- function(ego, label, top_n = 10) {
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  ego@result %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(group = label)
}

df_deg_combined <- bind_rows(
  combine_go_for_bar(ego_deg_up,  "Downregulated in MOR–SAL"),
  combine_go_for_bar(ego_deg_down, "Upregulated in MOR–SAL")
) %>%
  mutate(
    Description = str_wrap(Description, 50),
    GO = factor(Description, levels = rev(unique(Description)))
  )

unique(df_deg_combined$group)

p_deg_bar <- ggplot(df_deg_combined,
                    aes(x = reorder(Description, -log10(p.adjust)),
                        y = -log10(p.adjust),
                        fill = group)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Upregulated in MOR–SAL" = "#56B4E9",
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
    axis.title = element_text(face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
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

# ===============================================================
# Compact bidirectional GO bar (S2C alt)
# ===============================================================

library(ggrepel)
library(patchwork)

combine_go_for_bar <- function(ego, label, direction = "up", top_n = 8) {
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  ego@result %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(group = label,
           direction = direction)
}

df_up <- combine_go_for_bar(ego_deg_up, "Upregulated", direction = "Up")
df_down <- combine_go_for_bar(ego_deg_down, "Downregulated", direction = "Down")

df_bidirectional <- bind_rows(df_up, df_down) %>%
  mutate(
    log10FDR = -log10(p.adjust),
    log10FDR = if_else(direction == "Down", -log10FDR, log10FDR),
    Description = str_wrap(Description, 55),
    GO = factor(Description, levels = unique(Description))
  )

p_bidirectional <- ggplot(df_bidirectional,
       aes(x = log10FDR, y = reorder(Description, log10FDR), fill = direction)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = c("Up" = "#56B4E9", "Down" = "#E64B35")) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  labs(
    title = "Top GO terms for DEGs (MOR–SAL vs SAL–SAL)",
    x = expression(-log[10]~adjusted~italic(P)),
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "top"
  )

ggsave(
  "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/FigS2C_GOterms_DEG_updown_bidirectional.pdf",
  p_bidirectional, width = 7, height = 5, device = cairo_pdf
)

ggsave(
  "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/FigS2C_GOterms_DEG_updown_bidirectional.png",
  p_bidirectional, width = 7, height = 5, dpi = 600
)

# ===============================================================
# EVEN MORE COMPACT bidirectional version (final S2C)
# ===============================================================

shorten_terms <- function(x) {
  x <- str_replace_all(x, "regulation of ", "reg. of ")
  x <- str_replace_all(x, "positive ", "pos. ")
  x <- str_replace_all(x, "negative ", "neg. ")
  x <- str_replace_all(x, "process", "proc.")
  x <- str_replace_all(x, "response to ", "resp. to ")
  x <- str_replace_all(x, "cellular ", "cell ")
  x
}

combine_go_for_bar <- function(ego, label, direction = "up", top_n = 8) {
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  ego@result %>%
    arrange(p.adjust) %>%
    slice_head(n = top_n) %>%
    mutate(group = label,
           direction = direction)
}

df_up <- combine_go_for_bar(ego_deg_up, "Upregulated", direction = "Up")
df_down <- combine_go_for_bar(ego_deg_down, "Downregulated", direction = "Down")

df_bidirectional <- bind_rows(df_up, df_down) %>%
  mutate(
    Description = shorten_terms(Description),
    Description = str_wrap(Description, 35),
    log10FDR = -log10(p.adjust),
    log10FDR = if_else(direction == "Down", -log10FDR, log10FDR),
    GO = factor(Description, levels = unique(Description))
  )

theme_pub_compact <- theme_classic(base_size = 11, base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(size = 9.5, color = "black", lineheight = 0.9),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.position = "top",
    plot.margin = margin(5, 10, 5, 5)
  )

p_bidirectional <- ggplot(df_bidirectional,
       aes(x = log10FDR, y = reorder(Description, log10FDR), fill = direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Up" = "#56B4E9", "Down" = "#E64B35")) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
  labs(
    title = "Top GO terms for DEGs (MOR–SAL vs SAL–SAL)",
    x = expression(-log[10]~adjusted~italic(P)),
    y = NULL,
    fill = NULL
  ) +
  theme_pub_compact

ggsave(
  "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/FigS2C_GOterms_DEG_updown_compact.pdf",
  p_bidirectional, width = 7, height = 7, device = cairo_pdf
)

ggsave(
  "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/FigS2C_GOterms_DEG_updown_compact.png",
  p_bidirectional, width = 7, height = 7, dpi = 600
)

# ===============================================================
# Immune GO extraction + Fig3E plots
# ===============================================================

library(dplyr)
library(stringr)
library(readr)

immune_down <- ego_deg_up@result %>%
  filter(
    str_detect(Description, regex("immune|inflamm|cytokine|defense|leukocyte|macrophage|microglia|chemokine|interferon", ignore_case = TRUE))
  ) %>%
  arrange(p.adjust)

write_tsv(immune_down,
          "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_downregulated_immune_GO.tsv")

# === Figure 3E — Downregulated immune GO terms (Hypothalamus, MOR–SAL vs SAL–SAL) ===

library(tidyverse)
library(ggplot2)
library(stringr)

go_file <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_downregulated_immune_GO.tsv"

immune_go <- read_tsv(go_file, show_col_types = FALSE) %>%
  arrange(p.adjust) %>%
  slice_head(n = 12) %>%                                  # top 12 terms
  mutate(
    Description = str_wrap(Description, 55),
    GO = factor(Description, levels = rev(unique(Description)))
  )

p_immune <- ggplot(immune_go, aes(x = -log10(p.adjust), y = GO)) +
  geom_col(fill = "#E64B35", width = 0.7) +
  geom_text(
    aes(label = Count), hjust = -0.2, size = 3, color = "black"
  ) +
  labs(
    title = "Downregulated immune pathways in adult hypothalamus (MOR–SAL < SAL–SAL)",
    x = expression(-log[10]~adjusted~italic(P)),
    y = NULL
  ) +
  theme_classic(base_size = 12, base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    plot.margin = margin(10, 15, 10, 10)
  ) +
  coord_cartesian(clip = "off")

out_dir <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo"
ggsave(file.path(out_dir, "Fig3E_downregulated_immune_GOterms.pdf"),
       p_immune, width = 6, height = 4, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig3E_downregulated_immune_GOterms.png"),
       p_immune, width = 6, height = 4, dpi = 600)

p_immune

# === Figure 3E — Sankey + dot (immune GO) ===
library(tidyverse)
library(ggalluvial)
library(ggplot2)
library(stringr)

go_file <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/MOR_SAL_vs_SAL_SAL_downregulated_immune_GO.tsv"
deg_file <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/DEGs_treatment_only/SAL_SAL_vs_MOR_SAL_sig.tsv"  # if you want to link genes

go_df <- read_tsv(go_file, show_col_types = FALSE) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10) %>%
  mutate(
    Description = str_wrap(Description, 45),
    GO = factor(Description, levels = rev(unique(Description)))
  )

go_gene <- go_df %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::select(GO = Description, Gene = geneID)

gene_top <- go_gene %>%
  count(Gene, sort = TRUE) %>%
  slice_head(n = 18) %>%
  pull(Gene)

go_gene <- go_gene %>% filter(Gene %in% gene_top)

p_sankey <- ggplot(go_gene, aes(axis1 = Gene, axis2 = GO)) +
  geom_alluvium(fill = "#E64B35", alpha = 0.6, width = 1/12) +
  geom_stratum(fill = "grey95", color = "grey60") +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum)),
            size = 3, family = "Helvetica", color = "black") +
  theme_void(base_size = 12) +
  labs(title = "Downregulated Immune GO Terms (MOR–SAL < SAL–SAL)") +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5, margin = margin(5,0,5,0))
  )

ggsave("/Users/yanmiaodu/Downloads/MOR_SAL/hypo/temp_sankey.pdf",
       p_sankey, width = 6, height = 5, device = cairo_pdf)

p_dot <- go_df %>%
  mutate(GeneCount = str_count(geneID, "/") + 1) %>%
  ggplot(aes(x = GeneCount, y = GO, size = GeneCount, color = -log10(p.adjust))) +
  geom_point() +
  scale_color_gradient(low = "#56B4E9", high = "#E64B35") +
  scale_size_continuous(range = c(2, 6)) +
  labs(
    title = "Top Downregulated immune GO term in MOR.SAL",
    x = "Gene Count",
    y = NULL,
    color = expression(-log[10]~italic(P))
  ) +
  theme_classic(base_size = 8, base_family = "Helvetica") +
  theme(
    axis.text.y = element_text(size = 7, color = "black"),
    axis.text.x = element_text(size =7, color = "black"),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
    legend.position = "right",
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)
  )

ggsave("/Users/yanmiaodu/Downloads/MOR_SAL/hypo/temp_dot.pdf",
       p_dot, width = 4.5, height = 3.2, device = cairo_pdf)
ggsave("/Users/yanmiaodu/Downloads/MOR_SAL/hypo/temp_dot.png",
       p_dot, width = 4.5, height = 3.2, dpi = 600)

p_combined <- p_sankey + p_dot +
  plot_layout(widths = c(0.45, 0.55)) +
  plot_annotation(
    title = "Downregulated immune response pathways in adult hypothalamus",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  )

ggsave(file.path(out_dir, "Fig3E_downregulated_immune_GO_sankeydot_balanced.pdf"),
       p_combined, width = 9, height = 4.5, device = cairo_pdf)
ggsave(file.path(out_dir, "Fig3E_downregulated_immune_GO_sankeydot_balanced.png"),
       p_combined, width = 9, height = 4.5, dpi = 600)

p_combined
