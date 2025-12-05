#!/usr/bin/env Rscript
# ============================================================
# Hypothalamus H3K27ac – Early-life morphine & LPS analysis
#
# - Uses DiffBind DESeq2 results (all contrasts) from:
#       hypo_diffbind_DESeq2_allContrasts.RData
# - Produces:
#   * 4-panel scatter figure of baseline vs LPS responses
#   * Classification of morphine-distorted LPS peaks
#   * GO BP enrichment (clusterProfiler) for peak classes
#   * GREAT GO:BP enrichment for peak classes and directions
#   * Summary dot-heatmaps for GREAT results

# ============================================================

## ================== SETUP & DATA ============================
suppressPackageStartupMessages({
  library(DiffBind)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ChIPseeker)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(enrichplot)
  library(rGREAT)
  library(readr)
  library(stringr)
  library(forcats)
})

## ---- base directory & input RData ----
base_dir <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo"
setwd(base_dir)

load("hypo_diffbind_DESeq2_allContrasts.RData")

## ================== FUNCTIONS ========================

## ---- find DBA object in loaded .RData ----
get_dba_object <- function() {
  dba_name <- ls()[sapply(ls(), function(x) inherits(get(x), "DBA"))][1]
  if (is.na(dba_name)) stop("No DBA object found in workspace.")
  get(dba_name)
}

## ---- extract contrast as data.frame with peak_id ----
get_res <- function(db, idx, name) {
  df <- as.data.frame(dba.report(db, contrast = idx, th = 1))
  df$peak_id  <- rownames(df)
  df$contrast <- name
  df
}

## ---- write significant GO BP table (clusterProfiler) ----
save_go_tsv <- function(go_obj, file) {
  if (is.null(go_obj)) return(invisible(NULL))
  df <- as.data.frame(go_obj)
  if (!nrow(df)) return(invisible(NULL))

  df_sig <- df %>%
    dplyr::filter(p.adjust < 0.05) %>%
    arrange(p.adjust)

  if (!nrow(df_sig)) return(invisible(NULL))

  write.table(df_sig, file = file,
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("✅ Saved: ", file)
}

## ---- get significant peaks for a DiffBind contrast as GRanges ----
get_sig_peaks_gr <- function(db, contrast_idx,
                             fdr_cut = 0.05, lfc_cut = 0.25) {
  df <- as.data.frame(dba.report(db, contrast = contrast_idx, th = 1))
  df$peak_id <- rownames(df)

  df_sig <- df %>%
    filter(FDR < fdr_cut & abs(Fold) > lfc_cut)

  message("Contrast ", contrast_idx, ": ", nrow(df_sig), " significant peaks.")

  GRanges(
    seqnames = df_sig$seqnames,
    ranges   = IRanges(start = df_sig$start, end = df_sig$end),
    peak_id  = df_sig$peak_id
  )
}

## ---- run GREAT and return full & cleaned tables ----
run_great_clean <- function(gr, out_prefix,
                            size_max = 3000,
                            fe_min   = 1.5,
                            padj_cut = 0.05) {
  res <- great(
    gr         = gr,
    gene_sets  = "GO:BP",
    tss_source = "TxDb.Mmusculus.UCSC.mm10.knownGene"
  )

  tb <- getEnrichmentTable(res)

  write.table(tb,
              paste0("GREAT_full_", out_prefix, "_BP.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("✅ Saved GREAT_full_", out_prefix, "_BP.tsv")

  tb_clean <- tb %>%
    dplyr::filter(
      gene_set_size   < size_max,
      fold_enrichment > fe_min,
      p_adjust_hyper  < padj_cut
    ) %>%
    arrange(p_adjust_hyper)

  write.table(tb_clean,
              paste0("GREAT_clean_", out_prefix, "_BP.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message("✅ Saved GREAT_clean_", out_prefix, "_BP.tsv")

  invisible(list(full = tb, clean = tb_clean))
}

## ---- helper: make GRanges from df ----
make_gr <- function(df) {
  GRanges(
    seqnames = df$seqnames,
    ranges   = IRanges(start = df$start, end = df$end)
  )
}

## ================== DIFFBIND MERGE & SCATTER FIGURE =========

db <- get_dba_object()
dba.show(db, bContrasts = TRUE)
# 1 = MOR-LPS vs MOR-SAL
# 2 = MOR-LPS vs SAL-LPS
# 5 = MOR-SAL vs SAL-SAL
# 6 = SAL-LPS vs SAL-SAL

## ---- extract relevant contrasts ----
base     <- get_res(db, 5, "MOR-SAL_vs_SAL-SAL")
lps      <- get_res(db, 2, "MOR-LPS_vs_SAL-LPS")
sal_lps  <- get_res(db, 6, "SAL-LPS_vs_SAL-SAL")
mor_only <- get_res(db, 1, "MOR-LPS_vs_MOR-SAL")

## ---- merge log2FCs + FDRs ----
base2 <- base %>% dplyr::select(
  peak_id, seqnames, start, end,
  log2FC_base = Fold, FDR_base = FDR
)

lps2 <- lps %>% dplyr::select(
  peak_id,
  log2FC_lps = Fold, FDR_lps = FDR
)

sal2 <- sal_lps %>%
  dplyr::mutate(log2FC_SAL_LPS = -Fold) %>%  # flip sign so positive = normal LPS up
  dplyr::select(peak_id, log2FC_SAL_LPS, FDR_SAL_LPS = FDR)

mor2 <- mor_only %>% dplyr::select(
  peak_id,
  log2FC_MOR_LPS_direct = Fold,
  FDR_MOR_LPS_direct    = FDR
)

merged4 <- base2 %>%
  inner_join(lps2,  by = "peak_id") %>%
  inner_join(sal2,  by = "peak_id") %>%
  inner_join(mor2,  by = "peak_id")

## ---- correlation + slope for each panel ----
r_A     <- cor(merged4$log2FC_base,    merged4$log2FC_SAL_LPS,        use = "complete")
slope_A <- coef(lm(log2FC_SAL_LPS      ~ log2FC_base, merged4))[2]

r_B     <- cor(merged4$log2FC_base,    merged4$log2FC_MOR_LPS_direct, use = "complete")
slope_B <- coef(lm(log2FC_MOR_LPS_direct ~ log2FC_base, merged4))[2]

r_C     <- cor(merged4$log2FC_base,    merged4$log2FC_lps,            use = "complete")
slope_C <- coef(lm(log2FC_lps         ~ log2FC_base, merged4))[2]

r_D     <- cor(merged4$log2FC_SAL_LPS, merged4$log2FC_lps,            use = "complete")
slope_D <- coef(lm(log2FC_lps         ~ log2FC_SAL_LPS, merged4))[2]

## ---- axis limits & theme ----
xlim_all <- c(-0.50, 0.60)
ylim_all <- c(-1, 1)

lab_x <- xlim_all[1] + 0.02 * diff(xlim_all)
lab_y <- ylim_all[2] - 0.05 * diff(ylim_all)

base_theme <- theme_classic(base_size = 10) +
  theme(
    plot.title   = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.text    = element_text(size = 8)
  )

## ---- Panel A: normal LPS response vs baseline ----
p_sal <- ggplot(merged4, aes(log2FC_base, log2FC_SAL_LPS)) +
  geom_point(alpha = 0.25, size = 0.3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
  annotate("text", x = lab_x, y = lab_y,
           label = sprintf("r = %.2f\nslope = %.2f", r_A, slope_A),
           hjust = 0, vjust = 1, size = 3) +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +
  base_theme +
  labs(
    title = "Normal immune enhancer activation",
    x = "Baseline (MOR–SAL vs SAL–SAL)",
    y = "Normal LPS response (SAL–LPS vs SAL–SAL)"
  )

## ---- Panel B: direct MOR-LPS vs MOR-SAL ----
p_mor <- ggplot(merged4, aes(log2FC_base, log2FC_MOR_LPS_direct)) +
  geom_point(alpha = 0.25, size = 0.3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
  annotate("text", x = lab_x, y = lab_y,
           label = sprintf("r = %.2f\nslope = %.2f", r_B, slope_B),
           hjust = 0, vjust = 1, size = 3) +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +
  base_theme +
  labs(
    title = "Early-life morphine alters LPS response amplitude",
    x = "Baseline (MOR–SAL vs SAL–SAL)",
    y = "LPS response in early-life morphine group\n(MOR–LPS vs MOR–SAL)"
  )

## ---- Panel C: MOR-LPS vs SAL-LPS vs baseline ----
p_direct <- ggplot(merged4, aes(log2FC_base, log2FC_lps)) +
  geom_point(alpha = 0.25, size = 0.3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
  annotate("text", x = lab_x, y = lab_y,
           label = sprintf("r = %.2f\nslope = %.2f", r_C, slope_C),
           hjust = 0, vjust = 1, size = 3) +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +
  base_theme +
  labs(
    title = "Early-life morphine reprograms LPS responsiveness",
    x = "Baseline (MOR–SAL vs SAL–SAL)",
    y = "LPS response after early-life morphine\n(MOR–LPS vs SAL–LPS)"
  )

## ---- Panel D: MOR-LPS vs SAL-LPS LPS responses ----
p_compare <- ggplot(merged4, aes(log2FC_SAL_LPS, log2FC_lps)) +
  geom_point(alpha = 0.25, size = 0.3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
  annotate("text", x = lab_x, y = lab_y,
           label = sprintf("r = %.2f\nslope = %.2f", r_D, slope_D),
           hjust = 0, vjust = 1, size = 3) +
  coord_cartesian(xlim = xlim_all, ylim = ylim_all) +
  base_theme +
  labs(
    title = "Early-life morphine shifts LPS response pattern",
    x = "Normal LPS response (SAL–LPS vs SAL–SAL)",
    y = "LPS response after early-life morphine\n(MOR–LPS vs SAL–LPS)"
  )

p_all <- (p_sal | p_mor) /
         (p_direct | p_compare) +
  plot_annotation(tag_levels = "a")

ggsave("hypo_H3K27ac_scatter_4panel_NatureStyle.pdf",
       p_all, width = 10, height = 6, useDingbats = FALSE)
ggsave("hypo_H3K27ac_scatter_4panel_NatureStyle.png",
       p_all, width = 10, height = 6, dpi = 300)

## ================== DEFINE SIGNIFICANT LPS PEAKS ==============

lfc_cut   <- 0.25
fdr_cut   <- 0.05
delta_cut <- 0.30   # morphine-induced change in LPS response

merged4 <- merged4 %>%
  mutate(
    delta_LPS = log2FC_lps - log2FC_SAL_LPS,  # MOR LPS – normal LPS
    LPS_sig = (FDR_SAL_LPS < fdr_cut | FDR_lps < fdr_cut) &
              (abs(log2FC_SAL_LPS) > lfc_cut | abs(log2FC_lps) > lfc_cut),
    MOR_distorted_LPS = LPS_sig & abs(delta_LPS) > delta_cut
  )

merged4_sigLPS <- merged4 %>%
  filter(LPS_sig) %>%
  mutate(
    quadrant_D = case_when(
      log2FC_SAL_LPS >  0 & log2FC_lps >  0 ~ "hyper-activated (up→more up)",
      log2FC_SAL_LPS <  0 & log2FC_lps <  0 ~ "hyper-suppressed (down→more down)",
      log2FC_SAL_LPS >  0 & log2FC_lps <  0 ~ "inverted_up_to_down",
      log2FC_SAL_LPS <  0 & log2FC_lps >  0 ~ "inverted_down_to_up",
      TRUE                                  ~ "weak_or_no_LPS"
    )
  )

table(merged4_sigLPS$quadrant_D)

## ================== DISTORTION SCATTER + BARPLOT ==============

p_sig_scatter <- ggplot(merged4_sigLPS,
                        aes(x = log2FC_SAL_LPS,
                            y = log2FC_lps,
                            color = quadrant_D)) +
  geom_point(alpha = 0.8, size = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3) +
  coord_equal() +
  theme_classic(base_size = 10) +
  labs(
    x = "Normal LPS response (SAL–LPS vs SAL–SAL)",
    y = "LPS response after early-life morphine (MOR–LPS vs SAL–LPS)",
    title = "Significant LPS-responsive peaks\ncolored by morphine-induced distortion"
  )

count_quadrant <- merged4_sigLPS %>%
  as_tibble() %>%
  count(quadrant_D)

p_quadrant_bar <- ggplot(count_quadrant,
                         aes(x = fct_reorder(quadrant_D, n, .desc = TRUE),
                             y = n)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.2, size = 3) +
  theme_classic(base_size = 10) +
  labs(
    x = NULL,
    y = "Number of significant peaks",
    title = "Classes of morphine-altered LPS responses"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("hypo_H3K27ac_sigLPS_quadrant_scatter.pdf",
       p_sig_scatter, width = 5.5, height = 5)
ggsave("hypo_H3K27ac_sigLPS_quadrant_bar.pdf",
       p_quadrant_bar, width = 5.5, height = 4)

## ================== ANNOTATION ========

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

peak_gr <- GRanges(
  seqnames = merged4$seqnames,
  ranges   = IRanges(start = merged4$start, end = merged4$end),
  peak_id  = merged4$peak_id
)

anno <- annotatePeak(
  peak_gr,
  TxDb      = txdb,
  tssRegion = c(-1000, 1000),
  annoDb    = "org.Mm.eg.db"
)

anno_df <- as.data.frame(anno)

merged4_anno <- merged4 %>%
  left_join(anno_df %>% dplyr::select(peak_id, geneId, SYMBOL, annotation),
            by = "peak_id")

bg_genes <- merged4_anno %>%
  filter(!is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

merged4_anno_sig <- merged4_anno %>%
  left_join(merged4_sigLPS %>% select(peak_id, quadrant_D, MOR_distorted_LPS),
            by = "peak_id")

genes_baseline_MOR <- merged4_anno %>%
  filter(FDR_base < fdr_cut, abs(log2FC_base) > lfc_cut, !is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

genes_MOR_LPS <- merged4_anno %>%
  filter(FDR_lps < fdr_cut, abs(log2FC_lps) > lfc_cut, !is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

genes_MOR_distorted <- merged4_anno_sig %>%
  filter(MOR_distorted_LPS, !is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

genes_inverted <- merged4_anno_sig %>%
  filter(MOR_distorted_LPS,
         quadrant_D %in% c("inverted_up_to_down", "inverted_down_to_up"),
         !is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

genes_hyper <- merged4_anno_sig %>%
  filter(MOR_distorted_LPS,
         quadrant_D %in% c("hyper-activated (up→more up)",
                           "hyper-suppressed (down→more down)"),
         !is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

## ================== GREAT ANALYSES ===========================

# All distorted peaks
mor_distorted_LPS_peaks <- merged4_sigLPS %>%
  filter(MOR_distorted_LPS) %>%
  select(seqnames, start, end, peak_id)

gr_distorted <- make_gr(mor_distorted_LPS_peaks)

distorted_res <- run_great_clean(gr_distorted, "MOR_distorted")

# Baseline imprint (contrast 5)
gr_baseline <- get_sig_peaks_gr(db, contrast_idx = 5)
baseline_res <- run_great_clean(gr_baseline,
                                "baseline_MOR_SAL_vs_SAL_SAL")

# LPS response with morphine history (contrast 2)
gr_MOR_LPS <- get_sig_peaks_gr(db, contrast_idx = 2)
MOR_LPS_res <- run_great_clean(gr_MOR_LPS,
                               "MOR_LPS_vs_SAL_LPS")

# Distortion subclasses:
q_hyper_up   <- merged4_anno_sig %>%
  filter(quadrant_D == "hyper-activated (up→more up)")
q_hyper_down <- merged4_anno_sig %>%
  filter(quadrant_D == "hyper-suppressed (down→more down)")

gr_hyper_up   <- make_gr(q_hyper_up)
gr_hyper_down <- make_gr(q_hyper_down)

hyper_up_res   <- run_great_clean(gr_hyper_up,   "MOR_hyperUp")
hyper_down_res <- run_great_clean(gr_hyper_down, "MOR_hyperDown")

dist_up_to_down_peaks <- merged4_sigLPS %>%
  filter(MOR_distorted_LPS,
         quadrant_D == "inverted_up_to_down") %>%
  select(seqnames, start, end, peak_id) %>%
  distinct()

dist_down_to_up_peaks <- merged4_sigLPS %>%
  filter(MOR_distorted_LPS,
         quadrant_D == "inverted_down_to_up") %>%
  select(seqnames, start, end, peak_id) %>%
  distinct()

gr_dist_up_to_down   <- make_gr(dist_up_to_down_peaks)
gr_dist_down_to_up   <- make_gr(dist_down_to_up_peaks)

dist_up_to_down_res   <- run_great_clean(gr_dist_up_to_down,
                                         "MOR_distorted_upToDown")
dist_down_to_up_res   <- run_great_clean(gr_dist_down_to_up,
                                         "MOR_distorted_downToUp")

## ================== GREAT SUMMARY FIGURES ====================

# -------- 1) Summary across all peak sets --------------------
great_files <- c(
  "Baseline_MOR-SAL_vs_SAL-SAL"          = "GREAT_clean_baseline_MOR_SAL_vs_SAL_SAL_BP.tsv",
  "MOR-LPS_vs_SAL-LPS"                   = "GREAT_clean_MOR_LPS_vs_SAL_LPS_BP.tsv",
  "Distorted_all"                        = "GREAT_clean_MOR_distorted_BP.tsv",
  "Distorted_up→down"                    = "GREAT_clean_MOR_distorted_upToDown_BP.tsv",
  "Distorted_down→up"                    = "GREAT_clean_MOR_distorted_downToUp_BP.tsv",
  "Hyper_activated (up→more up)"         = "GREAT_clean_MOR_hyperUp_BP.tsv",
  "Hyper_suppressed (down→more down)"   = "GREAT_clean_MOR_hyperDown_BP.tsv"
)

great_all <- lapply(names(great_files), function(set_name) {
  fn <- great_files[[set_name]]
  if (!file.exists(fn)) {
    message("⚠️  File not found, skipping: ", fn)
    return(NULL)
  }
  df <- read_tsv(fn, show_col_types = FALSE)
  if (!nrow(df)) return(NULL)
  df$set <- set_name
  df
}) %>%
  bind_rows()

if (!nrow(great_all)) stop("No GREAT tables were loaded — check file names/paths.")

great_all <- great_all %>%
  mutate(
    minus_log10_p    = -log10(p_adjust_hyper),
    description_short = str_trunc(description, width = 60)
  )

n_top <- 20

great_top <- great_all %>%
  group_by(set) %>%
  slice_max(order_by = minus_log10_p, n = n_top, with_ties = FALSE) %>%
  ungroup()

set_levels <- c(
  "Baseline_MOR-SAL_vs_SAL-SAL",
  "MOR-LPS_vs_SAL-LPS",
  "Distorted_all",
  "Distorted_up→down",
  "Distorted_down→up",
  "Hyper_activated (up→more up)",
  "Hyper_suppressed (down→more down)"
)
great_top$set <- factor(great_top$set, levels = set_levels)

term_order <- great_top %>%
  group_by(description_short) %>%
  summarise(median_sig = median(minus_log10_p, na.rm = TRUE)) %>%
  arrange(desc(median_sig)) %>%
  pull(description_short)

great_top$description_short <- factor(great_top$description_short,
                                      levels = rev(term_order))

p_great_summary <- ggplot(
  great_top,
  aes(x = set,
      y = description_short,
      size = fold_enrichment_hyper,
      color = minus_log10_p)
) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(name = "Fold enrichment\n(hypergeometric)",
                        range = c(1.5, 6)) +
  scale_color_gradient(name = "-log10(adj. p)",
                       low = "grey80", high = "red3") +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank()
  ) +
  labs(
    title = "Summary of GREAT enrichment across morphine-sensitive peak sets",
    subtitle = "Top enriched GO:BP terms per H3K27ac peak class"
  )

ggsave("sup_fig_GREAT_summary.pdf",
       p_great_summary, width = 7.5, height = 11, useDingbats = FALSE)
ggsave("sup_fig_GREAT_summary.png",
       p_great_summary, width = 7.5, height = 11, dpi = 300)

# -------- 2) Shared terms across sets + categories ------------

PAT_IMMUNE <- paste(
  c("immune"," cytokine ","leukocyte"," T cell ","B cell","myeloid",
    "macrophage","microglia","interferon","chemokine",
    "interleukin","IL-","TNF","NF-κB","Nfkb","inflammatory",
    "inflammation","p38 MAPK","p38MAPK","Rho protein","Rap protein","Cdc42 protein"),
  collapse = "|"
)

PAT_METAB <- paste(
  c("lipid","glucose","insulin","thermogen","carbohydrate","triglyceride",
    "temperature","body weight","adipose","fat","energy","thermogenesis",
    " heat "," cold ","glucocorticoid","liposaccharide"),
  collapse = "|"
)

PAT_SLEEP <- paste(
  c("sleep","circadian","rhythmic"),
  collapse = "|"
)

great_shared_cat <- great_all %>%
  mutate(
    description_short = str_trunc(description, width = 60)
  ) %>%
  group_by(description_short) %>%
  filter(dplyr::n_distinct(set) >= 1) %>%  # keep if present in ≥1 set; categories do filtering
  ungroup() %>%
  mutate(
    Category = case_when(
      str_detect(description_short, regex(PAT_IMMUNE, ignore_case = TRUE)) ~ "Immune",
      str_detect(description_short, regex(PAT_METAB,  ignore_case = TRUE)) ~ "Metabolic",
      str_detect(description_short, regex(PAT_SLEEP,  ignore_case = TRUE)) ~ "Sleep / Circadian",
      TRUE ~ "Other"
    ),
    minus_log10_p = -log10(p_adjust_hyper)
  ) %>%
  filter(Category != "Other")

n_terms <- 50
top_shared_terms <- great_shared_cat %>%
  group_by(description_short) %>%
  summarise(median_sig = median(minus_log10_p, na.rm = TRUE)) %>%
  arrange(desc(median_sig)) %>%
  slice_head(n = n_terms) %>%
  pull(description_short)

great_shared_cat <- great_shared_cat %>%
  filter(description_short %in% top_shared_terms)

great_shared_cat$set <- factor(great_shared_cat$set, levels = set_levels)

term_order_cat <- great_shared_cat %>%
  group_by(Category, description_short) %>%
  summarise(median_sig = median(minus_log10_p, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(Category, desc(median_sig))

great_shared_cat$description_short <- factor(
  great_shared_cat$description_short,
  levels = rev(unique(term_order_cat$description_short))
)

p_great_shared_cat <- ggplot(
  great_shared_cat,
  aes(x = set,
      y = description_short,
      size  = fold_enrichment_hyper,
      color = minus_log10_p)
) +
  geom_point(alpha = 0.9) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  scale_size_continuous(name = "Fold enrichment\n(hypergeometric)",
                        range = c(1.5, 6)) +
  scale_color_gradient(name = "-log10(adj. p)",
                       low = "grey80", high = "red3") +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid   = element_blank(),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text.y = element_text(size = 9, face = "bold")
  ) +
  labs(
    title    = "Shared GREAT GO:BP terms across morphine-sensitive peak sets",
    subtitle = paste0(
      "GO terms grouped by category (Immune / Metabolic / Circadian); top ",
      n_terms, " shown"
    )
  )

ggsave("sup_fig_GREAT_shared_terms_byCategory.pdf",
       p_great_shared_cat, width = 7, height = 10, useDingbats = FALSE)
ggsave("sup_fig_GREAT_shared_terms_byCategory.png",
       p_great_shared_cat, width = 7, height = 10, dpi = 300)

message("Analysis complete.")
