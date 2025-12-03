suppressPackageStartupMessages({
  library(DESeq2)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(IRanges)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(tibble)
  library(readr)
  library(stringr)
})

ROOT <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo"

## ============================================================
## PATTERN DEFINITIONS (used everywhere, keep in sync)
## ============================================================
PAT_IMMUNE <- paste(
  c("immune","cytokine","leukocyte"," T cell ","B cell","myeloid",
    "macrophage","microglia","interferon","chemokine",
    "interleukin","IL-","TNF","NF-κB","Nfkb","inflammatory",
    "inflammation"),
  collapse = "|"
)

PAT_METAB <- paste(
  c("lipid","glucose","insulin","thermogen","carbohydrate","triglyceride"
    "temperature","body weight","adipose","fat","energy","thermogenesis","heat","cold"),
  collapse = "|"
)

PAT_SLEEP <- paste(
  c("sleep","circadian","rhythmic"),
  collapse = "|"
)

## ============================================================
## 1) Load DESeq2 object & compute TPM from raw counts
## ============================================================
dds <- readRDS(file.path(ROOT, "hypo_deg_qcpass_dds.rds"))

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
ex_by_gene <- exonsBy(txdb, by = "gene")  # GRangesList

gene_len_vec <- sapply(ex_by_gene, function(gr) {
  sum(width(IRanges::reduce(gr)))
})

gene_len_df <- tibble(
  ENTREZID  = names(gene_len_vec),
  length_bp = as.numeric(gene_len_vec)
)

map_df <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = gene_len_df$ENTREZID,
  keytype = "ENTREZID",
  columns = "SYMBOL"
)

gene_len_sym <- gene_len_df %>%
  left_join(map_df, by = "ENTREZID") %>%
  filter(!is.na(SYMBOL)) %>%
  distinct(SYMBOL, .keep_all = TRUE)

idx <- match(rownames(dds), gene_len_sym$SYMBOL)
rowData(dds)$gene_length <- gene_len_sym$length_bp[idx]

# raw counts
cts <- counts(dds, normalized = FALSE)
gene_len_kb <- rowData(dds)$gene_length / 1000

keep <- !is.na(gene_len_kb) & gene_len_kb > 0
cts_keep <- cts[keep, ]
len_keep <- gene_len_kb[keep]

# RPK and TPM
rpk <- sweep(cts_keep, 1, len_keep, "/")
scaling_factor <- colSums(rpk)
tpm <- sweep(rpk, 2, scaling_factor, "/") * 1e6

tpm_df <- tpm %>%
  as.data.frame() %>%
  rownames_to_column("gene")

## ============================================================
## 2) Sample metadata with 4 conditions
## ============================================================
sample_info <- as.data.frame(colData(dds)) %>%
  rownames_to_column("sample") %>%
  mutate(
    Condition = case_when(
      grepl("SAL[-_.]SAL", Group) ~ "SAL-SAL",
      grepl("MOR[-_.]SAL", Group) ~ "MOR-SAL",
      grepl("SAL[-_.]LPS", Group) ~ "SAL-LPS",
      grepl("MOR[-_.]LPS", Group) ~ "MOR-LPS",
      TRUE                        ~ NA_character_
    ),
    Condition = factor(
      Condition,
      levels = c("SAL-SAL","MOR-SAL","SAL-LPS","MOR-LPS")
    )
  )

save.image(file.path(ROOT, "hypo_deg_qcpass_dds_withTPM_andGeneLength.RData"))

## ============================================================
## 3) MOR–SAL vs SAL–SAL DEGs & GO
## ============================================================
deg_morsal <- read_tsv(
  file.path(ROOT, "DEGs_treatment_only/SAL_SAL_vs_MOR_SAL_sig.tsv"),
  show_col_types = FALSE
) %>%
  rename(gene = gene_name, FDR = padj)

go_down_deg <- read_tsv(
  file.path(ROOT, "SAL_SAL_vs_MOR_SAL_DEG_down_GO_BP_enrichment.tsv"),
  show_col_types = FALSE
)
go_up_deg <- read_tsv(
  file.path(ROOT, "SAL_SAL_vs_MOR_SAL_DEG_up_GO_BP_enrichment.tsv"),
  show_col_types = FALSE
)

go_all_deg <- bind_rows(
  mutate(go_down_deg, Direction = "Down"),
  mutate(go_up_deg,   Direction = "Up")
)

## ============================================================
## 4) Helpers: GO → genes, DEG selection
## ============================================================
extract_genes_from_GO <- function(go_tbl, pattern) {
  go_tbl %>%
    filter(str_detect(Description, regex(pattern, ignore_case = TRUE))) %>%
    separate_rows(geneID, sep = "/") %>%
    pull(geneID) %>%
    unique()
}

select_theme_degs <- function(genes, deg_tbl) {
  deg_tbl %>%
    filter(FDR < 0.05, gene %in% genes) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    pull(gene) %>%
    unique()
}

## Theme-specific gene pools from MOR–SAL DEG GO
immune_genes_all    <- extract_genes_from_GO(go_all_deg, PAT_IMMUNE)
metabolic_genes_all <- extract_genes_from_GO(go_all_deg, PAT_METAB)
sleep_genes_all     <- extract_genes_from_GO(go_all_deg, PAT_SLEEP)

genes_immune_MORSAL    <- select_theme_degs(immune_genes_all,    deg_morsal)
genes_metabolic_MORSAL <- select_theme_degs(metabolic_genes_all, deg_morsal)
genes_sleep_MORSAL     <- select_theme_degs(sleep_genes_all,     deg_morsal)

## ============================================================
## 5) MOR–LPS vs SAL–LPS GO (diff-peaks) → genes
## ============================================================
go_diff_LPS_enhancer <- read_tsv(
  file.path(ROOT, "Fig3F_GO_BP_Enhancer_MOR-LPS_vs_SAL-LPS.tsv"),
  show_col_types = FALSE
)
go_diff_LPS_promoter <- read_tsv(
  file.path(ROOT, "Fig3F_GO_BP_Promoter_MOR-LPS_vs_SAL-LPS.tsv"),
  show_col_types = FALSE
)

go_diff_LPS <- bind_rows(go_diff_LPS_enhancer, go_diff_LPS_promoter)

immune_genes_LPS_all    <- extract_genes_from_GO(go_diff_LPS, PAT_IMMUNE)
metabolic_genes_LPS_all <- extract_genes_from_GO(go_diff_LPS, PAT_METAB)
sleep_genes_LPS_all     <- extract_genes_from_GO(go_diff_LPS, PAT_SLEEP)

# For LPS, we just require the gene to exist in TPM
genes_immune_LPS    <- intersect(immune_genes_LPS_all,    tpm_df$gene)
genes_metabolic_LPS <- intersect(metabolic_genes_LPS_all, tpm_df$gene)
genes_sleep_LPS     <- intersect(sleep_genes_LPS_all,     tpm_df$gene)

## ============================================================
## 6) TPM boxplots (pairwise, gene+P in strip, auto height)
## ============================================================
plot_tpm_box_pair <- function(genes, cond_pair, theme_label, out_prefix,
                              ncol_facets = 4) {
  stopifnot(length(cond_pair) == 2)

  keep_genes <- intersect(genes, tpm_df$gene)
  if (length(keep_genes) == 0L) {
    message("No genes for ", out_prefix, " present in TPM matrix.")
    return(invisible(NULL))
  }

  meta_sub <- sample_info %>%
    filter(Condition %in% cond_pair)

  dat_long <- tpm_df %>%
    filter(gene %in% keep_genes) %>%
    pivot_longer(-gene, names_to = "sample", values_to = "TPM") %>%
    inner_join(meta_sub, by = "sample") %>%
    mutate(
      Condition2 = factor(Condition, levels = cond_pair)
    )

  ## per-gene P-values (Wilcoxon)
  pvals <- dat_long %>%
    group_by(gene) %>%
    summarise(
      p = tryCatch(
        wilcox.test(TPM ~ Condition2)$p.value,
        error = function(e) NA_real_
      ),
      .groups = "drop"
    ) %>%
    mutate(
      p_label = ifelse(is.na(p), "n.s.", paste0("P=", signif(p, 2)))
    )

  ## keep only significant genes (P <= 0.05)
  pvals_sig <- pvals %>%
    filter(!is.na(p) & p <= 0.05)

  if (nrow(pvals_sig) == 0L) {
    message("No genes with P <= 0.05 for ", out_prefix, "; skipping plot.")
    return(invisible(NULL))
  }

  ## join back only sig genes and build facet label gene + P
  dat_long <- dat_long %>%
    inner_join(pvals_sig, by = "gene") %>%
    mutate(gene_p = paste0(gene, "\n", p_label))

  n_panels <- length(unique(dat_long$gene_p))
  n_rows   <- ceiling(n_panels / ncol_facets)
  fig_h    <- max(2 * n_rows, 3)

  p <- ggplot(dat_long,
              aes(x = Condition2, y = TPM, fill = Condition2)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1.2, alpha = 0.7, colour = "black") +
    scale_y_log10() +
    scale_fill_manual(values = c(
      "SAL-SAL" = "#999999",
      "MOR-SAL" = "#E64B35",
      "SAL-LPS" = "#999999",
      "MOR-LPS" = "#E64B35"
    )) +
    facet_wrap(~ gene_p, scales = "free_y", ncol = ncol_facets) +
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(size = 8),
      plot.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = paste0("Supplementary: TPM of ", theme_label,
                     " (", cond_pair[2], " vs ", cond_pair[1], ")"),
      x     = NULL,
      y     = "TPM (log10)"
    )

  ggsave(file.path(ROOT, paste0(out_prefix, "_box.pdf")),
         p, width = 12, height = fig_h)
  ggsave(file.path(ROOT, paste0(out_prefix, "_box.png")),
         p, width = 12, height = fig_h, dpi = 600)

  p
}


## ---- Boxplot calls ----
# MOR–SAL vs SAL–SAL (DEGs)
plot_tpm_box_pair(
  genes       = genes_immune_MORSAL,
  cond_pair   = c("SAL-SAL","MOR-SAL"),
  theme_label = "immune-related MOR–SAL vs SAL–SAL DEGs",
  out_prefix  = "Supp_Fig2_immune_MORSAL"
)

plot_tpm_box_pair(
  genes       = genes_metabolic_MORSAL,
  cond_pair   = c("SAL-SAL","MOR-SAL"),
  theme_label = "metabolic/body-weight MOR–SAL vs SAL–SAL DEGs",
  out_prefix  = "Supp_Fig2_metabolic_MORSAL"
)

plot_tpm_box_pair(
  genes       = genes_sleep_MORSAL,
  cond_pair   = c("SAL-SAL","MOR-SAL"),
  theme_label = "sleep/circadian MOR–SAL vs SAL–SAL DEGs",
  out_prefix  = "Supp_Fig2_sleep_MORSAL"
)

# # MOR–LPS vs SAL–LPS (diff-peak genes)
# plot_tpm_box_pair(
#   genes       = genes_immune_LPS,
#   cond_pair   = c("SAL-LPS","MOR-LPS"),
#   theme_label = "immune MOR–LPS vs SAL–LPS diff-peak genes",
#   out_prefix  = "Supp_Fig3_immune_MORLPS"
# )

# plot_tpm_box_pair(
#   genes       = genes_metabolic_LPS,
#   cond_pair   = c("SAL-LPS","MOR-LPS"),
#   theme_label = "metabolic/thermogenic MOR–LPS vs SAL–LPS genes",
#   out_prefix  = "Supp_Fig3_metabolic_MORLPS"
# )

# plot_tpm_box_pair(
#   genes       = genes_sleep_LPS,
#   cond_pair   = c("SAL-LPS","MOR-LPS"),
#   theme_label = "sleep/circadian MOR–LPS vs SAL–LPS genes",
#   out_prefix  = "Supp_Fig3_sleep_MORLPS"
# )

## ============================================================
## 7) Pairwise TPM heatmaps (no col clustering, pair only)
## ============================================================
plot_tpm_heatmap_pair <- function(genes, cond_pair, title, out_prefix) {

  stopifnot(length(cond_pair) == 2)

  keep_genes <- intersect(genes, tpm_df$gene)
  if (length(keep_genes) == 0L) {
    message("No genes for ", out_prefix, " present in TPM matrix.")
    return(invisible(NULL))
  }

  mat_tpm <- tpm_df %>%
    filter(gene %in% keep_genes) %>%
    column_to_rownames("gene") %>%
    as.matrix()

  meta_sub <- sample_info %>%
    filter(Condition %in% cond_pair) %>%
    mutate(Condition = droplevels(Condition)) %>%
    arrange(Condition, sample)

  sample_order <- meta_sub$sample
  sample_order <- sample_order[sample_order %in% colnames(mat_tpm)]

  if (length(sample_order) == 0L) {
    message("No samples for ", paste(cond_pair, collapse = " vs "),
            " in TPM matrix; skipping heatmap.")
    return(invisible(NULL))
  }

  mat_tpm <- mat_tpm[, sample_order, drop = FALSE]
  mat_log <- log2(mat_tpm + 1)

  annot_col <- meta_sub %>%
    filter(sample %in% sample_order) %>%
    dplyr::select(sample, Condition)
  rownames(annot_col) <- annot_col$sample
  annot_col$sample <- NULL
  annot_col$Condition <- droplevels(annot_col$Condition)

  stopifnot(identical(colnames(mat_log), rownames(annot_col)))

  cond_levels <- levels(annot_col$Condition)
  base_cols <- c(
    "SAL-SAL" = "#D55E00",
    "MOR-SAL" = "#0072B2",
    "SAL-LPS" = "#009E73",
    "MOR-LPS" = "#E69F00"
  )

  ann_cols <- list(Condition = base_cols[cond_levels])

  pdf_width  <- 3 + 0.1 * ncol(mat_log)
  pdf_height <- 2 + 0.1 * nrow(mat_log)

  pdf(file.path(ROOT, paste0(out_prefix, "_heatmap.pdf")),
      width = pdf_width, height = pdf_height)

  pheatmap(
    mat_log,
    annotation_col    = annot_col,
    annotation_colors = ann_cols,
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,   # keep controls then MOR, no reordering
    scale             = "row",
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    color             = colorRampPalette(c("navy","white","firebrick3"))(50),
    main              = title
  )

  dev.off()
}

## Heatmap calls – MOR–SAL vs SAL–SAL
plot_tpm_heatmap_pair(
  genes      = genes_immune_MORSAL,
  cond_pair  = c("SAL-SAL","MOR-SAL"),
  title      = "Immune MOR–SAL vs SAL–SAL (log2 TPM+1, row z-score)",
  out_prefix = "Supp_Fig2_immune_MORSAL"
)

plot_tpm_heatmap_pair(
  genes      = genes_metabolic_MORSAL,
  cond_pair  = c("SAL-SAL","MOR-SAL"),
  title      = "Metabolic/body-weight MOR–SAL vs SAL–SAL (log2 TPM+1, row z-score)",
  out_prefix = "Supp_Fig2_metabolic_MORSAL"
)

plot_tpm_heatmap_pair(
  genes      = genes_sleep_MORSAL,
  cond_pair  = c("SAL-SAL","MOR-SAL"),
  title      = "Sleep/circadian MOR–SAL vs SAL–SAL (log2 TPM+1, row z-score)",
  out_prefix = "Supp_Fig2_sleep_MORSAL"
)

# ## Heatmap calls – MOR–LPS vs SAL–LPS
# plot_tpm_heatmap_pair(
#   genes      = genes_immune_LPS,
#   cond_pair  = c("SAL-LPS","MOR-LPS"),
#   title      = "Immune MOR–LPS vs SAL–LPS (log2 TPM+1, row z-score)",
#   out_prefix = "Supp_Fig3_immune_MORLPS"
# )

# plot_tpm_heatmap_pair(
#   genes      = genes_metabolic_LPS,
#   cond_pair  = c("SAL-LPS","MOR-LPS"),
#   title      = "Metabolic/thermogenic MOR–LPS vs SAL–LPS (log2 TPM+1, row z-score)",
#   out_prefix = "Supp_Fig3_metabolic_MORLPS"
# )

# plot_tpm_heatmap_pair(
#   genes      = genes_sleep_LPS,
#   cond_pair  = c("SAL-LPS","MOR-LPS"),
#   title      = "Sleep/circadian MOR–LPS vs SAL–LPS (log2 TPM+1, row z-score)",
#   out_prefix = "Supp_Fig3_sleep_MORLPS"
# )

## ============================================================
## 8) GO barplots colored by category (Immune/Metabolic/Sleep)
## ============================================================
annot_go_category <- function(go_tbl) {
  go_tbl %>%
    mutate(
      Category = case_when(
        str_detect(Description, regex(PAT_IMMUNE, ignore_case = TRUE)) ~ "Immune",
        str_detect(Description, regex(PAT_METAB,  ignore_case = TRUE)) ~ "Metabolic/thermo",
        str_detect(Description, regex(PAT_SLEEP,  ignore_case = TRUE)) ~ "Circadian",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(!is.na(Category))
}

plot_go_bar_category <- function(go_tbl, title, out_prefix,
                                 facet_by_category = TRUE) {

  go_cat <- annot_go_category(go_tbl) %>%
    mutate(minus_log10_FDR = -log10(p.adjust))

  if (nrow(go_cat) == 0L) {
    message("No GO terms with categories for ", out_prefix)
    return(invisible(NULL))
  }

  go_cat <- go_cat %>%
    mutate(
      Category    = factor(Category,
                           levels = c("Immune","Metabolic/thermo","Circadian")),
      Description = factor(Description,
                           levels = rev(unique(Description)))
    )

  p <- ggplot(go_cat,
              aes(x = minus_log10_FDR, y = Description, fill = Category)) +
    geom_col() +
    theme_bw(base_size = 11) +
    theme(
      axis.title.y = element_blank(),
      plot.title   = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = title,
      x     = expression(-log[10]("FDR"))
    )

  if (facet_by_category) {
    p <- p + facet_grid(Category ~ ., scales = "free_y", space = "free_y")
  }

  h <- max(3, 0.2 * nrow(go_cat))

  ggsave(file.path(ROOT, paste0(out_prefix, "_GObar.pdf")),
         p, width = 10, height = h)
  ggsave(file.path(ROOT, paste0(out_prefix, "_GObar.png")),
         p, width = 10, height = h, dpi = 600)

  p
}

## MOR–SAL (DEG-based GO)
plot_go_bar_category(
  go_tbl    = go_all_deg,
  title     = "GO categories for MOR–SAL vs SAL–SAL DEGs",
  out_prefix= "Supp_Fig2_MORSAL_GO_categories"
)

## MOR–LPS (diff-peak GO)
plot_go_bar_category(
  go_tbl    = go_diff_LPS,
  title     = "GO categories for MOR–LPS vs SAL–LPS diff-peak genes",
  out_prefix= "Supp_Fig3_MORLPS_GO_categories"
)

## ============================================================
## 8b) Alternative GO dot plots (like cytokine panel)
##      – size = gene count, colour = -log10(FDR)
## ============================================================

plot_go_dot_category <- function(go_tbl,
                                 category_filter  = "Immune",
                                 direction_filter = NULL,   # e.g. "Down"
                                 top_n            = 10,
                                 title,
                                 out_prefix,
                                 collapse_prefix  = 60) {   # <- NEW

  # 1) annotate categories with same patterns as bar plot
  go_cat <- annot_go_category(go_tbl)

  if (nrow(go_cat) == 0L) {
    message("No GO terms with categories for ", out_prefix)
    return(invisible(NULL))
  }

  # 2) filter by category
  go_cat <- go_cat %>%
    dplyr::filter(Category == category_filter) %>%
    dplyr::mutate(
      minus_log10_FDR = -log10(p.adjust),
      Category        = factor(Category,
                               levels = c("Immune","Metabolic/thermo","Circadian"))
    )

  # 3) optional filter by direction (robust / fuzzy match)
  if (!is.null(direction_filter) && "Direction" %in% colnames(go_cat)) {
    go_cat <- go_cat %>%
      dplyr::filter(stringr::str_detect(
        Direction,
        stringr::regex(direction_filter, ignore_case = TRUE)
      ))
  }

  if (nrow(go_cat) == 0L) {
    message("No GO terms after filtering for ", out_prefix)
    return(invisible(NULL))
  }

  # 4) "simplify": collapse very similar descriptions by prefix
  #    (e.g. all those CD4+ alpha-beta T-cell terms)
  go_cat_simplified <- go_cat %>%
    dplyr::mutate(
      Desc_prefix = stringr::str_sub(Description, 1, collapse_prefix)
    ) %>%
    dplyr::group_by(Desc_prefix) %>%
    dplyr::slice_min(p.adjust, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  # 5) take top_n of these simplified terms
  go_top <- go_cat_simplified %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(
      Description = factor(Description,
                           levels = rev(unique(Description)))
    )

  if (nrow(go_top) == 0L) {
    message("No GO terms left after simplification for ", out_prefix)
    return(invisible(NULL))
  }

  max_count <- max(go_top$Count, na.rm = TRUE)

  p <- ggplot(go_top,
              aes(x = Count, y = Description)) +
    geom_point(aes(size = Count, colour = minus_log10_FDR)) +
    scale_size_continuous(range = c(3, 8)) +
    scale_colour_gradient(low = "#91bfdb", high = "#fc8d59") +
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.25)),
      limits = c(1, max_count * 1.1)
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.title.y = element_blank(),
      plot.title   = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9)
    ) +
    labs(
      title  = title,
      x      = "Gene count",
      colour = expression(-log[10]("FDR")),
      size   = "Gene count"
    )

  h <- max(3, 0.35 * nrow(go_top))

  ggsave(file.path(ROOT, paste0(out_prefix, "_GOdot.pdf")),
         p, width = 8, height = h)
  ggsave(file.path(ROOT, paste0(out_prefix, "_GOdot.png")),
         p, width = 8, height = h, dpi = 600)

  p
}

## ============================
## Example calls
## ============================

## 1) MOR–SAL vs SAL–SAL: immune DOWN only (like your example)
##    (go_all_deg already has Direction = "Down"/"Up")
plot_go_dot_category(
  go_tbl          = go_all_deg,
  category_filter = "Immune",
  top_n           = 15,
  title           = "Top downregulated immune GO terms (MOR–SAL vs SAL–SAL))))",
  out_prefix      = "Supp_Fig2_MORSAL_ImmuneDown"
)

## 2) MOR–SAL vs SAL–SAL: metabolic / circadian (no direction filter)
plot_go_dot_category(
  go_tbl          = go_all_deg,
  category_filter = "Metabolic/thermo",
  top_n           = 10,
  title           = "Metabolic/thermogenic GO terms (MOR–SAL vs SAL–SAL)",
  out_prefix      = "Supp_Fig2_MORSAL_Metabolic"
)

plot_go_dot_category(
  go_tbl          = go_all_deg,
  category_filter = "Circadian",
  top_n           = 10,
  title           = "Circadian GO terms (MOR–SAL vs SAL–SAL)",
  out_prefix      = "Supp_Fig2_MORSAL_Circadian"
)

## 3) MOR–LPS vs SAL–LPS diff-peak GO (no Direction column)
plot_go_dot_category(
  go_tbl          = go_diff_LPS,
  category_filter = "Immune",
  top_n           = 10,
  title           = "Immune GO terms (MOR–LPS vs SAL–LPS diff-peak genes)",
  out_prefix      = "Supp_Fig3_MORLPS_Immune"
)

plot_go_dot_category(
  go_tbl          = go_diff_LPS,
  category_filter = "Metabolic/thermo",
  top_n           = 10,
  title           = "Metabolic/thermogenic GO terms (MOR–LPS vs SAL–LPS diff-peak genes)",
  out_prefix      = "Supp_Fig3_MORLPS_Metabolic"
)

plot_go_dot_category(
  go_tbl          = go_diff_LPS,
  category_filter = "Circadian",
  top_n           = 10,
  title           = "Circadian GO terms (MOR–LPS vs SAL–LPS diff-peak genes)",
  out_prefix      = "Supp_Fig3_MORLPS_Circadian"
)