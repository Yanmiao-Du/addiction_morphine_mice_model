#!/usr/bin/env Rscript
# 02_hypo_H3K27ac_DiffBind_DESeq2.R
#
# Hypothalamus H3K27ac DiffBind pipeline:
#  - load DiffBind metadata (diff4_input_hypo.csv)
#  - use consensus_peaks.bed, remove mm10 blacklist
#  - build count matrix and filter low-coverage peaks
#  - run DESeq2 contrasts for all 4 conditions
#  - export DiffBind results + annotated peaks (ChIPseeker)
#  - save RData objects for downstream figures (Fig2, Fig3)

library(DiffBind)
library(rtracklayer)
library(GenomicRanges)
# --- Load metadata and build DiffBind object ---
meta_file <- "diff4_input_hypo.csv"
dbObj <- dba(sampleSheet = meta_file)
dbObj$config$RunParallel <- FALSE
dbObj$config$AnalysisMethod <- DBA_DESEQ2
dbObj$config$DataType <- DBA_DATA_GRANGES
dbObj$config$Genome <- "mm10"

cat("Loaded", length(dbObj$samples$SampleID), "samples.\n")

# --- Load external consensus and blacklist ---
consensus_file <- "chip/consensus_peaks.bed"
mm10_blacklist_file <- "/stg3/data1/yanmiao/yanmiao2/addition_model_Mouse_Mor/hypothamus/Hypothalamus/mm10-blacklist.v2.bed"

consensus <- import(consensus_file)
blacklist <- import(mm10_blacklist_file)

cat("Loaded consensus peaks:", length(consensus), "regions.\n")
cat("Loaded mm10 blacklist regions:", length(blacklist), "regions.\n")

# --- Remove blacklist overlaps ---
consensus_clean <- subsetByOverlaps(consensus, blacklist, invert = TRUE)
cat("Removed", length(consensus) - length(consensus_clean), "blacklisted regions.\n")

# --- Count reads using cleaned consensus as the reference ---
dbObj <- dba.count(
  dbObj,
  peaks = consensus_clean,
  bUseSummarizeOverlaps = TRUE,
  filter = 1,
  summits = FALSE,
  bScaleControl = TRUE,
  score = DBA_SCORE_TMM_MINUS_FULL,
  minCount = 1
)
cat("Count matrix built successfully using blacklist-filtered consensus.\n")

save(dbObj, file = "hypo_diffbind_count_ready.RData")

dbObj <- get(load("/stg3/data1/yanmiao/yanmiao2/addition_model_Mouse_Mor/hypothamus/Hypothalamus/hypo_diffbind_count_ready.RData"))
count_matrix <- dbObj$binding
cat("Retrieved count matrix with", nrow(count_matrix), "peaks and", ncol(count_matrix), "samples.\n")
library(GenomicRanges)

# --- Convert all peaksets to proper GRanges ---
dbObj$peaks <- lapply(dbObj$peaks, function(p) {
  if (is.data.frame(p) && nrow(p) > 0) {
    names(p) <- tolower(names(p))  # lowercase all names for safety
    
    # Rename the key columns
    if ("chr" %in% names(p)) names(p)[names(p) == "chr"] <- "seqnames"
    if ("start" %in% names(p)) names(p)[names(p) == "start"] <- "start"
    if ("end" %in% names(p)) names(p)[names(p) == "end"] <- "end"
    
    # Ensure the essential columns are present
    if (!all(c("seqnames", "start", "end") %in% names(p))) {
      warning("⚠️ Skipping malformed peakset with missing seqnames/start/end.")
      return(NULL)
    }
    
    # Convert to GRanges
    GRanges(
      seqnames = p$seqnames,
      ranges   = IRanges(start = p$start, end = p$end),
      score    = if ("score" %in% names(p)) p$score else NA_real_
    )
  } else if (inherits(p, "GRanges")) {
    p
  } else {
    NULL
  }
})

# --- Drop NULL or empty entries ---
nonempty <- sapply(dbObj$peaks, function(p) !is.null(p) && length(p) > 0)
dbObj <- dbObj[nonempty]
cat("Converted", sum(nonempty), "peaksets to GRanges successfully.\n")

# --- Verify ---
print(class(dbObj$peaks[[1]]))
print(dbObj$peaks[[1]][1:3])
# --- Ensure all peaksets are valid before filtering ---
nonempty <- sapply(dbObj$peaks, function(p) !is.null(p) && length(p) > 0)
cat("Found", sum(nonempty), "non-empty peaksets out of", length(dbObj$peaks), "\n")

dbObj_clean <- dbObj
dbObj_clean$peaks <- dbObj$peaks[nonempty]

# --- Apply high-coverage filtering only to nonempty GRanges ---
dbObj_filtered <- dbObj_clean
dbObj_filtered$peaks <- lapply(dbObj_clean$peaks, function(p) {
  if (inherits(p, "GRanges")) {
    p[keep_idx]
  } else {
    NULL
  }
})

cat("Filtered low-coverage peaks successfully.\n")

save(dbObj_filtered, file = "/stg3/data1/yanmiao/yanmiao2/addition_model_Mouse_Mor/hypothamus/Hypothalamus/hypo_diffbind_filtered.RData")
meta <- read.csv("diff4_input_hypo.csv")
consensus_filtered <- dbObj_filtered$peaks[[1]] 
# Keep peaks with >=20 reads in ≥50% samples
keep_idx <- rowSums(count_matrix >= 20) >= (0.5 * ncol(count_matrix))
dbObj_filtered <- dbObj
dbObj_filtered$peaks <- lapply(dbObj$peaks, function(p) p[keep_idx])

cat("Filtered peaks with low coverage (", sum(keep_idx), "kept).\n")

library(DiffBind)
library(GenomicRanges)
library(rtracklayer)

# Reload metadata and filtered GRanges ------------------------------------
meta <- read.csv("diff4_input_hypo.csv")

# Retrieve the filtered GRanges peakset (from your filtering step)
consensus_filtered <- dbObj_filtered$peaks[[1]]  # they are all the same consensus

# Sanity check
stopifnot(inherits(consensus_filtered, "GRanges"))
cat("Using filtered consensus of", length(consensus_filtered), "regions.\n")

# Rebuild new DiffBind object ---------------------------------------------
dbObj_clean <- dba(sampleSheet = meta)

dbObj_clean <- dba.count(
  dbObj_clean,
  peaks = consensus_filtered,
  bUseSummarizeOverlaps = TRUE,
  filter = 1,
  summits = FALSE,
  bScaleControl = TRUE,
  score = DBA_SCORE_TMM_MINUS_FULL,
  minCount = 1
)

cat("Count matrix successfully rebuilt.\n")

# Save clean object -------------------------------------------------------
save(dbObj_clean, file = "hypo_diffbind_clean_ready.RData")
dbObj_clean <- dba.contrast(dbObj_clean, categories = DBA_CONDITION)
dbObj_clean <- dba.analyze(dbObj_clean, method = DBA_DESEQ2)

cat("Differential analysis complete.\n")

# ---------------------------------------------------------
#  Step 1. Verify design and contrast order
# ---------------------------------------------------------
dba.show(dbObj_clean, bContrasts = TRUE)
# (should show 6 contrasts with Condition as factor)

# ---------------------------------------------------------
#  Step 2. Extract each relevant comparison
# ---------------------------------------------------------

# MOR-LPS vs SAL-LPS  → morphine effect under immune challenge
mor_lps_vs_sal_lps <- dba.report(dbObj_clean,
                                 contrast = 2,  # from your table
                                 th = 0.05)

# MOR-LPS vs MOR-SAL  → morphine × LPS interaction
mor_lps_vs_mor_sal <- dba.report(dbObj_clean,
                                 contrast = 1,
                                 th = 0.05)

# MOR-SAL vs SAL-SAL  → baseline morphine effect
mor_sal_vs_sal_sal <- dba.report(dbObj_clean,
                                 contrast = 5,
                                 th = 0.05)

# SAL-LPS vs SAL-SAL  → immune response without morphine
sal_lps_vs_sal_sal <- dba.report(dbObj_clean,
                                 contrast = 6,
                                 th = 0.05)

# ---------------------------------------------------------
#  Step 3. Export to tab-delimited files
# ---------------------------------------------------------
write.table(as.data.frame(mor_lps_vs_sal_lps),
            "DiffBind_MORLPS_vs_SALLPS.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(mor_lps_vs_mor_sal),
            "DiffBind_MORLPS_vs_MORSAL.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(mor_sal_vs_sal_sal),
            "DiffBind_MORSAL_vs_SALSAL.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(as.data.frame(sal_lps_vs_sal_sal),
            "DiffBind_SALLPS_vs_SALSAL.tsv", sep="\t", quote=FALSE, row.names=FALSE)


save(dbObj_clean,
     mor_lps_vs_sal_lps,
     mor_lps_vs_mor_sal,
     mor_sal_vs_sal_sal,
     sal_lps_vs_sal_sal,
     file = "hypo_diffbind_DESeq2_allContrasts.RData")

# ============================================================
# Step 1. Load required packages
# ============================================================
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

# ============================================================
# Step 2. Annotate each contrast
# ============================================================

cat("Annotating MOR-LPS vs SAL-LPS...\n")
anno_mor_lps_vs_sal_lps <- annotatePeak(
  mor_lps_vs_sal_lps,
  TxDb   = TxDb.Mmusculus.UCSC.mm10.knownGene,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Mm.eg.db"
)

cat("Annotating MOR-LPS vs MOR-SAL...\n")
anno_mor_lps_vs_mor_sal <- annotatePeak(
  mor_lps_vs_mor_sal,
  TxDb   = TxDb.Mmusculus.UCSC.mm10.knownGene,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Mm.eg.db"
)

cat("Annotating MOR-SAL vs SAL-SAL...\n")
anno_mor_sal_vs_sal_sal <- annotatePeak(
  mor_sal_vs_sal_sal,
  TxDb   = TxDb.Mmusculus.UCSC.mm10.knownGene,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Mm.eg.db"
)

cat("Annotating SAL-LPS vs SAL-SAL...\n")
anno_sal_lps_vs_sal_sal <- annotatePeak(
  sal_lps_vs_sal_sal,
  TxDb   = TxDb.Mmusculus.UCSC.mm10.knownGene,
  tssRegion = c(-1000, 1000),
  annoDb = "org.Mm.eg.db"
)

cat("✅ All four contrasts annotated.\n")

# ============================================================
# Step 3. Export annotation tables
# ============================================================
write.table(as.data.frame(anno_mor_lps_vs_sal_lps),
            "Annotated_MORLPS_vs_SALLPS.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(anno_mor_lps_vs_mor_sal),
            "Annotated_MORLPS_vs_MORSAL.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(anno_mor_sal_vs_sal_sal),
            "Annotated_MORSAL_vs_SALSAL.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(as.data.frame(anno_sal_lps_vs_sal_sal),
            "Annotated_SALLPS_vs_SALSAL.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


# ============================================================
# Step 4. Save everything for reuse
# ============================================================
save(
  anno_mor_lps_vs_sal_lps,
  anno_mor_lps_vs_mor_sal,
  anno_mor_sal_vs_sal_sal,
  anno_sal_lps_vs_sal_sal,
  file = "hypo_diffbind_annotated_all.RData"
)