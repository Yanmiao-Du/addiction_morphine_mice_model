suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
# suppressMessages(library(xlsx))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(org.Mm.eg.db))
suppressMessages(library(clusterProfiler))
suppressMessages(library(AnnotationDbi))
suppressMessages(library(openxlsx))
# suppressMessages(library(readxl))

maindir <- "/Users/yanmiaodu/Downloads/MOR_SAL/hypo/hypo_Taiji_heatmaps"
if (!dir.exists(maindir)) {
  dir.create(maindir, recursive = TRUE)
}
setwd("/Users/yanmiaodu/Downloads/MOR_SAL/hypo/hypo_Taiji_heatmaps")
set.seed(42)
fl.sources <- list.files("/Users/yanmiaodu/Downloads/MOR_SAL/Taiji2/scripts/utils", full.names = T)
tmp <- sapply(fl.sources,source)


Data <- read.csv("/Users/yanmiaodu/Downloads/MOR_SAL/hypo/merged_taiji_generanks.csv", header = TRUE, row.names = 1)
head(Data)
# remove SAL.LPS.M_1.4
Data <- Data[, !(colnames(Data) == "SAL.LPS.M_1.4_output")]
colnames(Data) <- gsub("_output$", "", colnames(Data))


first_upper <- function(string) {
  sapply(string, function(x) {
    paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
  })
}
row.names(Data) <- first_upper(rownames(Data))

meta <- as.data.frame(t(as.data.frame(strsplit(colnames(Data), "_"))))
colnames(meta) <- c("group", "id")
meta$group_no_sex <- sub("\\.F$|\\.M$", "", meta$group)
split_values <- as.data.frame(do.call(rbind, strsplit(as.character(meta$group), "\\.")))
colnames(split_values) <- c("condition", "treatment", "sex")

# Add the split values to the group_sorted dataframe
meta$condition <- split_values$condition
meta$treatment <- split_values$treatment
meta$sex <- split_values$sex

# Set row names
row.names(meta) <- colnames(Data)
meta <- meta[, !colnames(meta) %in% c("NA")]
#only plot group_no_gender, condition and treatment 
mycolors <- list(
    group = c("#A6A6A6", "#FF0000", "#0007F4", "#FFA75C", "#00AAFE", 
                             "#009051", "#FEB5B5", "#E6CFFF"),
  group_no_sex = c("MOR.LPS" = "#FF9999", "MOR.SAL" = "#FFCC99", "SAL.LPS" = "#99CCFF", "SAL.SAL" = "#99FF99"),
  sex = c("F" = "#9370DB", "M" = "#4682B4"),
  condition = c("MOR" = "#FFD700", "SAL" = "#CD5C5C"),
  treatment = c("LPS" = "#8FBC8F", "SAL" = "#F08080")
)
names(mycolors$group) <- unique(meta$group)

group_sorted <- meta[, colnames(meta) %in% c("group_no_sex", "condition","treatment","id")]

#heatmap func #plot heatmap - House keeping genes
breaksList = seq(-2, 2, by = 0.1)
# function
get_heatmap <- function(df, size = 3, ...) {
  # Ensure column names match the row names of group_sorted
  g2 <- group_sorted[match(colnames(df), rownames(group_sorted)),]
  
  # Create annotation for the columns
  annotation_col <- g2
  rownames(annotation_col) <- colnames(df)
  
  # Generate the heatmap
  p1 <- pheatmap(df, 
                 fontsize = size, 
                 angle_col = 90, 
                 cellwidth = 2 * size, 
                 cellheight = size,
                 cluster_cols = FALSE, 
                 annotation_col = annotation_col,
                 annotation_colors = mycolors, 
                 clustering_distance_cols = 'correlation', 
                 clustering_distance_rows = 'correlation', 
                 clustering_method = 'average',
                 breaks = breaksList,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
                 annotation_names_row = FALSE, 
                 annotation_names_col = FALSE,
                 border_color = NA, ...)
  
  return(p1)
}

get_clean_heatmap <- function(df, size = 3, ...) {
  g2 <- group_sorted[match(colnames(df), rownames(group_sorted)), ]
  
  # Ensure this is a proper data frame with rownames
  annotation_col <- data.frame(group_no_sex = g2$group_no_sex)
  rownames(annotation_col) <- colnames(df)

  p1 <- pheatmap(df, 
  cluster_cols = FALSE,
                 fontsize = size, 
                 angle_col = 90, 
                 cellwidth = 8, 
                 cellheight = 4,
                 # Don't define cluster_rows or cluster_cols here if using ...
                 annotation_col = annotation_col,
                 annotation_colors = mycolors[c("group_no_sex")],
                 show_rownames = TRUE,
                 show_colnames = FALSE,
                 legend = TRUE,
                 legend_breaks = c(-2, -1, 0, 1, 2),
                 breaks = seq(-2, 2, by = 0.1),
                 color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(41),
                 border_color = NA,
                 treeheight_row = 0,
                 treeheight_col = 10,
                 annotation_names_col = FALSE,
                 annotation_names_row = FALSE,
                 ...)
  return(p1)
}


#pageRank score Normalization
Data <- Data[,rownames(group_sorted)]
Data_normed <- zscore(Data) #zscore
Data_normed_2 <- scaleData(Data) #min_max
tmp <- apply(Data_normed,1:2, function(x) min(x, 2))
Data_normed_3 <- as.data.frame(apply(tmp,1:2, function(x) max(x, -2))) #-2 to 2
#identify non-specific TFs
tmp <- mostConserveTFs(highlyRankedTFs(Data,0.1),0.5) 
housekeeping <- rownames(tmp)
writeLines(housekeeping, paste0('housekeeping_sample',ncol(Data),'top01_CV05.txt'))
# Adjusting the heatmap size
p <- get_heatmap(Data_normed_3[housekeeping,], cluster_rows = FALSE, size = 8) # Increased size from 5 to 8
# Saving the PDF with larger dimensions
pdf(paste0('hp_sample',ncol(Data), '_',length(housekeeping), '_non_specific_TFs.pdf'), width = 15, height = 15)
print(p)
dev.off()

#gsea(x = "housekeeping_top01_CV05_37_new.txt", show.item = 20, is.plot = TRUE)

#Vent mouse morphine RNA seq 
# rna data 48 Sample
rna_seq <- read.csv("/Users/yanmiaodu/Downloads/MOR_SAL/hypo/combined_by_gene_filtered.tsv", header = TRUE, row.names = 1, sep = "\t")
dim(rna_seq)
colnames(rna_seq) <- rownames(meta)[match(gsub("X", "", colnames(rna_seq)), meta$id)]
# selected sample
rna_seq <- rna_seq[,rownames(group_sorted)]

rna_normed_2 <- as.data.frame(t(apply(as.matrix(rna_seq), MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))))
rna_normed_2 <- rna_normed_2[complete.cases(rna_normed_2), ]
zero_counts <- rowSums(rna_normed_2 == 0)
rna_normed_2_filtered <- rna_normed_2[zero_counts <= 1, ]

meta <- meta[, !colnames(meta) %in% c("id")]

annotation_col = data.frame(
   group_sorted$group_no_sex)

p2 <- pheatmap(rna_normed_2_filtered, fontsize = 5, show_rownames = F,
               angle_col = 90, show_colnames = T,
               cellwidth = 7, cellheight = 0.1,
               clustering_distance_cols = 'correlation', 
               clustering_distance_rows = 'correlation', 
               clustering_method = 'average',
               treeheight_row = 0, treeheight_col = 20,
               annotation_col = meta, cutree_cols = 1,
               annotation_colors = mycolors, border_color = NA,
               filename = paste0("hp_",ncol(rna_normed_2_filtered),"_samples_",nrow(rna_normed_2_filtered),"_RNA_2.pdf"),
               width = 12, height = 20)

ordered_RNA <- rna_normed_2_filtered[, order(colnames(rna_normed_2_filtered))]

p1.2 <- pheatmap(ordered_RNA, fontsize = 5, show_rownames = F,
               angle_col = 90, show_colnames = T,
               cellwidth = 8, cellheight = 0.4,
               cluster_cols = F,
               clustering_distance_cols = 'correlation',
               clustering_distance_rows = 'correlation',
               clustering_method = 'average',
               treeheight_row = 0, treeheight_col = 0,
               annotation_col = meta,
               annotation_colors = mycolors, border_color = NA,
               filename = paste0("RNA_all_",ncol(ordered_RNA),"order_samples_",nrow(ordered_RNA),".pdf"))


TFnames <- intersect(rownames(rna_normed_2),rownames(Data))
rna_normed_2_TF <- rna_normed_2[TFnames,]


rna_normed_2 <- as.data.frame(t(apply(as.matrix(rna_normed_2_TF), MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))))

p2 <- pheatmap(rna_normed_2, fontsize = 5, show_rownames = F,
               angle_col = 90, show_colnames = T,
               cellwidth = 7, cellheight = 0.23,
               clustering_distance_cols = 'correlation', 
               clustering_distance_rows = 'correlation', 
               clustering_method = 'average',
               treeheight_row = 0, treeheight_col = 20,
               annotation_col = meta, cutree_cols = 1,
               annotation_colors = mycolors, border_color = NA,
               filename = paste0("hp_rna_",ncol(rna_normed_2),"_samples_",nrow(rna_normed_2),"_TFs_2.pdf"),
               width = 12, height = 12)

RNA_selected <- mostVariableTFs(rna_normed_2_filtered,cv.low = 1,cv.high = 2) 
dim(RNA_selected)
p1.3 <- pheatmap(RNA_selected, fontsize = 5, show_rownames = F,
               angle_col = 90, show_colnames = T,
               cellwidth = 7, cellheight = 0.1,
               clustering_distance_cols = 'correlation', 
               clustering_distance_rows = 'correlation', 
               clustering_method = 'average',
               treeheight_row = 0, treeheight_col = 20,
               annotation_col = meta, cutree_cols = 1,
               annotation_colors = mycolors, border_color = NA,
               filename = paste0("hp_",ncol(RNA_selected),"_samples_",nrow(RNA_selected),"_RNA_cv_1_2.pdf"),
               width = 12, height = 20)

TFnames <- intersect(rownames(rna_seq),rownames(Data))
rna_seq_TF <- rna_seq[TFnames,]
rna_seq_TF <- rna_seq_TF[rowSums(rna_seq_TF)>0,]
Data_TF <- Data[TFnames,]

#no need to split samples to different part
tfs <- lapply(unique(group_sorted$group_no_sex),
            function(x) {
             pr <- Data_TF[!(row.names(rna_seq_TF) %in% housekeeping),]
             n <- which(group_sorted[["group_no_sex"]]==x)
             x <- unpairedtTest(pr,n,p=0.05,d=0.5)
             t <- rownames(x)
             t <- t[unlist(lapply(t, function(x) sum(rna_seq[x, n]>0) >= ceiling(0.2*length(n))))]
            }) 
names(tfs) <- unique(group_sorted$group_no_sex)
for (i in seq_along(tfs)) {
  writeLines(tfs[[i]], paste0(names(tfs)[i], "_CS_filter.txt"))
}

ps <- lapply(names(tfs), function(x) {
  annotation_row <- data.frame("subsets"=rep(x,length(tfs[[x]])))
  rownames(annotation_row) <- tfs[[x]]
  p <- get_heatmap(Data_normed_3[rownames(annotation_row),], 
                 cluster_rows=F, 
                 annotation_row=annotation_row, 
                 size = 4)  # Removed cellwidth and cellheight here
  pdf(paste0('hp_',nrow(annotation_row),'_',x,'_TFs.pdf'))
  print(p)
  dev.off()
})

#identify multi-states specific TFs and single-state specific TFs
getRowAnnotation <- function(TFs){
  tmp <- unlist(TFs)
  redundants <- unname(tmp[duplicated(tmp)])
  uniques <- unname(tmp[!duplicated(tmp)])
  tmp2 <- tibble("subsets" = "multi","TFs"=unique(redundants))
  annotation_row <- tibble("subsets" = unlist(lapply(1:length(TFs), function(x) {rep(names(TFs[x]), length(TFs[[x]]))})),
                           "TFs" = unlist(TFs)) %>% 
    filter(TFs %in% setdiff(uniques, redundants)) %>%
    bind_rows(tmp2) %>%
    tibble::column_to_rownames("TFs")
  return(annotation_row)
  
}

annotation_row <- getRowAnnotation(tfs)
annotation_row %>% filter(subsets=="multi") %>% rownames() %>% writeLines('multi_TFs_group.txt')
annotation_row %>% filter(subsets!="multi") %>% rownames() %>% writeLines('single_TFs_group.txt')

df_sgl <- Data_normed_3[readLines('single_TFs_group.txt'),]

# add gap between cell state groups
annotation_row_single <- annotation_row |> filter(subsets!='multi')
last_element_per_group <- aggregate(rownames(annotation_row_single), by=list(annotation_row_single$subsets), FUN=tail, n=1) |> pull(x)
gaps <- sort(unlist(lapply(last_element_per_group, function(x) which(x==rownames(df_sgl)))))

p1 <- get_heatmap(df_sgl, size=4, gaps_row = gaps, cluster_rows=F, annotation_row = annotation_row_single)
pdf(paste0('hp_',nrow(df_sgl),'_morVSsal_TFs.pdf'), width = 8, height=25)
print(p1)
dev.off()



#in same treatment test sal vs lps
LPS_cols <- grepl("LPS", colnames(Data))  # Matches columns contain LPS"
sal_cols <- grepl("SAL", colnames(Data)) & !grepl("LPS", colnames(Data))
DATA_LPS <- Data[, LPS_cols]
DATA_sal <- Data[, sal_cols]

#for test specific effect base on MOR on LPS response , U wanan fine CS tf in MOR.LPS vs SAL.LPS
group_sorted_LPS <- group_sorted[grepl("LPS", group_sorted$treatment), ]
# for Treatment LPS 
tfs <- lapply(unique(group_sorted_LPS$condition), function(x) {
  # Filter out housekeeping genes from the Data matrix
  pr <- DATA_LPS[!(row.names(DATA_LPS) %in% housekeeping),]
  # Identify columns corresponding to the group
  n <- which(group_sorted_LPS[["condition"]] == x)
  # Apply unpaired t-test (adjust the p-value and effect size as needed)
  significant_genes <- unpairedtTest(pr, n, p = 0.05, d = 0.5)
  # Get the rownames (gene names) of the significant genes
  sig_gene_names <- rownames(significant_genes)
  # Further filter for genes expressed in at least 20% of the cells in that group
  filtered_genes <- sig_gene_names[unlist(lapply(sig_gene_names, function(gene) {
    sum(pr[gene, n] > 0) >= ceiling(0.2 * length(n))
  }))]
  return(filtered_genes)
})
names(tfs) <- unique(group_sorted_LPS$condition)
for (i in seq_along(tfs)) {
  writeLines(tfs[[i]], paste0(names(tfs)[i], "_LPS_tre_filter_P001.txt"))
}

annotation_row <- getRowAnnotation(tfs)
DF = Data_normed_3[rownames(annotation_row),rownames(group_sorted_LPS)]
ps <- lapply(names(tfs), function(x) {
  p <- get_clean_heatmap(DF, 
                 cluster_rows=F, 
                 annotation_row=annotation_row, 
                 size = 4)  # Removed cellwidth and cellheight here
  pdf(paste0('hp_',nrow(annotation_row),'_',x,'cleaner_LPS_TFs.pdf'))
  print(p)
  dev.off()
})

group_sorted_SAL <- group_sorted[grepl("SAL", group_sorted$treatment), ]
tfs <- lapply(unique(group_sorted_SAL$condition), function(x) {
  # Filter out housekeeping genes from the Data matrix
  pr <- DATA_sal[!(row.names(DATA_sal) %in% housekeeping),]
  # Identify columns corresponding to the group
  n <- which(group_sorted_SAL[["condition"]] == x)
  # Apply unpaired t-test (adjust the p-value and effect size as needed)
  significant_genes <- unpairedtTest(pr, n, p = 0.05, d = 0.5)
  # Get the rownames (gene names) of the significant genes
  sig_gene_names <- rownames(significant_genes)
  # Further filter for genes expressed in at least 20% of the cells in that group
  filtered_genes <- sig_gene_names[unlist(lapply(sig_gene_names, function(gene) {
    sum(pr[gene, n] > 0) >= ceiling(0.2 * length(n))
  }))]
  return(filtered_genes)
})
names(tfs) <- unique(group_sorted_SAL$condition)
for (i in seq_along(tfs)) {
  writeLines(tfs[[i]], paste0(names(tfs)[i], "_sal_CS_filter.txt"))
}
annotation_row <- getRowAnnotation(tfs)
DF = Data_normed_3[rownames(annotation_row),rownames(group_sorted_SAL)]
ps <- lapply(names(tfs), function(x) {
  p <- get_clean_heatmap(DF, 
                 cluster_rows=F, 
                 annotation_row=annotation_row, 
                 size = 4)  # Removed cellwidth and cellheight here
  pdf(paste0('hp_',nrow(annotation_row),'_',x,'_TFs.pdf'))
  print(p)
  dev.off()
})

# 3. GO term analysis
# L <- list.files(path = "./", pattern = "_filter.txt")
# lapply(L, function(x) gsea(x = x, show.item = 20))

files <- list.files(pattern = "^[^~$].*_GO_KEGG\\.xlsx$")

# Function to create a dot plot for each Excel file
for (file in files) {
  # Extract basename (without extension) to use for naming plot files
  basename <- sub("\\.xlsx$", "", file)
  
  # Read the data from the first sheet of the Excel file
  data <- read_excel(file, sheet = 1)
  
  # Filter for rows where the 'ID' column contains "GO:"
  if ("ID" %in% colnames(data)) {
    go_data <- data[grep("GO:", data$ID), ]
  } else {
    # Skip the file if 'ID' column is not present
    print(paste("Skipping file", basename, "- 'ID' column not found"))
    next
  }
  
  # Check if the necessary columns are present
  if (!all(c("p.adjust", "GeneRatio", "Description") %in% colnames(go_data))) {
    # Print a message if required columns are missing
    print(paste("Skipping file", basename, "- required columns are missing"))
    next
  }
  
  # Select the top 20 GO terms based on p.adjust (q-value)
  top_go_data <- go_data[order(go_data$p.adjust), ][1:20, ]
  
  # Create the dot plot
  dot_plot <- ggplot(top_go_data, aes(x = p.adjust, y = reorder(Description, p.adjust), size = GeneRatio)) +
    geom_point() +
    theme_bw() +
    labs(
      title = paste("Top 20 GO Term Enrichment for", basename),
      x = "p.adjust",
      y = "GO Terms",
      size = "GeneRatio"
    ) +
    theme(axis.text.y = element_text(size = 8)) +
    theme(plot.title = element_text(hjust = 0.5))

  # Save the plot as PNG file
  png_filename <- paste0(basename, "_GO_dotplot.png")
  ggsave(filename = png_filename, plot = dot_plot, width = 8, height = 6)
}


for (file in files) {
  # Extract basename (without extension) to use for naming plot files
  basename <- sub("\\.xlsx$", "", file)
  
  # Read the data from the first sheet of the Excel file
  data <- read_excel(file, sheet = 1)
  
  # Filter for rows where the 'ID' column contains "GO:"
  if ("ID" %in% colnames(data)) {
    go_data <- data[grep("GO:", data$ID), ]
  } else {
    # Skip the file if 'ID' column is not present
    print(paste("Skipping file", basename, "- 'ID' column not found"))
    next
  }
  
  # Check if the necessary columns are present
  if (!all(c("p.adjust", "GeneRatio", "Description") %in% colnames(go_data))) {
    # Print a message if required columns are missing
    print(paste("Skipping file", basename, "- required columns are missing"))
    next
  }
  
  # Select the top 20 GO terms based on p.adjust (q-value)
  top_go_data <- go_data[order(go_data$p.adjust), ][1:20, ]
  
  # Create the dot plot
  dot_plot <- ggplot(top_go_data, aes(x = p.adjust, y = reorder(Description, -p.adjust), size = GeneRatio)) +
  geom_point() +
  theme_bw() +
  labs(
    title = paste("Top 20 GO Term Enrichment for", basename),
    x = "p.adjust",
    y = "GO Terms",
    size = "GeneRatio"
  ) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(plot.title = element_text(hjust = 0.5))

  # Save the plot as PNG file
  png_filename <- paste0(basename, "_GO_dotplot.png")
  ggsave(filename = png_filename, plot = dot_plot, width = 8, height = 6)
}

#maybe also some general when MOR vs LPS regarding lps or sal or F/M

# select most variable TFs
Data_selected <- mostVariableTFs(Data,cv.low = 1.0,cv.high = 2) 
Data_selected_normed <- as.data.frame(t(apply(as.matrix(Data_selected), MARGIN = 1, FUN = function(X) (X - min(X))/diff(range(X)))))
annotation_col = group_sorted
p3 <- pheatmap(Data_selected_normed, fontsize = 5, show_rownames = T,
               angle_col = 90, show_colnames = T,
               cellwidth = 5, cellheight = 5,
               cluster_cols = F,
               clustering_distance_cols = 'correlation', 
               clustering_distance_rows = 'correlation', 
               clustering_method = 'average',
               treeheight_row = 0, treeheight_col = 0,
               annotation_col = annotation_col,
               annotation_colors = mycolors, border_color = NA,
               filename = paste0("hp_mv_",ncol(Data_selected_normed),"_08_15_cv_samples_",nrow(Data_selected_normed),"_TFs.pdf"))
high_cv_tf = as.data.frame(rownames(Data_selected))
write.table(high_cv_tf,"highcv_TF_pagerank.txt",quote = F,row.names = F,col.names = F)

#cluster by variable genens
L <- list.files(path = "./", pattern = "*_filter.txt")
gene_list <- unlist(lapply(L, readLines))
# Get unique gene names
unique_genes <- unique(gene_list)
# Subset the pagerank scores data frame to only include the key TFs

pagerank_scores_key_tfs <- Data_normed_2[rownames(Data_normed_2) %in% unique_genes, , drop = FALSE]
# Plot using pheatmap with clustering enabled for columns and rows using the subset of key TFs
p2 <- pheatmap(pagerank_scores_key_tfs, fontsize = 5, show_rownames = T,
               angle_col = 90, show_colnames = T,
               cellwidth = 8, cellheight = 4,
               clustering_cols = TRUE, # Enable clustering for columns
               clustering_rows = TRUE, # Enable clustering for rows
               clustering_distance_cols = "correlation", 
               clustering_distance_rows = "correlation",
               clustering_method = 'average',
               treeheight_row = 50, treeheight_col = 20, # Set height for dendrograms
               annotation_col = annotation_col, 
               annotation_colors = mycolors, border_color = NA,
               filename = "pagerank_scores_new_key_tfs_clustering_average_30.pdf",
               width = 12, height = 12)


#Order the annotation columns by group information if necessary
ordered_data <- Data_normed_2[, order(colnames(Data_normed_2))]
ordered_cols <- annotation_col[order(annotation_col$group), , drop = FALSE]
ordered_sample_names <- rownames(ordered_cols)

# Reorder the columns in the pagerank scores data frame based on the ordered sample names
pagerank_scores_ordered <- ordered_data[, ordered_sample_names, drop = FALSE]
colnames(pagerank_scores_ordered) <- ordered_sample_names

p1.2 <- pheatmap(pagerank_scores_ordered, fontsize = 5, show_rownames = F,
               angle_col = 90, show_colnames = T,
               cellwidth = 8, cellheight = 0.4,
               cluster_cols = F,
               clustering_distance_cols = 'correlation',
               clustering_distance_rows = 'correlation',
               clustering_method = 'average',
               treeheight_row = 0, treeheight_col = 0,
               annotation_col = annotation_col,
               annotation_colors = mycolors, border_color = NA,
               filename = paste0("hp_all_",ncol(Data_normed),"order_samples_",nrow(Data_normed),"_TFs_2.pdf"))

p1 <- pheatmap(Data_normed_2, fontsize = 5, show_rownames = F,
               angle_col = 90, show_colnames = T,
               cellwidth = 6, cellheight = 0.2,
               clustering_distance_cols = 'correlation', 
               clustering_distance_rows = 'correlation', 
               clustering_method = 'average',
               treeheight_row = 0, treeheight_col = 25,
               annotation_col = annotation_col, cutree_cols = 1,
               annotation_colors = mycolors, border_color = NA,
               filename = paste0("hp_LPS_",ncol(Data_normed_2),"_samples_",nrow(Data_normed_2),"_TFs_2.pdf"),
               width = 8, height = 12)



group_sorted_LPS = group_sorted[group_sorted$treatment == "LPS",]
meta_sorted <- meta[rownames(group_sorted), ]
meta_LPS = meta_sorted[meta_sorted$treatment == "LPS",]
group_sorted_LPS$sex = meta_LPS$sex
group_sorted_LPS$group_no_sex = NULL
Data_normed_2_LPS = Data_normed_2[,rownames(group_sorted_LPS)]
p1 <- pheatmap(Data_normed_2_LPS, fontsize = 5, show_rownames = F,
               angle_col = 90, show_colnames = T,
               cellwidth = 6, cellheight = 0.2,
               clustering_distance_cols = 'correlation', 
               clustering_distance_rows = 'correlation', 
               clustering_method = 'average',
               treeheight_row = 0, treeheight_col = 25,
               annotation_col = group_sorted_LPS, cutree_cols = 1,
               annotation_colors = mycolors, border_color = NA,
               filename = paste0("hp_LPS_",ncol(Data_normed_2_LPS),"_samples_",nrow(Data_normed_2_LPS),"_TFs_2.pdf"),
               width = 8, height = 12)


group_sorted_SAL = group_sorted[group_sorted$treatment == "SAL",]
meta_sorted <- meta[rownames(group_sorted), ]
meta_SAL = meta_sorted[meta_sorted$treatment == "SAL",]
group_sorted_SAL$sex = meta_SAL$sex
group_sorted_SAL$group_no_sex = NULL
Data_normed_2_SAL = Data_normed_2[,rownames(group_sorted_SAL)]

p1 <- pheatmap(Data_normed_2_SAL, fontsize = 5, show_rownames = F,
               angle_col = 90, show_colnames = T,
               cellwidth = 6, cellheight = 0.2,
               clustering_distance_cols = 'correlation', 
               clustering_distance_rows = 'correlation', 
               clustering_method = 'average',
               treeheight_row = 0, treeheight_col = 25,
               annotation_col = group_sorted_SAL, cutree_cols = 1,
               annotation_colors = mycolors, border_color = NA,
               filename = paste0("hp_LPS_",ncol(Data_normed_2_SAL),"_samples_",nrow(Data_normed_2_SAL),"_TFs_2.pdf"),
               width = 8, height = 12)


