#=============================================================
# Load Required Packages 
#=============================================================
library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(cowplot)
library(patchwork)
library(clustree)
library(pheatmap)
library(RColorBrewer)

#=============================================================
# Set Working and Data Directories
#=============================================================
# Use forward slashes for cross-platform compatibility
work_dir <- "C:/Users/xieru/Desktop/GEO10X"
setwd(work_dir)

data_dir <- "C:/Users/xieru/Desktop/GEO10X/RAW"

if (!dir.exists(data_dir)) {
  stop("Data directory does not exist: ", data_dir, "\nPlease create it and place 10X data inside.")
}

cat("Working Directory:", work_dir, "\n")
cat("Data Directory:", data_dir, "\n\n")

#=============================================================
# Function: Process Single Sample 
#=============================================================
process_single_sample <- function(folder_path, sample_name, condition) {
  counts <- Read10X(data.dir = folder_path)
  
  seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  
  seurat_obj$batch <- sample_name
  seurat_obj$condition <- condition
  
  # Calculate QC metrics
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^HB[^(P)]")
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  
  return(seurat_obj)
}

#=============================================================
# Function: Determine Sample Condition
#=============================================================
get_condition <- function(sample_name) {
  if (grepl("Normal|Control|Ctrl", sample_name, ignore.case = TRUE)) {
    return("Control")
  } else if (grepl("TLE|Epilepsy|Patient|Disease", sample_name, ignore.case = TRUE)) {
    return("TLE")
  } else {
    return("Unknown")
  }
}

#=============================================================
# Batch Read All Samples
#=============================================================
all_folders <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)

cat("========== Scanning Data Directory ==========\n")
cat("Found", length(all_folders), "folders\n\n")

valid_folders <- c()
for (folder in all_folders) {
  folder_full_path <- file.path(data_dir, folder)
  
  has_matrix <- file.exists(file.path(folder_full_path, "matrix.mtx.gz")) || 
                file.exists(file.path(folder_full_path, "matrix.mtx"))
  has_barcodes <- file.exists(file.path(folder_full_path, "barcodes.tsv.gz")) || 
                  file.exists(file.path(folder_full_path, "barcodes.tsv"))
  has_features <- file.exists(file.path(folder_full_path, "features.tsv.gz")) || 
                  file.exists(file.path(folder_full_path, "genes.tsv.gz")) ||
                  file.exists(file.path(folder_full_path, "features.tsv"))
  
  if (has_matrix && has_barcodes && has_features) {
    valid_folders <- c(valid_folders, folder)
    cat("✓ Valid sample:", folder, "\n")
  } else {
    cat("✗ Skipped:", folder, "(Missing required files)\n")
  }
}

cat("\n========== Detected", length(valid_folders), "valid 10X sample folders ==========\n\n")

if (length(valid_folders) == 0) {
  stop("No valid 10X data folders found! Please check the data directory and file formats.")
}

#=============================================================
# Process All Valid Samples
#=============================================================
sample_list <- list()

for (folder in valid_folders) {
  sample_name <- folder
  condition <- get_condition(sample_name)
  
  cat("Processing:", sample_name, "| Condition:", condition, "\n")
  
  folder_full_path <- file.path(data_dir, folder)
  
  sample_list[[sample_name]] <- process_single_sample(
    folder_path = folder_full_path,
    sample_name = sample_name,
    condition = condition
  )
  
  cat("  -> Cells:", ncol(sample_list[[sample_name]]), "\n")
  cat("  -> Genes:", nrow(sample_list[[sample_name]]), "\n\n")
}

#=============================================================
# Sample Information Summary
#=============================================================
cat("========== Sample Information Summary ==========\n")
sample_summary <- data.frame(
  Sample = names(sample_list),
  Condition = sapply(sample_list, function(x) unique(x$condition)),
  Cells = sapply(sample_list, ncol),
  Genes = sapply(sample_list, nrow)
)
print(sample_summary)

unknown_samples <- sample_summary$Sample[sample_summary$Condition == "Unknown"]
if (length(unknown_samples) > 0) {
  cat("\n⚠️ WARNING: The following samples could not be automatically grouped. Please assign manually:\n")
  print(unknown_samples)
}

#=============================================================
# Step 1: Pre-QC Visualization
#=============================================================
cat("\n========== Generating Pre-QC Plots ==========\n")

merged_pre_qc <- merge(sample_list[[1]], 
                       y = sample_list[2:length(sample_list)], 
                       add.cell.ids = names(sample_list))

pre_qc_plots <- VlnPlot(merged_pre_qc, 
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"),
                        group.by = "batch", ncol = 4, pt.size = 0)

pdf("01_pre_QC_violin.pdf", width = 20, height = 6)
print(pre_qc_plots)
dev.off()

scatter1 <- FeatureScatter(merged_pre_qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "batch", pt.size = 0.5) + NoLegend()
scatter2 <- FeatureScatter(merged_pre_qc, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "batch", pt.size = 0.5) + NoLegend()
scatter3 <- FeatureScatter(merged_pre_qc, feature1 = "nCount_RNA", feature2 = "percent.hb", group.by = "batch", pt.size = 0.5) + NoLegend()

pdf("01_pre_QC_scatter.pdf", width = 18, height = 5)
print(scatter1 + scatter2 + scatter3)
dev.off()

cat("\nPre-QC Statistics per Sample:\n")
pre_qc_stats <- data.frame(
  Sample = names(sample_list),
  Cells = sapply(sample_list, ncol),
  Median_nFeature = sapply(sample_list, function(x) median(x$nFeature_RNA)),
  Median_nCount = sapply(sample_list, function(x) median(x$nCount_RNA)),
  Median_MT = sapply(sample_list, function(x) round(median(x$percent.mt), 2)),
  Median_HB = sapply(sample_list, function(x) round(median(x$percent.hb), 2))
)
print(pre_qc_stats)

rm(merged_pre_qc)
gc()

#=============================================================
# Step 2: Quality Control Filtering
#=============================================================
cat("\n========== Quality Control Filtering ==========\n")

qc_before <- data.frame(
  Sample = names(sample_list),
  Cells_Before = sapply(sample_list, ncol)
)

for (sample_name in names(sample_list)) {
  cells_before <- ncol(sample_list[[sample_name]])
  
  sample_list[[sample_name]] <- subset(
    sample_list[[sample_name]], 
    subset = nFeature_RNA > 200 & 
             nFeature_RNA < 6000 & 
             percent.mt < 20 &
             percent.hb < 1
  )
  
  cells_after <- ncol(sample_list[[sample_name]])
  cat("Sample", sample_name, ":", cells_before, "->", cells_after, 
      "cells (Removed", cells_before - cells_after, ")\n")
}

qc_after <- data.frame(
  Sample = names(sample_list),
  Cells_Before = qc_before$Cells_Before,
  Cells_After_QC = sapply(sample_list, ncol)
)
qc_after$Cells_Removed <- qc_after$Cells_Before - qc_after$Cells_After_QC
qc_after$Percent_Kept <- round(qc_after$Cells_After_QC / qc_after$Cells_Before * 100, 1)

cat("\n========== QC Filtering Summary ==========\n")
print(qc_after)
write.csv(qc_after, "QC_summary_no_doublet.csv", row.names = FALSE)

#=============================================================
# Step 3: Post-QC Visualization
#=============================================================
merged_seurat <- merge(sample_list[[1]], 
                       y = sample_list[2:length(sample_list)], 
                       add.cell.ids = names(sample_list))

post_qc_plots <- VlnPlot(merged_seurat, 
                         features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.hb"),
                         group.by = "batch", ncol = 4, pt.size = 0)

pdf("02_post_QC_violin.pdf", width = 20, height = 6)
print(post_qc_plots)
dev.off()

scatter1_post <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "batch", pt.size = 0.5) + NoLegend()
scatter2_post <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = "batch", pt.size = 0.5) + NoLegend()
scatter3_post <- FeatureScatter(merged_seurat, feature1 = "nCount_RNA", feature2 = "percent.hb", group.by = "batch", pt.size = 0.5) + NoLegend()

pdf("02_post_QC_scatter.pdf", width = 18, height = 5)
print(scatter1_post + scatter2_post + scatter3_post)
dev.off()

cat("\nTotal cells after merging:", ncol(merged_seurat), "\n")

rm(sample_list)
gc()

#=============================================================
# Step 4: Normalization & Harmony Integration
#=============================================================
cat("\n========== Starting Normalization and Integration ==========\n")

merged_seurat <- merged_seurat %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData(vars.to.regress = c("percent.mt", "percent.hb"))

merged_seurat <- RunPCA(merged_seurat, npcs = 30)

merged_seurat <- IntegrateLayers(
  object = merged_seurat,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)

saveRDS(merged_seurat, file = "merged_seurat_QC_filtered_no_doublet.rds")

pdf("03_pca_elbow_plot.pdf", width = 8, height = 6)
print(ElbowPlot(merged_seurat, reduction = "pca", ndims = 30) +
  ggtitle("PCA Elbow Plot") +
  theme_classic())
dev.off()

cat("\n========== QC and Integration Complete! ==========\n")

#=============================================================
# Step 5: Clustering Analysis (Multiple Resolutions)
#=============================================================
cat("\n========== Starting Clustering Analysis ==========\n")

merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:20)

resolutions <- seq(0.1, 1.0, by = 0.1)

for (res in resolutions) {
  cat("Calculating resolution =", res, "\n")
  merged_seurat <- FindClusters(merged_seurat, resolution = res, algorithm = 1, verbose = FALSE)
}

merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:20)

#=============================================================
# Step 6: Clustree Analysis
#=============================================================
pdf("04_clustree_analysis.pdf", width = 12, height = 16)
print(clustree(merged_seurat, prefix = "RNA_snn_res.") +
  theme(legend.position = "right") +
  ggtitle("Clustree Analysis - Resolution Selection"))
dev.off()

cat("Clustree plot saved: 04_clustree_analysis.pdf\n")

#=============================================================
# Step 7: UMAP Visualization at Different Resolutions
#=============================================================
cat("\n========== Plotting UMAPs for All Resolutions ==========\n")

plot_list <- list()
for (res in resolutions) {
  p <- DimPlot(merged_seurat, 
               reduction = "umap", 
               group.by = paste0("RNA_snn_res.", res),
               label = TRUE, label.size = 4, repel = TRUE) +
    ggtitle(paste("Resolution:", res)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  pdf(paste0("05_UMAP_res", res, ".pdf"), width = 10, height = 8)
  print(p)
  dev.off()
  
  # Store for combined plot
  p_combined <- p + NoLegend() + theme(plot.title = element_text(size = 10))
  plot_list[[as.character(res)]] <- p_combined
}

pdf("05_UMAP_all_resolutions.pdf", width = 25, height = 20)
print(wrap_plots(plot_list, ncol = 4))
dev.off()

cat("UMAP plots saved.\n")

#=============================================================
# Step 8: UMAP by Sample and Condition
#=============================================================
p_batch <- DimPlot(merged_seurat, reduction = "umap", group.by = "batch", label = FALSE) +
  ggtitle("UMAP by Sample")

p_condition <- DimPlot(merged_seurat, reduction = "umap", group.by = "condition",
                       cols = c("Control" = "#4DAF4A", "TLE" = "#E41A1C")) +
  ggtitle("UMAP by Condition")

p_split <- DimPlot(merged_seurat, reduction = "umap", split.by = "batch", group.by = "condition",
                   cols = c("Control" = "#4DAF4A", "TLE" = "#E41A1C"), ncol = 4) +
  ggtitle("UMAP Split by Sample")

pdf("06_UMAP_by_sample_condition.pdf", width = 16, height = 6)
print(p_batch + p_condition)
dev.off()

pdf("06_UMAP_split_by_sample.pdf", width = 20, height = 10)
print(p_split)
dev.off()

#=============================================================
# Step 9: Cluster Count Statistics
#=============================================================
cat("\n========== Cluster Counts per Resolution ==========\n")

cluster_stats <- data.frame(
  Resolution = resolutions,
  Num_Clusters = sapply(resolutions, function(res) {
    length(unique(merged_seurat[[paste0("RNA_snn_res.", res)]][,1]))
  })
)
print(cluster_stats)

write.csv(cluster_stats, "cluster_stats_by_resolution.csv", row.names = FALSE)

#=============================================================
# Step 10: Save Pre-Annotation Results
#=============================================================
saveRDS(merged_seurat, file = "merged_seurat_clustered_no_doublet.rds")

cat("\n========== Clustering Analysis Complete! ==========\n")
cat("Please review '04_clustree_analysis.pdf' to select the optimal resolution before proceeding with cell annotation.\n")

#=============================================================
# ↓↓↓ CELL ANNOTATION PIPELINE ↓↓↓
#=============================================================

# Join layers first (Required for Seurat v5)
merged_seurat <- JoinLayers(merged_seurat)

# Set resolution (Modify this based on Clustree results)
chosen_res <- 0.1  

merged_seurat <- FindClusters(merged_seurat, resolution = chosen_res)
Idents(merged_seurat) <- paste0("RNA_snn_res.", chosen_res)

cat("\nCurrent active resolution:", chosen_res, "\n")
cat("Number of clusters:", length(unique(Idents(merged_seurat))), "\n")
cat("Cells per cluster:\n")
print(table(Idents(merged_seurat)))

#=============================================================
# Step 11: Find Marker Genes
#=============================================================
cat("\n========== Identifying Marker Genes ==========\n")

markers <- FindAllMarkers(
  merged_seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

# Filter for significant markers
markers_sorted <- markers %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.csv(markers, paste0("all_markers_res", chosen_res, ".csv"), row.names = FALSE)
write.csv(markers_sorted, paste0("top10_markers_res", chosen_res, ".csv"), row.names = FALSE)

cat("Marker genes saved.\n")

#=============================================================
# Step 12: Marker Gene Visualization
#=============================================================
# Heatmap - Top 5 genes per cluster
top5_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

pdf(paste0("07_markers_heatmap_res", chosen_res, ".pdf"), width = 14, height = 10)
print(DoHeatmap(merged_seurat, features = top5_markers,
          group.colors = scales::hue_pal()(length(unique(Idents(merged_seurat))))) +
  theme(axis.text.y = element_text(size = 7)) +
  ggtitle(paste("Top 5 Markers per Cluster (Resolution:", chosen_res, ")")))
dev.off()

# Dotplot - Top 3 genes per cluster
top3_markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC) %>%
  pull(gene) %>%
  unique()

pdf(paste0("07_markers_dotplot_res", chosen_res, ".pdf"), width = 16, height = 8)
print(DotPlot(merged_seurat, features = top3_markers, cols = c("lightgrey", "red")) +
  RotatedAxis() +
  ggtitle(paste("Top 3 Markers per Cluster (Resolution:", chosen_res, ")")))
dev.off()

#=============================================================
# Step 13: Automated Cell Type Annotation
#=============================================================
Idents(merged_seurat) <- "RNA_snn_res.0.1"

print("Current active clustering:")
print("Resolution: RNA_snn_res.0.1")
print(paste("Number of clusters:", length(levels(Idents(merged_seurat)))))

# Define cell type markers (Merged duplicated Endothelial_Cell)
cell_type_markers <- list(
  "Choroid_Plexus_Epithelial_Cell"= c("LRP2", "ST18", "FOLH1"),
  "Neurons"                       = c("NEFL"),
  "Endothelial_Cell"              = c("SLC22A25", "SLC22A24", "CLDN5", "MECOM", "EBF1"),
  "Fibroblast"                    = c("PRDM16", "COL1A1", "DCN", "LUM", "GLI2", "GLI3"), 
  "Microglia"                     = c("CD86", "SYK"), 
  "Interneuron"                   = c("RELN", "GAD1", "SST", "PVALB"), 
  "Oligodendrocyte_precursor_cells"= c("PDGFRA", "CSPG4", "SOX10")
)

# Check which markers exist in the dataset
all_markers <- unique(unlist(cell_type_markers))
present_markers <- all_markers[all_markers %in% rownames(merged_seurat)]
missing_markers <- all_markers[!all_markers %in% rownames(merged_seurat)]

cat("\n=== Marker Gene Check ===\n")
cat("Total markers provided:", length(all_markers), "\n")
cat("Markers present in data:", length(present_markers), "\n")
cat("Missing markers:", length(missing_markers), "\n")

if(length(missing_markers) > 0) {
  print("Missing markers:")
  print(missing_markers)
}

# Filter out missing genes from the marker list
cell_type_markers_filtered <- lapply(cell_type_markers, function(genes) {
  genes[genes %in% rownames(merged_seurat)]
})

# Remove empty cell types
cell_type_markers_filtered <- cell_type_markers_filtered[sapply(cell_type_markers_filtered, length) > 0]

# Calculate average expression of present markers
cat("\n=== Calculating Average Expression ===\n")
cluster_avg_exp <- AverageExpression(merged_seurat, 
                                     features = present_markers,
                                     group.by = "RNA_snn_res.0.1")$RNA

# Calculate cell type scores for each cluster
cluster_scores <- data.frame(cluster = colnames(cluster_avg_exp))
rownames(cluster_scores) <- cluster_scores$cluster

for(cell_type in names(cell_type_markers_filtered)) {
  markers_list <- cell_type_markers_filtered[[cell_type]]
  if(length(markers_list) > 0) {
    if(length(markers_list) == 1) {
      scores <- cluster_avg_exp[markers_list, ]
    } else {
      scores <- colMeans(cluster_avg_exp[markers_list, , drop = FALSE])
    }
    cluster_scores[[cell_type]] <- scores
  }
}

# Assign the highest scoring cell type to each cluster
cluster_annotations <- apply(cluster_scores[, -1, drop = FALSE], 1, function(x) {
  if(all(is.na(x))) return("Unknown")
  return(names(which.max(x)))
})

cat("\n=== Automated Annotation Results ===\n")
print(cluster_annotations)

# Clean cluster names (remove 'g' prefix if present)
cluster_to_celltype <- cluster_annotations
names(cluster_to_celltype) <- gsub("^g", "", names(cluster_to_celltype))

current_clusters <- levels(Idents(merged_seurat))
mapping_clusters <- names(cluster_to_celltype)

missing_in_mapping <- setdiff(current_clusters, mapping_clusters)
extra_in_mapping <- setdiff(mapping_clusters, current_clusters)

if(length(missing_in_mapping) > 0) {
  cat("\nWarning: Clusters missing in mapping, assigning 'Unknown':\n")
  print(missing_in_mapping)
  for(missing_cluster in missing_in_mapping) {
    cluster_to_celltype[missing_cluster] <- "Unknown"
  }
}

if(length(extra_in_mapping) > 0) {
  cluster_to_celltype <- cluster_to_celltype[names(cluster_to_celltype) %in% current_clusters]
}

# Apply annotations
cat("\n=== Applying Cell Type Annotations ===\n")
merged_seurat <- RenameIdents(merged_seurat, cluster_to_celltype)
merged_seurat$cell_type_auto <- Idents(merged_seurat)

cat("\nCell Counts by Type:\n")
print(table(merged_seurat$cell_type_auto))

# Generate Annotation Summary
annotation_summary <- data.frame(
  Cluster = names(cluster_to_celltype),
  Cell_Type = cluster_to_celltype,
  Cell_Count = as.numeric(table(merged_seurat$RNA_snn_res.0.1)[names(cluster_to_celltype)]),
  stringsAsFactors = FALSE
)

cell_type_summary <- aggregate(Cell_Count ~ Cell_Type, data = annotation_summary, sum)
cell_type_summary <- cell_type_summary[order(cell_type_summary$Cell_Count, decreasing = TRUE), ]

#=============================================================
# Step 14: Export Annotation Visualizations
#=============================================================
pdf("umap_auto_annotation.pdf", width = 16, height = 6)
p1 <- DimPlot(merged_seurat, reduction = "umap", group.by = "RNA_snn_res.0.1", label = TRUE, label.size = 3) +
  ggtitle("Original Clusters (Resolution 0.1)") + theme_classic()
p2 <- DimPlot(merged_seurat, reduction = "umap", group.by = "cell_type_auto", label = TRUE, label.size = 3) +
  ggtitle("Auto-annotated Cell Types") + theme_classic()
print(p1 + p2)
dev.off()

pdf("cluster_celltype_scores_heatmap.pdf", width = 12, height = 8)
score_matrix <- as.matrix(cluster_scores[, -1])
rownames(score_matrix) <- gsub("^g", "", rownames(score_matrix))
pheatmap(t(score_matrix), cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
         main = "Cell Type Scores for Each Cluster", fontsize = 8, cellwidth = 15, cellheight = 15)
dev.off()

# Validation Plots
pdf("marker_validation.pdf", width = 15, height = 20)
for(cell_type in names(cell_type_markers_filtered)) {
  markers_list <- cell_type_markers_filtered[[cell_type]]
  if(length(markers_list) > 0) {
    markers_to_plot <- markers_list[1:min(4, length(markers_list))]
    p <- FeaturePlot(merged_seurat, features = markers_to_plot, ncol = 2, pt.size = 0.5) +
         plot_annotation(title = paste("Markers for", gsub("_", " ", cell_type)))
    print(p)
  }
}
dev.off()

# Clean DotPlot
dotplot_pretty <- DotPlot(merged_seurat, features = present_markers, group.by = "cell_type_auto") + 
  coord_flip() + theme_bw() +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = NULL, y = NULL) +
  scale_color_gradientn(values = seq(0, 1, 0.2), colours = c('#330066', '#336699', '#66CC66', '#FFCC33')) +
  ggtitle("Cell Type Marker Genes Expression")

ggsave("dotplot_pretty.pdf", dotplot_pretty, width = 12, height = 10)

# Save Seurat Object and CSVs
saveRDS(merged_seurat, "merged_seurat_auto_annotated.rds")
write.csv(cluster_scores, "cluster_celltype_scores.csv", row.names = TRUE)
write.csv(annotation_summary, "cluster_annotation_summary.csv", row.names = FALSE)
write.csv(cell_type_summary, "celltype_summary.csv", row.names = FALSE)

#=============================================================
# Step 15: Deep Validation Heatmaps
#=============================================================
cat("\n=== Generating Cell Annotation Validation Heatmaps ===\n")

if(!exists("markers") || !is.data.frame(markers)) {
  cat("⚠️ Warning: Markers object not found. Recalculating...\n")
  markers <- FindAllMarkers(merged_seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox", verbose = FALSE)
}

celltype_avg_exp <- AverageExpression(merged_seurat, features = present_markers, group.by = "cell_type_auto")$RNA
heatmap_matrix <- as.matrix(celltype_avg_exp)

pdf("08_celltype_annotation_validation_heatmap.pdf", width = 14, height = 12)
color_palette <- colorRampPalette(c("#00008B", "#0066CC", "#66CCFF", "#FFFFFF", "#FFCC66", "#FF6600", "#8B0000"))(100)
pheatmap(heatmap_matrix, color = color_palette, cluster_rows = TRUE, cluster_cols = FALSE, scale = "row",
         main = "Cell Type Annotation Validation", fontsize_row = 8, cellwidth = 25, cellheight = 12)
dev.off()

cat("\n=== Automated Annotation and Validation Complete! ===\n")
