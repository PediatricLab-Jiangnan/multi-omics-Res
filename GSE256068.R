# GEO data download, sample annotation, expression preprocessing, and PCA
# =====================================================================

library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(tibble)
library(org.Hs.eg.db)
library(AnnotationDbi)

work_dir <- "C:/Users/xieru/Desktop/Res/GEO/GSE256068"
setwd(work_dir)

cat("Downloading GEO dataset...\n")
gset <- getGEO("GSE256068", destdir = ".", getGPL = TRUE)

cat("Extracting expression set and sample metadata...\n")
eset <- if (is.list(gset)) gset[[1]] else gset
sample_info <- pData(eset)

sample_info$group <- ifelse(
  grepl("Control", sample_info$characteristics_ch1.3),
  "Control",
  "Seizures"
)

cat("Group summary:\n")
print(table(sample_info$group))

cat("Reading TPM expression matrix...\n")
expr_data <- read.table(
  "Norm_counts_TPM.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

original_sample_names <- colnames(expr_data)

cat("Annotating genes by ENTREZ ID...\n")
gene_ids <- rownames(expr_data)

gene_anno <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_ids,
  columns = "SYMBOL",
  keytype = "ENTREZID"
)

gene_anno_unique <- gene_anno[!duplicated(gene_anno$ENTREZID), ]

expr_data_anno <- expr_data
expr_data_anno$ENTREZID <- rownames(expr_data_anno)

expr_data_anno <- merge(
  gene_anno_unique,
  expr_data_anno,
  by.x = "ENTREZID",
  by.y = "ENTREZID",
  all.y = TRUE,
  sort = FALSE
)

rownames(expr_data_anno) <- expr_data_anno$ENTREZID

sample_columns <- setdiff(colnames(expr_data_anno), c("ENTREZID", "SYMBOL"))

sample_mapping <- sample_info[, c("geo_accession", "group")]
target_samples <- sample_mapping[sample_mapping$group %in% c("Control", "Seizures"), ]

available_samples <- intersect(target_samples$geo_accession, sample_columns)

cat("Matched samples:", length(available_samples), "\n")

expr_filtered <- expr_data_anno[, c("SYMBOL", available_samples)]

final_groups <- target_samples[target_samples$geo_accession %in% available_samples, ]
group_vector <- setNames(final_groups$group, final_groups$geo_accession)

cat("Final group summary:\n")
print(table(group_vector))

cat("Filtering low-expression genes...\n")
expr_matrix <- as.matrix(expr_filtered[, -1])
gene_symbols <- expr_filtered$SYMBOL

keep <- rowSums(expr_matrix > 1) >= 3
expr_matrix_filtered <- expr_matrix[keep, ]
gene_symbols_filtered <- gene_symbols[keep]

cat("Genes retained:", nrow(expr_matrix_filtered), "\n")

cat("Checking data distribution...\n")
max_value <- max(expr_matrix_filtered, na.rm = TRUE)

if (max_value > 50) {
  cat("Applying log2(TPM + 1) transformation...\n")
  expr_log <- log2(expr_matrix_filtered + 1)
  transformation_applied <- TRUE
} else {
  cat("Using original expression values...\n")
  expr_log <- expr_matrix_filtered
  transformation_applied <- FALSE
}

sample_groups <- group_vector[colnames(expr_log)]

cat("Generating expression distribution plot...\n")
set.seed(123)
sample_genes <- sample(nrow(expr_matrix_filtered), min(1000, nrow(expr_matrix_filtered)))

png("expression_distribution.png", width = 2400, height = 1200, res = 300)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

sample_data <- as.vector(expr_matrix_filtered[sample_genes, ])

hist(
  sample_data,
  breaks = 50,
  main = "Raw TPM Distribution",
  xlab = "TPM",
  ylab = "Frequency",
  col = "lightblue",
  border = "white"
)

hist(
  log2(sample_data + 1),
  breaks = 50,
  main = "Log2(TPM + 1) Distribution",
  xlab = "Log2(TPM + 1)",
  ylab = "Frequency",
  col = "lightcoral",
  border = "white"
)

dev.off()

cat("Running PCA analysis...\n")
expr_for_pca <- t(expr_log)
pca_result <- prcomp(expr_for_pca, scale. = TRUE)
variance_explained <- summary(pca_result)$importance[2, 1:2] * 100

pca_data <- data.frame(
  Sample = rownames(expr_for_pca),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = sample_groups,
  stringsAsFactors = FALSE
)

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  stat_ellipse(geom = "polygon", alpha = 0.2, color = NA) +
  stat_ellipse(size = 1, alpha = 0.8) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("Control" = "#2E86AB", "Seizures" = "#A23B72")) +
  scale_fill_manual(values = c("Control" = "#2E86AB", "Seizures" = "#A23B72")) +
  labs(
    title = "PCA Analysis of Gene Expression Data",
    subtitle = paste0("Samples: ", nrow(pca_data), " | Genes: ", nrow(expr_log)),
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)")
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray60"),
    legend.position = "bottom"
  )

print(pca_plot)

ggsave("PCA_analysis.png", pca_plot, width = 10, height = 8, dpi = 300)
ggsave("PCA_analysis.pdf", pca_plot, width = 10, height = 8)

expr_output <- data.frame(
  ENTREZID = rownames(expr_log),
  SYMBOL = gene_symbols_filtered,
  expr_log,
  stringsAsFactors = FALSE
)

write.csv(expr_output, "expression_matrix_processed.csv", row.names = FALSE)

sample_summary <- data.frame(
  Sample_ID = names(sample_groups),
  Group = sample_groups,
  stringsAsFactors = FALSE
)

write.csv(sample_summary, "sample_information.csv", row.names = FALSE)

saveRDS(
  list(
    expr_log = expr_log,
    gene_symbols_filtered = gene_symbols_filtered,
    sample_groups = sample_groups,
    transformation_applied = transformation_applied
  ),
  "preprocessed_expression_data.rds"
)

cat("GEO preprocessing completed.\n")
