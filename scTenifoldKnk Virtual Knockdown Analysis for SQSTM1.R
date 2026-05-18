# scTenifoldKnk Virtual Knockdown Analysis for SQSTM1
# =============================================================

library(Seurat)
library(scTenifoldKnk)
library(dplyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(ggrepel)

# Set working directory
work_dir <- "C:/Users/xieru/Desktop/Xiaoxieandyueying"
setwd(work_dir)

# =============================================================
# Step 1: Load Seurat Object
# =============================================================
cat("========== Loading Seurat Object ==========\n")

if (!file.exists("merged_seurat_auto_annotated.rds")) {
    stop("Error: Input RDS file not found.")
}

merged_seurat <- readRDS("merged_seurat_auto_annotated.rds")
cat("Total Cells:", ncol(merged_seurat), "\n")
cat("Total Genes:", nrow(merged_seurat), "\n")

# =============================================================
# Step 2: Target Gene Validation
# =============================================================
cat("\n========== Validating Target Gene ==========\n")

# Standardize gene names to uppercase
rownames(merged_seurat) <- toupper(rownames(merged_seurat))

target_gene <- "SQSTM1"

if (target_gene %in% rownames(merged_seurat)) {
    cat("Target gene identified:", target_gene, "\n")
} else {
    stop("Critical Error: SQSTM1 not found in the dataset!")
}

# =============================================================
# Step 3: Feature Selection Using Highly Variable Genes
# =============================================================
cat("\n========== Feature Selection: Highly Variable Genes ==========\n")

# Use 3000 highly variable genes to improve regulatory network density
N_HVG <- 3000

if (length(VariableFeatures(merged_seurat)) == 0) {
    cat("Computing highly variable genes...\n")
    merged_seurat <- FindVariableFeatures(
        merged_seurat,
        nfeatures = N_HVG,
        verbose = FALSE
    )
}

hvgs <- VariableFeatures(merged_seurat)

# Force inclusion of SQSTM1 even if it is not highly variable
if (target_gene %in% hvgs) {
    selected_hvgs <- head(hvgs, N_HVG)
} else {
    cat("Note: Target gene is not in the top HVGs; adding it manually.\n")
    selected_hvgs <- c(head(hvgs, N_HVG - 1), target_gene)
}

# =============================================================
# Step 4: Global Virtual SQSTM1 Knockdown Analysis
# =============================================================
cat("\n========== Global SQSTM1 Knockdown Analysis ==========\n")

set.seed(2026)

# Sample up to 1000 cells for network construction
n_cells_to_sample <- min(1000, ncol(merged_seurat))
global_cells <- sample(colnames(merged_seurat), n_cells_to_sample)
seurat_global <- subset(merged_seurat, cells = global_cells)

# Check SQSTM1 expression sparsity
sqstm1_expr_count <- sum(
    GetAssayData(seurat_global, layer = "counts")[target_gene, ] > 0
)

cat("Cells expressing", target_gene, "in global sample:", sqstm1_expr_count, "\n")

if (sqstm1_expr_count < 10) {
    warning("Low SQSTM1 expression detected. Network inference may be unstable.")
}

# Construct count matrix
expr_matrix_global <- as.matrix(
    GetAssayData(seurat_global[selected_hvgs, ], layer = "counts")
)

cat("Building global regulatory network...\n")

perturbation_global <- scTenifoldKnk(
    countMatrix = expr_matrix_global,
    gKO = target_gene,
    qc = TRUE,
    nc_nNet = 10,
    nc_nCells = 500,
    nCores = parallel::detectCores()
)

write.csv(
    perturbation_global,
    "SQSTM1_global_knockdown_result.csv"
)

# =============================================================
# Step 5: Neuron-Specific Virtual SQSTM1 Knockdown Analysis
# =============================================================
cat("\n========== Neuron-Specific SQSTM1 Knockdown Analysis ==========\n")

# Use previously annotated cell types
if (!"cell_type_auto" %in% colnames(merged_seurat@meta.data)) {
    stop("Cell type annotations missing in metadata.")
}

neuron_cells <- WhichCells(
    merged_seurat,
    expression = cell_type_auto %in% c("Neurons", "Neuron", "Interneuron")
)

cat("Total neurons available:", length(neuron_cells), "\n")

if (length(neuron_cells) < 100) {
    stop("Insufficient neurons for robust network construction.")
}

n_neuron_sample <- min(1000, length(neuron_cells))

seurat_neuron <- subset(
    merged_seurat,
    cells = sample(neuron_cells, n_neuron_sample)
)

expr_matrix_neuron <- as.matrix(
    GetAssayData(seurat_neuron[selected_hvgs, ], layer = "counts")
)

cat("Building neuron-specific regulatory network...\n")

perturbation_neuron <- scTenifoldKnk(
    countMatrix = expr_matrix_neuron,
    gKO = target_gene,
    qc = TRUE,
    nc_nNet = 10,
    nc_nCells = min(300, ncol(expr_matrix_neuron)),
    nCores = parallel::detectCores()
)

write.csv(
    perturbation_neuron,
    "SQSTM1_neuron_specific_knockdown_result.csv"
)

# =============================================================
# Visualization Function
# =============================================================

generate_knockdown_plots <- function(result_df, prefix) {
    cat("\nGenerating report for:", prefix, "\n")
    
    # Ensure consistent data frame structure
    df <- as.data.frame(result_df)
    if (!"gene" %in% colnames(df)) df$gene <- rownames(df)
    
    # Standardize perturbation metrics
    if (!"FC" %in% colnames(df)) {
        df$FC <- if ("delta" %in% colnames(df)) 2^(df$delta) else 1
    }
    
    if (!"p.value" %in% colnames(df)) {
        df$p.value <- if ("zscore" %in% colnames(df)) {
            2 * pnorm(-abs(df$zscore))
        } else {
            1
        }
    }
    
    df <- df %>%
        filter(gene != target_gene) %>%
        mutate(log2FC = log2(FC))
    
    # Bar plot: top 20 perturbed genes
    top_genes <- df %>%
        arrange(desc(abs(log2FC))) %>%
        head(20)
    
    p1 <- ggplot(top_genes, aes(x = reorder(gene, log2FC), y = log2FC, fill = log2FC > 0)) +
        geom_bar(stat = "identity", alpha = 0.8) +
        scale_fill_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB")) +
        coord_flip() +
        labs(
            title = paste(prefix, "Top 20 Affected Genes"),
            subtitle = paste("Knockdown Target:", target_gene),
            x = "Gene",
            y = "log2(Fold Change)"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
    
    ggsave(paste0(prefix, "_SQSTM1_Barplot.pdf"), p1, width = 8, height = 7)
    
    # Volcano plot
    p2 <- ggplot(df, aes(x = log2FC, y = -log10(p.value))) +
        geom_point(aes(color = p.value < 0.05), alpha = 0.5) +
        scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#E74C3C")) +
        geom_text_repel(
            data = head(arrange(df, p.value), 10),
            aes(label = gene),
            size = 3
        ) +
        labs(
            title = paste(prefix, "SQSTM1 Perturbation Magnitude"),
            x = "log2(FC)",
            y = "-log10(P-value)"
        ) +
        theme_bw() +
        theme(legend.position = "none")
    
    ggsave(paste0(prefix, "_SQSTM1_Volcano.pdf"), p2, width = 7, height = 6)
}

# =============================================================
# Run Visualization
# =============================================================

generate_knockdown_plots(perturbation_global, "Global")
generate_knockdown_plots(perturbation_neuron, "Neuron")

cat("\n=============================================================\n")
cat("scTenifoldKnk SQSTM1 knockdown analysis complete.\n")
cat("Reports saved to:", getwd(), "\n")
cat("=============================================================\n")