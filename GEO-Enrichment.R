# Differential expression, volcano plot, GO/KEGG enrichment, and selected gene boxplots
# ==================================================================================

library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(patchwork)

work_dir <- "C:/Users/xieru/Desktop/Res/GEO"
setwd(work_dir)

data_obj <- readRDS("preprocessed_expression_data.rds")

expr_log <- data_obj$expr_log
gene_symbols_filtered <- data_obj$gene_symbols_filtered
sample_groups <- data_obj$sample_groups
transformation_applied <- data_obj$transformation_applied

cat("Running differential expression analysis...\n")

group_factor <- factor(sample_groups, levels = c("Control", "Seizures"))

design <- model.matrix(~0 + group_factor)
colnames(design) <- c("Control", "Seizures")

contrast_matrix <- makeContrasts(
  Seizures_vs_Control = Seizures - Control,
  levels = design
)

fit <- lmFit(expr_log, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

deg_results <- topTable(
  fit2,
  coef = "Seizures_vs_Control",
  number = Inf,
  sort.by = "P"
)

deg_results$ENTREZID <- rownames(deg_results)

deg_results <- merge(
  deg_results,
  data.frame(
    ENTREZID = rownames(expr_log),
    SYMBOL = gene_symbols_filtered
  ),
  by = "ENTREZID",
  all.x = TRUE
)

deg_results <- deg_results[!is.na(deg_results$SYMBOL), ]
deg_results <- deg_results[order(deg_results$P.Value), ]

logFC_threshold <- 0.58
pvalue_threshold <- 0.05
padj_threshold <- 0.05

deg_results$Change <- "Not Significant"
deg_results$Change[
  deg_results$logFC > logFC_threshold & deg_results$P.Value < pvalue_threshold
] <- "Up-regulated"
deg_results$Change[
  deg_results$logFC < -logFC_threshold & deg_results$P.Value < pvalue_threshold
] <- "Down-regulated"

deg_results$Change_strict <- "Not Significant"
deg_results$Change_strict[
  deg_results$logFC > logFC_threshold & deg_results$adj.P.Val < padj_threshold
] <- "Up-regulated"
deg_results$Change_strict[
  deg_results$logFC < -logFC_threshold & deg_results$adj.P.Val < padj_threshold
] <- "Down-regulated"

cat("DEG summary:\n")
print(table(deg_results$Change))

write.csv(deg_results, "differential_expression_results.csv", row.names = FALSE)
write.csv(deg_results[deg_results$Change == "Up-regulated", ], "upregulated_genes.csv", row.names = FALSE)
write.csv(deg_results[deg_results$Change == "Down-regulated", ], "downregulated_genes.csv", row.names = FALSE)

cat("Drawing volcano plot...\n")

volcano_data <- deg_results
volcano_data$neg_log10_pval <- -log10(volcano_data$P.Value)

top_up_genes <- head(volcano_data[volcano_data$Change == "Up-regulated", ], 5)
top_down_genes <- head(volcano_data[volcano_data$Change == "Down-regulated", ], 5)
top_genes <- rbind(top_up_genes, top_down_genes)

volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_pval)) +
  geom_point(aes(color = Change), alpha = 0.6, size = 1.2) +
  scale_color_manual(
    values = c(
      "Up-regulated" = "#E31A1C",
      "Down-regulated" = "#1F78B4",
      "Not Significant" = "gray70"
    )
  ) +
  geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "gray40") +
  geom_text_repel(
    data = top_genes,
    aes(label = SYMBOL),
    size = 3,
    max.overlaps = 15
  ) +
  labs(
    title = "Volcano Plot: Seizures vs Control",
    subtitle = paste0(
      "Up-regulated: ", sum(volcano_data$Change == "Up-regulated"),
      " | Down-regulated: ", sum(volcano_data$Change == "Down-regulated")
    ),
    x = "Log2 Fold Change",
    y = "-Log10 P-value",
    color = "Regulation"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(volcano_plot)

ggsave("volcano_plot.png", volcano_plot, width = 10, height = 8, dpi = 300)
ggsave("volcano_plot.pdf", volcano_plot, width = 10, height = 8)

cat("Running GO and KEGG enrichment analysis...\n")

sig_genes <- deg_results[deg_results$Change != "Not Significant", ]
sig_entrez <- sig_genes$ENTREZID[!is.na(sig_genes$ENTREZID)]

go_bp <- enrichGO(sig_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)

go_mf <- enrichGO(sig_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)

go_cc <- enrichGO(sig_entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)

kegg <- enrichKEGG(
  gene = sig_entrez,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

if (nrow(as.data.frame(kegg)) > 0) {
  kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

save_enrichment_plot <- function(enrich_obj, prefix, title) {
  enrich_df <- as.data.frame(enrich_obj)
  
  if (nrow(enrich_df) == 0) {
    cat("No significant result for", prefix, "\n")
    return(NULL)
  }
  
  write.csv(enrich_df, paste0(prefix, "_enrichment.csv"), row.names = FALSE)
  
  p <- dotplot(enrich_obj, showCategory = 15, title = title) +
    theme(axis.text.y = element_text(size = 10))
  
  print(p)
  ggsave(paste0(prefix, "_dotplot.png"), p, width = 12, height = 8, dpi = 300)
  ggsave(paste0(prefix, "_dotplot.pdf"), p, width = 12, height = 8)
}

save_enrichment_plot(go_bp, "GO_BP", "GO Biological Process")
save_enrichment_plot(go_mf, "GO_MF", "GO Molecular Function")
save_enrichment_plot(go_cc, "GO_CC", "GO Cellular Component")
save_enrichment_plot(kegg, "KEGG", "KEGG Pathway Enrichment")

cat("Drawing selected gene expression boxplots...\n")

target_genes <- c("GPX4", "SQSTM1", "IL6", "IL6R", "KEAP1", "SLC7A11", "FTH1", "NFE2L2", "NFKB1", "NFKB2")

create_gene_boxplot <- function(gene_symbol) {
  gene_idx <- which(gene_symbols_filtered == gene_symbol)
  
  if (length(gene_idx) == 0) {
    cat("Gene not found:", gene_symbol, "\n")
    return(NULL)
  }
  
  gene_expr <- expr_log[gene_idx[1], ]
  
  plot_data <- data.frame(
    Sample = names(gene_expr),
    Expression = as.numeric(gene_expr),
    Group = sample_groups[names(gene_expr)]
  )
  
  t_result <- t.test(Expression ~ Group, data = plot_data)
  p_text <- ifelse(t_result$p.value < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", t_result$p.value)))
  
  p <- ggplot(plot_data, aes(x = Group, y = Expression, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    geom_jitter(aes(color = Group), width = 0.2, size = 2.5, alpha = 0.8) +
    scale_fill_manual(values = c("Control" = "#4DBBD5", "Seizures" = "#E64B35")) +
    scale_color_manual(values = c("Control" = "#0073C2", "Seizures" = "#BC3C29")) +
    annotate("text", x = 1.5, y = max(plot_data$Expression) * 1.05, label = p_text, size = 4, fontface = "bold") +
    labs(
      title = paste0("Expression of ", gene_symbol),
      x = "Group",
      y = ifelse(transformation_applied, "Log2(TPM + 1)", "TPM")
    ) +
    theme_bw() +
    theme(legend.position = "none")
  
  return(p)
}

gene_plots <- list()
valid_genes <- c()

for (gene in target_genes) {
  p <- create_gene_boxplot(gene)
  if (!is.null(p)) {
    gene_plots[[gene]] <- p
    valid_genes <- c(valid_genes, gene)
    ggsave(paste0("boxplot_", gene, ".png"), p, width = 8, height = 6, dpi = 300)
    ggsave(paste0("boxplot_", gene, ".pdf"), p, width = 8, height = 6)
  }
}

if (length(valid_genes) > 1) {
  combined_plot <- wrap_plots(gene_plots[valid_genes], ncol = 3) +
    plot_annotation(title = "Expression Analysis of Selected Genes")
  
  ggsave("combined_boxplots.png", combined_plot, width = 12, height = 10, dpi = 300)
  ggsave("combined_boxplots.pdf", combined_plot, width = 12, height = 10)
}

cat("Gene expression analysis completed.\n")
