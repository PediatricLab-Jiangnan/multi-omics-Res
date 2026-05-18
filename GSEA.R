# GSEA analysis based on limma-ranked genes
# =========================================

library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)

work_dir <- "C:/Users/xieru/Desktop/Res/GEO"
setwd(work_dir)

cat("Loading differential expression results...\n")
deg_results <- read.csv("differential_expression_results.csv", stringsAsFactors = FALSE)

# Remove duplicated ENTREZ IDs and missing values
deg_gsea <- deg_results[!is.na(deg_results$ENTREZID), ]
deg_gsea <- deg_gsea[!duplicated(deg_gsea$ENTREZID), ]

# Build ranked gene list using logFC
gene_list <- deg_gsea$logFC
names(gene_list) <- deg_gsea$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

cat("Running GO GSEA...\n")

gsea_go_bp <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

gsea_go_mf <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "MF",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

gsea_go_cc <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

cat("Running KEGG GSEA...\n")

gsea_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# Convert gene IDs to gene symbols where possible
if (nrow(as.data.frame(gsea_go_bp)) > 0) {
  gsea_go_bp <- setReadable(gsea_go_bp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

if (nrow(as.data.frame(gsea_go_mf)) > 0) {
  gsea_go_mf <- setReadable(gsea_go_mf, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

if (nrow(as.data.frame(gsea_go_cc)) > 0) {
  gsea_go_cc <- setReadable(gsea_go_cc, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

save_gsea_result <- function(gsea_obj, prefix, title) {
  gsea_df <- as.data.frame(gsea_obj)
  
  if (nrow(gsea_df) == 0) {
    cat("No significant GSEA result for", prefix, "\n")
    return(NULL)
  }
  
  write.csv(gsea_df, paste0(prefix, "_GSEA_results.csv"), row.names = FALSE)
  
  p_dot <- dotplot(gsea_obj, showCategory = 15, split = ".sign") +
    facet_grid(. ~ .sign) +
    labs(title = title)
  
  ggsave(paste0(prefix, "_GSEA_dotplot.png"), p_dot, width = 12, height = 8, dpi = 300)
  ggsave(paste0(prefix, "_GSEA_dotplot.pdf"), p_dot, width = 12, height = 8)
  
  p_ridge <- ridgeplot(gsea_obj, showCategory = 15) +
    labs(title = paste0(title, " Ridge Plot"))
  
  ggsave(paste0(prefix, "_GSEA_ridgeplot.png"), p_ridge, width = 12, height = 8, dpi = 300)
  ggsave(paste0(prefix, "_GSEA_ridgeplot.pdf"), p_ridge, width = 12, height = 8)
  
  top_id <- gsea_df$ID[1]
  
  p_curve <- gseaplot2(
    gsea_obj,
    geneSetID = top_id,
    title = gsea_df$Description[1]
  )
  
  ggsave(paste0(prefix, "_top_GSEA_curve.png"), p_curve, width = 10, height = 7, dpi = 300)
  ggsave(paste0(prefix, "_top_GSEA_curve.pdf"), p_curve, width = 10, height = 7)
}

save_gsea_result(gsea_go_bp, "GO_BP", "GO BP GSEA")
save_gsea_result(gsea_go_mf, "GO_MF", "GO MF GSEA")
save_gsea_result(gsea_go_cc, "GO_CC", "GO CC GSEA")
save_gsea_result(gsea_kegg, "KEGG", "KEGG GSEA")

saveRDS(
  list(
    gsea_go_bp = gsea_go_bp,
    gsea_go_mf = gsea_go_mf,
    gsea_go_cc = gsea_go_cc,
    gsea_kegg = gsea_kegg
  ),
  "GSEA_results.rds"
)

cat("GSEA analysis completed.\n")
