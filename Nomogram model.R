# Nomogram model construction using ferroptosis-related genes
# ==========================================================

library(rms)
library(regplot)
library(caret)

work_dir <- "C:/Users/xieru/Desktop/Res/GEO"
setwd(work_dir)

data_obj <- readRDS("preprocessed_expression_data.rds")

expr_log <- data_obj$expr_log
gene_symbols_filtered <- data_obj$gene_symbols_filtered
sample_groups <- data_obj$sample_groups

ferroptosis_genes <- c("SQSTM1", "KEAP1", "NFE2L2", "SLC7A11", "GPX4")

available_genes <- c()
gene_expr_matrix <- NULL

for (gene in ferroptosis_genes) {
  gene_idx <- which(gene_symbols_filtered == gene)
  
  if (length(gene_idx) > 0) {
    available_genes <- c(available_genes, gene)
    gene_expr <- expr_log[gene_idx[1], ]
    
    if (is.null(gene_expr_matrix)) {
      gene_expr_matrix <- data.frame(as.numeric(gene_expr))
      colnames(gene_expr_matrix) <- gene
    } else {
      gene_expr_matrix[[gene]] <- as.numeric(gene_expr)
    }
    
    cat("Found gene:", gene, "\n")
  } else {
    cat("Gene not found:", gene, "\n")
  }
}

if (length(available_genes) == 0) {
  stop("No target genes were found.")
}

model_data <- data.frame(
  Sample = names(sample_groups),
  Group = as.numeric(sample_groups == "Seizures"),
  gene_expr_matrix,
  stringsAsFactors = FALSE
)

rownames(model_data) <- model_data$Sample

set.seed(123)
train_index <- createDataPartition(model_data$Group, p = 0.7, list = FALSE)

train_data <- model_data[train_index, ]
test_data <- model_data[-train_index, ]

formula_lr <- as.formula(paste("Group ~", paste(available_genes, collapse = " + ")))

model_glm <- glm(formula_lr, data = train_data, family = binomial(link = "logit"))

dd <- datadist(train_data)
options(datadist = "dd")

model_lr <- lrm(formula_lr, data = train_data, x = TRUE, y = TRUE)

cat("Nomogram model summary:\n")
print(model_lr)

pdf("nomogram_regplot.pdf", width = 10, height = 8)

regplot(
  model_lr,
  observation = train_data[1, ],
  points = TRUE,
  odds = FALSE,
  showP = FALSE,
  rank = "sd",
  failtime = NULL,
  droplines = FALSE,
  interval = "confidence",
  clickable = FALSE,
  title = "Nomogram for Seizures Prediction"
)

dev.off()

png("nomogram_regplot.png", width = 3000, height = 2400, res = 300)

regplot(
  model_lr,
  observation = train_data[1, ],
  points = TRUE,
  odds = FALSE,
  showP = FALSE,
  rank = "sd",
  failtime = NULL,
  droplines = FALSE,
  interval = "confidence",
  clickable = FALSE,
  title = "Nomogram for Seizures Prediction"
)

dev.off()

saveRDS(
  list(
    model_data = model_data,
    train_data = train_data,
    test_data = test_data,
    model_glm = model_glm,
    model_lr = model_lr,
    available_genes = available_genes
  ),
  "nomogram_model_data.rds"
)

cat("Nomogram analysis completed.\n")
