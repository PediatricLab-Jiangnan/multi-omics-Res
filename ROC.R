# ROC curve analysis for overall model and individual genes
# ========================================================

library(pROC)
library(ggplot2)

work_dir <- "C:/Users/xieru/Desktop/Res/GEO"
setwd(work_dir)

model_obj <- readRDS("nomogram_model_data.rds")

model_data <- model_obj$model_data
model_glm <- model_obj$model_glm
available_genes <- model_obj$available_genes

cat("Calculating ROC for overall model...\n")

all_pred_prob <- predict(model_glm, newdata = model_data, type = "response")

roc_all <- roc(model_data$Group, all_pred_prob)
auc_all <- auc(roc_all)
ci_auc_all <- ci.auc(roc_all)

cat(sprintf(
  "Overall model AUC = %.3f, 95%% CI: %.3f-%.3f\n",
  auc_all, ci_auc_all[1], ci_auc_all[3]
))

single_gene_roc <- list()
single_gene_auc <- data.frame()

for (gene in available_genes) {
  single_model <- glm(
    as.formula(paste("Group ~", gene)),
    data = model_data,
    family = binomial
  )
  
  single_pred <- predict(single_model, type = "response")
  roc_single <- roc(model_data$Group, single_pred)
  auc_single <- auc(roc_single)
  ci_single <- ci.auc(roc_single)
  
  single_gene_roc[[gene]] <- roc_single
  
  single_gene_auc <- rbind(
    single_gene_auc,
    data.frame(
      Gene = gene,
      AUC = round(as.numeric(auc_single), 3),
      CI_lower = round(ci_single[1], 3),
      CI_upper = round(ci_single[3], 3)
    )
  )
}

single_gene_auc <- single_gene_auc[order(-single_gene_auc$AUC), ]
write.csv(single_gene_auc, "single_gene_ROC_AUC.csv", row.names = FALSE)

roc_plot_main <- ggroc(roc_all, color = "#A23B72", size = 1.5) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    title = "ROC Curve - Overall Model",
    subtitle = paste0(
      "AUC = ", round(auc_all, 3),
      " (95% CI: ", round(ci_auc_all[1], 3), "-",
      round(ci_auc_all[3], 3), ")"
    ),
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw()

ggsave("ROC_curve_overall.png", roc_plot_main, width = 8, height = 7, dpi = 300)
ggsave("ROC_curve_overall.pdf", roc_plot_main, width = 8, height = 7)

overall_label <- sprintf("Overall Model (AUC = %.3f)", auc_all)

roc_data <- data.frame(
  specificity = 1 - roc_all$specificities,
  sensitivity = roc_all$sensitivities,
  Model = overall_label
)

for (gene in available_genes) {
  roc_temp <- single_gene_roc[[gene]]
  gene_auc <- single_gene_auc$AUC[match(gene, single_gene_auc$Gene)]
  
  roc_data <- rbind(
    roc_data,
    data.frame(
      specificity = 1 - roc_temp$specificities,
      sensitivity = roc_temp$sensitivities,
      Model = sprintf("%s (AUC = %.3f)", gene, gene_auc)
    )
  )
}

roc_combined <- ggplot(roc_data, aes(x = specificity, y = sensitivity, color = Model)) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "gray50") +
  labs(
    title = "ROC Curves Comparison",
    x = "1 - Specificity",
    y = "Sensitivity",
    color = "Model"
  ) +
  theme_bw() +
  theme(legend.position = "right")

ggsave("ROC_curves_complete.png", roc_combined, width = 12, height = 7, dpi = 300)
ggsave("ROC_curves_complete.pdf", roc_combined, width = 12, height = 7)

saveRDS(
  list(
    roc_all = roc_all,
    auc_all = auc_all,
    ci_auc_all = ci_auc_all,
    single_gene_roc = single_gene_roc,
    single_gene_auc = single_gene_auc,
    all_pred_prob = all_pred_prob
  ),
  "ROC_results.rds"
)

cat("ROC analysis completed.\n")
