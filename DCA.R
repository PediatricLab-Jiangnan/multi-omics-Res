# Decision curve analysis for the overall model
# =============================================

library(rmda)

work_dir <- "C:/Users/xieru/Desktop/Res/GEO"
setwd(work_dir)

model_obj <- readRDS("nomogram_model_data.rds")

model_data <- model_obj$model_data
model_glm <- model_obj$model_glm

cat("Preparing DCA data...\n")

cal_dca_data <- model_data
cal_dca_data$Predicted_Prob <- predict(
  model_glm,
  newdata = cal_dca_data,
  type = "response"
)

dca_data <- cal_dca_data[, c("Group", "Predicted_Prob")]
dca_data$Group <- as.numeric(dca_data$Group)

cat("Running DCA analysis...\n")

set.seed(123)

dca_overall <- decision_curve(
  Group ~ Predicted_Prob,
  data = dca_data,
  family = binomial(link = "logit"),
  thresholds = seq(0.01, 0.99, by = 0.01),
  confidence.intervals = 0.95,
  bootstraps = 1000,
  study.design = "cohort",
  fitted.risk = TRUE
)

png("DCA_curve_overall_model.png", width = 2400, height = 2000, res = 300)

plot_decision_curve(
  dca_overall,
  curve.names = "Overall Model",
  xlab = "Threshold Probability",
  ylab = "Net Benefit",
  standardize = FALSE,
  confidence.intervals = FALSE,
  col = "#A23B72",
  lwd = 2,
  legend.position = "topright"
)

title(main = "Decision Curve Analysis - Overall Model", cex.main = 1.2, font.main = 2)

dev.off()

pdf("DCA_curve_overall_model.pdf", width = 8, height = 7)

plot_decision_curve(
  dca_overall,
  curve.names = "Overall Model",
  xlab = "Threshold Probability",
  ylab = "Net Benefit",
  standardize = FALSE,
  confidence.intervals = FALSE,
  col = "#A23B72",
  lwd = 2,
  legend.position = "topright"
)

title(main = "Decision Curve Analysis - Overall Model", cex.main = 1.2, font.main = 2)

dev.off()

png("DCA_curve_overall_model_standardized.png", width = 2400, height = 2000, res = 300)

plot_decision_curve(
  dca_overall,
  curve.names = "Overall Model",
  xlab = "Threshold Probability",
  ylab = "Standardized Net Benefit",
  standardize = TRUE,
  confidence.intervals = FALSE,
  col = "#A23B72",
  lwd = 2,
  legend.position = "topright"
)

title(main = "Standardized Decision Curve Analysis - Overall Model", cex.main = 1.2, font.main = 2)

dev.off()

pdf("DCA_curve_overall_model_standardized.pdf", width = 8, height = 7)

plot_decision_curve(
  dca_overall,
  curve.names = "Overall Model",
  xlab = "Threshold Probability",
  ylab = "Standardized Net Benefit",
  standardize = TRUE,
  confidence.intervals = FALSE,
  col = "#A23B72",
  lwd = 2,
  legend.position = "topright"
)

title(main = "Standardized Decision Curve Analysis - Overall Model", cex.main = 1.2, font.main = 2)

dev.off()

write.csv(dca_overall$derived.data, "DCA_overall_model_results.csv", row.names = FALSE)
write.csv(cal_dca_data, "Overall_model_predicted_probabilities.csv", row.names = FALSE)

saveRDS(dca_overall, "DCA_results.rds")

cat("DCA analysis completed.\n")
