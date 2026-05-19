# Set the working directory
setwd("C:\\Users\\xieru\\Desktop\\Res\\MR")

# Check and install required packages before loading them
required_packages <- c("TwoSampleMR", "dplyr", "ggplot2", "readr", "purrr", "tidyr")
installed_packages <- rownames(installed.packages())
install_missing <- required_packages[!required_packages %in% installed_packages]
if(length(install_missing) > 0) {
    install.packages(install_missing)
}

# Load required packages
library(TwoSampleMR)  # Core package for Mendelian randomization analysis
library(dplyr)        # Data manipulation
library(ggplot2)      # Plotting
library(readr)        # Fast file reading and writing
library(purrr)        # Functional programming helpers
library(tidyr)        # Data tidying

# ----------------------
# Part 2: Analysis parameters
# ----------------------

# Main parameters. Modify these according to your project.
analysis_settings <- list(
    exposure_path = "clumped_exposures",        # Directory containing clumped exposure CSV files
    outcome_file = "finngen_R10_FE_MODE_outcome.csv",     # Local outcome CSV file
    output_dir = "MR_Results_Local",            # Main output directory
    save_individual = TRUE,                     # Whether to save results for each exposure
    plot_width = 10,                            # Plot width in inches
    plot_height = 8,                            # Plot height in inches
    pval_threshold = 5e-8                       # P-value threshold used only for reporting significant results
)

# Create output directory structure
dir.create(analysis_settings$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(analysis_settings$output_dir, "Individual_Results"), showWarnings = FALSE)
dir.create(file.path(analysis_settings$output_dir, "Plots"), showWarnings = FALSE)
dir.create(file.path(analysis_settings$output_dir, "Summary_Tables"), showWarnings = FALSE)

# ----------------------
# Part 3: Load outcome data
# ----------------------
message("=========================================")
message("STEP 1: Loading outcome data")
message("=========================================")

# Read the preprocessed outcome data
message("Loading outcome data from: ", analysis_settings$outcome_file)

outcome_dat_full <- tryCatch({
    read_csv(analysis_settings$outcome_file, show_col_types = FALSE)
}, error = function(e) {
    stop("Failed to read outcome file: ", e$message)
})

# Check and standardize outcome column names
message("Checking outcome data format...")

# Required outcome columns for TwoSampleMR-compatible local analysis
required_outcome_cols <- c("SNP", "effect_allele.outcome", "other_allele.outcome", 
                           "beta.outcome", "se.outcome", "pval.outcome")

# If column names are not already standardized, try to detect and rename them automatically
if(!all(required_outcome_cols %in% colnames(outcome_dat_full))) {
    message("Standardizing outcome column names...")
    
    # Possible mappings from common raw column names to TwoSampleMR outcome column names
    col_mapping <- list(
        SNP = c("SNP", "rsids", "rsid", "snp"),
        effect_allele.outcome = c("effect_allele.outcome", "alt", "effect_allele", "A1"),
        other_allele.outcome = c("other_allele.outcome", "ref", "other_allele", "A2"),
        beta.outcome = c("beta.outcome", "beta", "b"),
        se.outcome = c("se.outcome", "sebeta", "se", "std_error"),
        pval.outcome = c("pval.outcome", "pval", "p_value", "p"),
        eaf.outcome = c("eaf.outcome", "af_alt_cases", "eaf", "af"),
        samplesize.outcome = c("samplesize.outcome", "samplesize", "n"),
        ncase.outcome = c("ncase.outcome", "ncase", "cases"),
        ncontrol.outcome = c("ncontrol.outcome", "ncontrol", "controls")
    )
    
    # Rename columns using the first matching candidate name
    for(new_name in names(col_mapping)) {
        for(old_name in col_mapping[[new_name]]) {
            if(old_name %in% colnames(outcome_dat_full)) {
                names(outcome_dat_full)[names(outcome_dat_full) == old_name] <- new_name
                break
            }
        }
    }
}

# Validate required columns
missing_cols <- required_outcome_cols[!required_outcome_cols %in% colnames(outcome_dat_full)]
if(length(missing_cols) > 0) {
    stop("Missing required columns in outcome data: ", paste(missing_cols, collapse = ", "))
}

# Add required identifier columns
outcome_dat_full <- outcome_dat_full %>%
    mutate(
        outcome = ifelse("outcome" %in% names(.), outcome, "FE"),
        id.outcome = outcome,
        mr_keep.outcome = TRUE
    )

# Print outcome data summary
message("Outcome data loaded successfully:")
message("  - Total SNPs: ", nrow(outcome_dat_full))
message("  - Outcome: ", unique(outcome_dat_full$outcome)[1])
if("samplesize.outcome" %in% names(outcome_dat_full)) {
    message("  - Sample size: ", unique(outcome_dat_full$samplesize.outcome)[1])
}
if(all(c("ncase.outcome", "ncontrol.outcome") %in% names(outcome_dat_full))) {
    message("  - Cases: ", unique(outcome_dat_full$ncase.outcome)[1])
    message("  - Controls: ", unique(outcome_dat_full$ncontrol.outcome)[1])
}

# ----------------------
# Part 4: Prepare exposure files
# ----------------------
message("\n=========================================")
message("STEP 2: Preparing exposure files")
message("=========================================")

# Get exposure file list. These files should already be clumped.
exposure_files <- list.files(
    path = analysis_settings$exposure_path,
    pattern = "\\.csv$",
    full.names = TRUE,
    ignore.case = TRUE
)

if(length(exposure_files) == 0) {
    stop("No exposure files found in: ", analysis_settings$exposure_path)
}

message("Found ", length(exposure_files), " exposure files")
message("Files to process:")
for(i in seq_along(exposure_files)) {
    message("  ", i, ". ", basename(exposure_files[i]))
}

# Initialize result containers
all_mr_results <- list()
all_heterogeneity <- list()
all_pleiotropy <- list()
analysis_summary <- data.frame()

# ----------------------
# Part 5: Main analysis loop
# ----------------------
message("\n=========================================")
message("STEP 3: Running MR analyses")
message("=========================================")

# Use safe execution so that one failed exposure does not stop the full batch
for(i in seq_along(exposure_files)) {
    file <- exposure_files[i]
    
    # Generate a unique analysis ID from the exposure file name
    analysis_id <- tools::file_path_sans_ext(basename(file))
    
    message("\n[", i, "/", length(exposure_files), "] Processing: ", analysis_id)
    
    # Create a subdirectory for the current exposure analysis
    indiv_dir <- file.path(analysis_settings$output_dir, "Individual_Results", analysis_id)
    dir.create(indiv_dir, showWarnings = FALSE, recursive = TRUE)
    
    # ----------------------
    # Step 1: Load clumped exposure data
    # ----------------------
    exposure_dat <- tryCatch({
        read_csv(file, show_col_types = FALSE)
    }, error = function(e) {
        message("  ERROR: Failed to read exposure file - ", e$message)
        return(NULL)
    })
    
    if(is.null(exposure_dat) || nrow(exposure_dat) == 0) {
        message("  SKIP: Empty exposure data")
        next
    }
    
    # Standardize exposure column names
    message("  - Standardizing exposure column names...")
    
    # Possible mappings from common exposure column names to TwoSampleMR exposure column names
    exp_col_mapping <- list(
        SNP = c("SNP", "rsids", "rsid", "snp"),
        effect_allele.exposure = c("effect_allele.exposure", "effect_allele", "alt", "A1"),
        other_allele.exposure = c("other_allele.exposure", "other_allele", "ref", "A2"),
        beta.exposure = c("beta.exposure", "beta", "b"),
        se.exposure = c("se.exposure", "se", "std_error"),
        pval.exposure = c("pval.exposure", "pval", "p_value", "p"),
        eaf.exposure = c("eaf.exposure", "eaf", "af", "af_alt")
    )
    
    for(new_name in names(exp_col_mapping)) {
        for(old_name in exp_col_mapping[[new_name]]) {
            if(old_name %in% colnames(exposure_dat)) {
                names(exposure_dat)[names(exposure_dat) == old_name] <- new_name
                break
            }
        }
    }
    
    # Check required exposure columns
    required_exp_cols <- c("SNP", "effect_allele.exposure", "other_allele.exposure", 
                           "beta.exposure", "se.exposure", "pval.exposure")
    missing_exp_cols <- required_exp_cols[!required_exp_cols %in% colnames(exposure_dat)]
    
    if(length(missing_exp_cols) > 0) {
        message("  ERROR: Missing required columns in exposure data: ", paste(missing_exp_cols, collapse = ", "))
        next
    }
    
    # Add required identifier columns
    exposure_dat <- exposure_dat %>%
        mutate(
            exposure = ifelse("exposure" %in% names(.), exposure, analysis_id),
            id.exposure = analysis_id,
            mr_keep.exposure = TRUE
        )
    
    message("  - Exposure SNPs: ", nrow(exposure_dat))
    
    # ----------------------
    # Step 2: Extract matching SNPs from outcome data
    # ----------------------
    message("  - Extracting matching SNPs from outcome data...")
    
    # Keep outcome rows for SNPs present in the exposure dataset
    outcome_dat <- outcome_dat_full %>%
        filter(SNP %in% exposure_dat$SNP)
    
    if(nrow(outcome_dat) == 0) {
        message("  WARNING: No matching SNPs found in outcome data")
        analysis_summary <- bind_rows(analysis_summary, data.frame(
            analysis_id = analysis_id,
            status = "Failed",
            reason = "No matching SNPs",
            n_snps_exposure = nrow(exposure_dat),
            n_snps_outcome = 0,
            n_snps_harmonised = 0,
            stringsAsFactors = FALSE
        ))
        next
    }
    
    message("  - Matching outcome SNPs: ", nrow(outcome_dat))
    
    # ----------------------
    # Step 3: Harmonise exposure and outcome data
    # ----------------------
    message("  - Harmonising data...")
    
    harmonised_dat <- tryCatch({
        harmonise_data(
            exposure_dat = exposure_dat,
            outcome_dat = outcome_dat
        )
    }, error = function(e) {
        message("  ERROR: Harmonisation failed - ", e$message)
        return(NULL)
    })
    
    if(is.null(harmonised_dat) || nrow(harmonised_dat) == 0) {
        message("  WARNING: Harmonisation produced no valid SNPs")
        analysis_summary <- bind_rows(analysis_summary, data.frame(
            analysis_id = analysis_id,
            status = "Failed",
            reason = "Harmonisation failed",
            n_snps_exposure = nrow(exposure_dat),
            n_snps_outcome = nrow(outcome_dat),
            n_snps_harmonised = 0,
            stringsAsFactors = FALSE
        ))
        next
    }
    
    # Keep only SNPs that TwoSampleMR marks as usable after harmonisation
    harmonised_dat <- harmonised_dat %>%
        filter(mr_keep == TRUE)
    
    if(nrow(harmonised_dat) == 0) {
        message("  WARNING: No valid SNPs after harmonisation")
        analysis_summary <- bind_rows(analysis_summary, data.frame(
            analysis_id = analysis_id,
            status = "Failed",
            reason = "No valid SNPs after harmonisation",
            n_snps_exposure = nrow(exposure_dat),
            n_snps_outcome = nrow(outcome_dat),
            n_snps_harmonised = 0,
            stringsAsFactors = FALSE
        ))
        next
    }
    
    message("  - Valid SNPs after harmonisation: ", nrow(harmonised_dat))
    
    # ----------------------
    # Step 4: Run MR analysis
    # ----------------------
    message("  - Running MR analysis...")
    
    # Select MR methods according to the number of available instruments
    if(nrow(harmonised_dat) == 1) {
        method_list <- c("mr_wald_ratio")
        message("  - Only 1 SNP available, using Wald ratio method")
    } else if(nrow(harmonised_dat) == 2) {
        method_list <- c("mr_ivw", "mr_egger_regression")
        message("  - 2 SNPs available, using IVW and MR-Egger")
    } else {
        method_list <- c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode")
        message("  - ", nrow(harmonised_dat), " SNPs available, using multiple methods")
    }
    
    mr_res <- tryCatch({
        mr(harmonised_dat, method_list = method_list)
    }, error = function(e) {
        message("  ERROR: MR analysis failed - ", e$message)
        return(NULL)
    })
    
    if(is.null(mr_res) || nrow(mr_res) == 0) {
        message("  WARNING: MR analysis produced no results")
        analysis_summary <- bind_rows(analysis_summary, data.frame(
            analysis_id = analysis_id,
            status = "Failed",
            reason = "MR analysis failed",
            n_snps_exposure = nrow(exposure_dat),
            n_snps_outcome = nrow(outcome_dat),
            n_snps_harmonised = nrow(harmonised_dat),
            stringsAsFactors = FALSE
        ))
        next
    }
    
    # Calculate odds ratios and 95% confidence intervals
    mr_res <- mr_res %>%
        mutate(
            or = exp(b),
            or_lci = exp(b - 1.96 * se),
            or_uci = exp(b + 1.96 * se),
            pval_formatted = format.pval(pval, digits = 3, eps = 0.001)
        )
    
    # ----------------------
    # Step 5: Run sensitivity analyses
    # ----------------------
    message("  - Running sensitivity analyses...")
    
    # Heterogeneity test
    het_res <- tryCatch({
        mr_heterogeneity(harmonised_dat)
    }, error = function(e) {
        message("    - Heterogeneity test failed: ", e$message)
        return(NULL)
    })
    
    # Horizontal pleiotropy test
    pleio_res <- tryCatch({
        mr_pleiotropy_test(harmonised_dat)
    }, error = function(e) {
        message("    - Pleiotropy test failed: ", e$message)
        return(NULL)
    })
    
    # Leave-one-out analysis requires at least 3 SNPs
    loo_res <- NULL
    if(nrow(harmonised_dat) >= 3) {
        loo_res <- tryCatch({
            mr_leaveoneout(harmonised_dat)
        }, error = function(e) {
            message("    - Leave-one-out analysis failed: ", e$message)
            return(NULL)
        })
    }
    
    # ----------------------
    # Step 6: Save results
    # ----------------------
    message("  - Saving results...")
    
    # Save harmonised data
    write_csv(harmonised_dat, file.path(indiv_dir, paste0(analysis_id, "_harmonised_data.csv")))
    
    # Save MR results
    write_csv(mr_res, file.path(indiv_dir, paste0(analysis_id, "_mr_results.csv")))
    
    # Save sensitivity analysis results
    if(!is.null(het_res)) {
        write_csv(het_res, file.path(indiv_dir, paste0(analysis_id, "_heterogeneity.csv")))
    }
    if(!is.null(pleio_res)) {
        write_csv(pleio_res, file.path(indiv_dir, paste0(analysis_id, "_pleiotropy.csv")))
    }
    if(!is.null(loo_res)) {
        write_csv(loo_res, file.path(indiv_dir, paste0(analysis_id, "_leaveoneout.csv")))
    }
    
    # Add results to combined containers
    all_mr_results[[analysis_id]] <- mr_res %>% mutate(analysis_id = analysis_id)
    if(!is.null(het_res)) {
        all_heterogeneity[[analysis_id]] <- het_res %>% mutate(analysis_id = analysis_id)
    }
    if(!is.null(pleio_res)) {
        all_pleiotropy[[analysis_id]] <- pleio_res %>% mutate(analysis_id = analysis_id)
    }
    
    # ----------------------
    # Step 7: Generate plots
    # ----------------------
    message("  - Generating plots...")
    
    # Scatter plot
    tryCatch({
        p1 <- mr_scatter_plot(mr_res, harmonised_dat)
        ggsave(
            file.path(indiv_dir, paste0(analysis_id, "_scatter.pdf")),
            p1[[1]],
            width = analysis_settings$plot_width,
            height = analysis_settings$plot_height,
            device = "pdf"
        )
    }, error = function(e) {
        message("    - Scatter plot failed: ", e$message)
    })
    
    # Forest plot
    if(nrow(mr_res) > 1) {
        tryCatch({
            p2 <- mr_forest_plot(mr_res)
            ggsave(
                file.path(indiv_dir, paste0(analysis_id, "_forest.pdf")),
                p2[[1]],
                width = analysis_settings$plot_width,
                height = max(4, nrow(mr_res) * 0.5),
                device = "pdf"
            )
        }, error = function(e) {
            message("    - Forest plot failed: ", e$message)
        })
    }
    
    # Funnel plot
    tryCatch({
        p3 <- mr_funnel_plot(mr_res)
        ggsave(
            file.path(indiv_dir, paste0(analysis_id, "_funnel.pdf")),
            p3[[1]],
            width = analysis_settings$plot_width,
            height = analysis_settings$plot_height,
            device = "pdf"
        )
    }, error = function(e) {
        message("    - Funnel plot failed: ", e$message)
    })
    
    # Leave-one-out plot
    if(!is.null(loo_res)) {
        tryCatch({
            p4 <- mr_leaveoneout_plot(loo_res)
            ggsave(
                file.path(indiv_dir, paste0(analysis_id, "_leaveoneout.pdf")),
                p4[[1]],
                width = analysis_settings$plot_width,
                height = analysis_settings$plot_height,
                device = "pdf"
            )
        }, error = function(e) {
            message("    - Leave-one-out plot failed: ", e$message)
        })
    }
    
    # Record successful analysis
    ivw_result <- mr_res %>% filter(method == "Inverse variance weighted")
    analysis_summary <- bind_rows(analysis_summary, data.frame(
        analysis_id = analysis_id,
        status = "Success",
        reason = NA,
        n_snps_exposure = nrow(exposure_dat),
        n_snps_outcome = nrow(outcome_dat),
        n_snps_harmonised = nrow(harmonised_dat),
        ivw_b = ifelse(nrow(ivw_result) > 0, ivw_result$b[1], NA),
        ivw_se = ifelse(nrow(ivw_result) > 0, ivw_result$se[1], NA),
        ivw_pval = ifelse(nrow(ivw_result) > 0, ivw_result$pval[1], NA),
        ivw_or = ifelse(nrow(ivw_result) > 0, ivw_result$or[1], NA),
        stringsAsFactors = FALSE
    ))
    
    message("  - Analysis completed successfully")
}

# ----------------------
# Part 6: Summarize results
# ----------------------
message("\n=========================================")
message("STEP 4: Summarizing results")
message("=========================================")

# Combine all available results
if(length(all_mr_results) > 0) {
    # Combined MR results
    combined_mr_results <- bind_rows(all_mr_results) %>%
        arrange(pval)
    write_csv(combined_mr_results, file.path(analysis_settings$output_dir, "Summary_Tables", "all_mr_results.csv"))
    
    # Combined heterogeneity results
    if(length(all_heterogeneity) > 0) {
        combined_heterogeneity <- bind_rows(all_heterogeneity)
        write_csv(combined_heterogeneity, file.path(analysis_settings$output_dir, "Summary_Tables", "all_heterogeneity.csv"))
    }
    
    # Combined pleiotropy results
    if(length(all_pleiotropy) > 0) {
        combined_pleiotropy <- bind_rows(all_pleiotropy)
        write_csv(combined_pleiotropy, file.path(analysis_settings$output_dir, "Summary_Tables", "all_pleiotropy.csv"))
    }
    
    # Analysis summary table
    write_csv(analysis_summary, file.path(analysis_settings$output_dir, "Summary_Tables", "analysis_summary.csv"))
    
    # ----------------------
    # Generate text report
    # ----------------------
    sink(file.path(analysis_settings$output_dir, "MR_Analysis_Report.txt"))
    
    cat("========================================\n")
    cat("Mendelian Randomization Analysis Summary Report - Local Version\n")
    cat("========================================\n\n")
    
    cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Outcome data:", analysis_settings$outcome_file, "\n")
    cat("Exposure data directory:", analysis_settings$exposure_path, "\n\n")
    
    cat("----------------------------------------\n")
    cat("Analysis overview\n")
    cat("----------------------------------------\n")
    cat("Total exposure files:", length(exposure_files), "\n")
    cat("Successful analyses:", sum(analysis_summary$status == "Success"), "\n")
    cat("Failed analyses:", sum(analysis_summary$status == "Failed"), "\n")
    cat("Success rate:", round(sum(analysis_summary$status == "Success") / length(exposure_files) * 100, 1), "%\n\n")
    
    # Summarize failure reasons
    if(sum(analysis_summary$status == "Failed") > 0) {
        cat("Failure reason summary:\n")
        failed_reasons <- analysis_summary %>%
            filter(status == "Failed") %>%
            group_by(reason) %>%
            summarise(count = n(), .groups = 'drop')
        for(i in 1:nrow(failed_reasons)) {
            cat("  - ", failed_reasons$reason[i], ": ", failed_reasons$count[i], "\n")
        }
        cat("\n")
    }
    
    # Significant results according to the reporting threshold
    sig_results <- combined_mr_results %>%
        filter(pval < analysis_settings$pval_threshold)
    
    cat("----------------------------------------\n")
    cat("Significant results (p < ", analysis_settings$pval_threshold, ")\n", sep="")
    cat("----------------------------------------\n")
    cat("Number of significant results:", nrow(sig_results), "\n\n")
    
    if(nrow(sig_results) > 0) {
        cat("Significant results by method:\n")
        sig_by_method <- sig_results %>%
            group_by(method) %>%
            summarise(
                count = n(),
                min_p = min(pval),
                .groups = 'drop'
            )
        print(sig_by_method, row.names = FALSE)
        cat("\n")
        
        cat("Top 20 most significant results:\n")
        print(
            sig_results %>%
                arrange(pval) %>%
                head(20) %>%
                select(analysis_id, method, nsnp, b, se, pval_formatted, or, or_lci, or_uci),
            row.names = FALSE
        )
    }
    
    cat("\n----------------------------------------\n")
    cat("IVW method result summary\n")
    cat("----------------------------------------\n")
    ivw_all <- combined_mr_results %>% filter(method == "Inverse variance weighted")
    if(nrow(ivw_all) > 0) {
        cat("Total IVW analyses:", nrow(ivw_all), "\n")
        cat("Significant IVW results (p < 0.05):", sum(ivw_all$pval < 0.05), "\n")
        cat("\nOR distribution:\n")
        cat("  - OR > 1 (risk factor):", sum(ivw_all$or > 1, na.rm = TRUE), "\n")
        cat("  - OR < 1 (protective factor):", sum(ivw_all$or < 1, na.rm = TRUE), "\n")
    }
    
    sink()
    
    # Generate summary visualization
    message("Generating summary plots...")
    
    # Volcano-style plot for IVW results
    if(nrow(ivw_all) > 0) {
        p_volcano <- ggplot(ivw_all, aes(x = b, y = -log10(pval))) +
            geom_point(aes(color = pval < 0.05), alpha = 0.6) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
            geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
            scale_color_manual(values = c("grey", "red")) +
            labs(x = "Beta (Effect size)", y = "-log10(P-value)", 
                 title = "MR Results Overview (IVW method)") +
            theme_minimal() +
            theme(legend.position = "none")
        
        ggsave(
            file.path(analysis_settings$output_dir, "Plots", "ivw_results_overview.pdf"),
            p_volcano,
            width = 10,
            height = 8,
            device = "pdf"
        )
    }
    
    message("\n=========================================")
    message("Analysis completed successfully!")
    message("=========================================")
    message("Results saved to: ", analysis_settings$output_dir)
    message("  - Individual results: ", file.path(analysis_settings$output_dir, "Individual_Results"))
    message("  - Summary tables: ", file.path(analysis_settings$output_dir, "Summary_Tables"))
    message("  - Plots: ", file.path(analysis_settings$output_dir, "Plots"))
    message("  - Full report: ", file.path(analysis_settings$output_dir, "MR_Analysis_Report.txt"))
    message("\nSuccessful analyses: ", sum(analysis_summary$status == "Success"), "/", length(exposure_files))
    
} else {
    message("\n=========================================")
    message("No analyses completed successfully!")
    message("=========================================")
    message("Please check:")
    message("  1. Exposure files are in the correct format")
    message("  2. Outcome file contains matching SNPs")
    message("  3. Allele coding is consistent between exposure and outcome")
}
