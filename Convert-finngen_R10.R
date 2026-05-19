etwd("C:\\Users\\xieru\\Desktop\\Res\\MR")

input_file <- "finngen_R10_FE_MODE.csv"
output_file <- "finngen_R10_FE_MODE_outcome.csv"
outcome_name <- "FE"

# Optional: set this to a vector of SNP IDs, or keep NULL to export all SNPs.
# Example: snp_filter <- c("rs123", "rs456")
snp_filter <- NULL

# Optional: keep only variants with non-missing beta, SE, P value, SNP, and alleles.
drop_incomplete_required_rows <- TRUE

# ----------------------
# Package setup
# ----------------------

required_packages <- c("dplyr", "readr", "stringr")
installed_packages <- rownames(installed.packages())
install_missing <- required_packages[!required_packages %in% installed_packages]
if(length(install_missing) > 0) {
    install.packages(install_missing)
}

library(dplyr)
library(readr)
library(stringr)

# ----------------------
# Helper functions
# ----------------------

rename_first_match <- function(dat, new_name, candidates) {
    if(new_name %in% names(dat)) {
        return(dat)
    }
    
    matched <- candidates[candidates %in% names(dat)]
    if(length(matched) > 0) {
        names(dat)[names(dat) == matched[1]] <- new_name
    }
    
    dat
}

standardize_finngen_columns <- function(dat) {
    column_map <- list(
        SNP = c("SNP", "rsids", "rsid", "variant_id", "markername"),
        effect_allele.outcome = c("effect_allele.outcome", "alt", "effect_allele", "A1", "ea"),
        other_allele.outcome = c("other_allele.outcome", "ref", "other_allele", "A2", "nea"),
        beta.outcome = c("beta.outcome", "beta", "b", "effect"),
        se.outcome = c("se.outcome", "sebeta", "se", "std_error", "stderr"),
        pval.outcome = c("pval.outcome", "pval", "p_value", "p", "p.value"),
        eaf.outcome = c("eaf.outcome", "af_alt", "af_alt_cases", "af", "eaf"),
        ncase.outcome = c("ncase.outcome", "ncase", "cases", "n_cases"),
        ncontrol.outcome = c("ncontrol.outcome", "ncontrol", "controls", "n_controls"),
        samplesize.outcome = c("samplesize.outcome", "samplesize", "n", "N")
    )
    
    for(new_name in names(column_map)) {
        dat <- rename_first_match(dat, new_name, column_map[[new_name]])
    }
    
    dat
}

# ----------------------
# Read and convert FinnGen summary statistics
# ----------------------

message("Reading FinnGen file: ", input_file)
finngen_raw <- read_csv(input_file, show_col_types = FALSE)

message("Standardizing column names...")
outcome_dat <- standardize_finngen_columns(finngen_raw)

required_cols <- c(
    "SNP",
    "effect_allele.outcome",
    "other_allele.outcome",
    "beta.outcome",
    "se.outcome",
    "pval.outcome"
)

missing_cols <- required_cols[!required_cols %in% names(outcome_dat)]
if(length(missing_cols) > 0) {
    stop(
        "Missing required FinnGen columns after standardization: ",
        paste(missing_cols, collapse = ", "),
        "\nAvailable columns are: ",
        paste(names(finngen_raw), collapse = ", ")
    )
}

message("Cleaning SNP IDs and allele columns...")
outcome_dat <- outcome_dat %>%
    mutate(
        SNP = str_trim(as.character(SNP)),
        effect_allele.outcome = toupper(str_trim(as.character(effect_allele.outcome))),
        other_allele.outcome = toupper(str_trim(as.character(other_allele.outcome))),
        beta.outcome = as.numeric(beta.outcome),
        se.outcome = as.numeric(se.outcome),
        pval.outcome = as.numeric(pval.outcome)
    )

if("eaf.outcome" %in% names(outcome_dat)) {
    outcome_dat <- outcome_dat %>%
        mutate(eaf.outcome = as.numeric(eaf.outcome))
}

if("ncase.outcome" %in% names(outcome_dat)) {
    outcome_dat <- outcome_dat %>%
        mutate(ncase.outcome = as.numeric(ncase.outcome))
}

if("ncontrol.outcome" %in% names(outcome_dat)) {
    outcome_dat <- outcome_dat %>%
        mutate(ncontrol.outcome = as.numeric(ncontrol.outcome))
}

if(!"samplesize.outcome" %in% names(outcome_dat)) {
    if(all(c("ncase.outcome", "ncontrol.outcome") %in% names(outcome_dat))) {
        outcome_dat <- outcome_dat %>%
            mutate(samplesize.outcome = ncase.outcome + ncontrol.outcome)
    } else {
        outcome_dat <- outcome_dat %>%
            mutate(samplesize.outcome = NA_real_)
    }
} else {
    outcome_dat <- outcome_dat %>%
        mutate(samplesize.outcome = as.numeric(samplesize.outcome))
}

if(!is.null(snp_filter)) {
    message("Filtering to requested SNP list...")
    outcome_dat <- outcome_dat %>%
        filter(SNP %in% snp_filter)
}

if(drop_incomplete_required_rows) {
    message("Dropping rows with missing required MR fields...")
    outcome_dat <- outcome_dat %>%
        filter(
            !is.na(SNP),
            SNP != "",
            !is.na(effect_allele.outcome),
            !is.na(other_allele.outcome),
            !is.na(beta.outcome),
            !is.na(se.outcome),
            !is.na(pval.outcome)
        )
}

message("Adding TwoSampleMR outcome metadata columns...")
outcome_dat <- outcome_dat %>%
    mutate(
        outcome = outcome_name,
        id.outcome = outcome_name,
        mr_keep.outcome = TRUE
    )

# Keep standard MR columns first, while preserving any extra raw columns at the end.
standard_first_cols <- c(
    "SNP",
    "effect_allele.outcome",
    "other_allele.outcome",
    "beta.outcome",
    "se.outcome",
    "pval.outcome",
    "eaf.outcome",
    "samplesize.outcome",
    "ncase.outcome",
    "ncontrol.outcome",
    "outcome",
    "id.outcome",
    "mr_keep.outcome"
)

standard_first_cols <- standard_first_cols[standard_first_cols %in% names(outcome_dat)]
outcome_dat <- outcome_dat %>%
    select(all_of(standard_first_cols), everything())

message("Writing MR outcome CSV: ", output_file)
write_csv(outcome_dat, output_file)

message("Conversion complete.")
message("Rows written: ", nrow(outcome_dat))
message("Output file: ", output_file)
