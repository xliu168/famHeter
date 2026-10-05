#!/usr/bin/env Rscript

# Run alternative tests and compute size-corrected empirical power.
# For Burden/SKAT, method-specific critical values are the 99.9th percentiles
# of the corresponding null test-statistic distributions when alpha = 0.001.
# For ACAT-O/famACAT-O, the equivalent lower-tail empirical threshold is the
# 0.1th percentile of the null ACAT p-value distribution.
#
# Usage:
#   Rscript simulation/04_run_power_tests.R <trait> <architecture> <hcr> <causal_sign> [nrep] [causal_r2] [alpha]
#
# Examples:
#   Rscript simulation/04_run_power_tests.R continuous balanced 0.3 1.00 10000 0.015 0.001
#   Rscript simulation/04_run_power_tests.R binary     balanced 0.3 1.00 10000 0.03  0.001

suppressPackageStartupMessages({
  library(seqMeta)
  library(GMMAT)
  library(Matrix)
  library(foreach)
  library(doParallel)
})

source("R/simulation_helpers.R")
source("R/prepScores_glmm_git.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(
    "Usage: Rscript simulation/04_run_power_tests.R ",
    "<trait> <architecture> <hcr> <causal_sign> [nrep] [causal_r2] [alpha]"
  )
}

trait <- validate_trait(args[1])
architecture <- validate_architecture(args[2])
hcr <- validate_hcr(args[3])
causal_sign <- as.numeric(args[4])
nrep <- if (length(args) >= 5) as.integer(args[5]) else 10000L
causal_r2 <- if (length(args) >= 6) as.numeric(args[6]) else default_causal_r2(trait)
alpha <- if (length(args) >= 7) as.numeric(args[7]) else 0.001

if (is.na(causal_sign) || causal_sign < 0.5 || causal_sign > 1) {
  stop("causal_sign must be between 0.5 and 1.")
}
if (is.na(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1.")

repo_root <- get_repo_root()
input_dir <- file.path(repo_root, "results", "simulation_inputs")
null_dir <- file.path(repo_root, "results", "null", trait, architecture)
power_dir <- file.path(repo_root, "results", "power", trait, architecture)
dir.create(power_dir, recursive = TRUE, showWarnings = FALSE)

null_stats_file <- file.path(
  null_dir,
  paste0(
    "null_stats_", trait,
    "_hcr_", hcr,
    "_nrep_", nrep,
    "_r2ref_", format(causal_r2, scientific = FALSE),
    ".csv"
  )
)

tag <- paste0(
  "hcr_", hcr,
  "_sign_", format(causal_sign, nsmall = 2),
  "_nrep_", nrep,
  "_r2_", format(causal_r2, scientific = FALSE)
)
alt_file <- file.path(power_dir, paste0("alternative_", trait, "_", tag, ".rds"))

if (!file.exists(null_stats_file)) {
  stop(
    "Null statistics file not found:\n", null_stats_file,
    "\nRun simulation/02_run_null_tests.R first."
  )
}
if (!file.exists(alt_file)) {
  stop(
    "Alternative phenotype file not found:\n", alt_file,
    "\nRun simulation/03_generate_power.R first."
  )
}

null_stats <- read.csv(null_stats_file, check.names = FALSE)
stat_cols <- statistic_columns()
acat_cols <- acat_p_columns()

if (!all(stat_cols %in% names(null_stats))) {
  stop("Null statistics file does not contain all expected Burden/SKAT statistic columns.")
}
if (!all(acat_cols %in% names(null_stats))) {
  stop(
    "Null statistics file does not contain ACAT-O columns. ",
    "Rerun simulation/02_run_null_tests.R with the current code."
  )
}

# Burden/SKAT: upper-tail empirical critical values.
empirical_critical_values <- vapply(
  null_stats[, stat_cols, drop = FALSE],
  quantile,
  numeric(1),
  probs = 1 - alpha,
  na.rm = TRUE,
  names = FALSE
)

# ACAT-O: smaller p-values are more significant, so use the lower alpha
# quantile of the empirical null p-value distribution.
empirical_acat_thresholds <- vapply(
  null_stats[, acat_cols, drop = FALSE],
  quantile,
  numeric(1),
  probs = alpha,
  na.rm = TRUE,
  names = FALSE
)

simdata <- as.matrix(readRDS(alt_file))
aaf <- readRDS(file.path(input_dir, paste0("heteroplasmy_hcr_", hcr, ".rds")))
cov_sim <- read.csv(file.path(input_dir, "cov_sim.csv"))
kin_n <- readRDS(file.path(input_dir, "kin_n.rds"))
kin_m2 <- readRDS(file.path(input_dir, "kin_m2.rds"))

IDs <- colnames(aaf)
dimnames(kin_n) <- list(IDs, IDs)
dimnames(kin_m2) <- list(IDs, IDs)

region <- prepare_region_inputs(aaf, gene = "MT-CYB")
seqMeta_cov_base <- as.data.frame(cov_sim)

# IMPORTANT: binary outcomes are fitted with binomial models here.
family_type <- if (trait == "continuous") "gaussian" else "binomial"

ncores <- get_ncores()
doParallel::registerDoParallel(cores = ncores)
cat("Using", ncores, "workers.\n")

res <- foreach::foreach(
  i = seq_len(nrep),
  .combine = "rbind",
  .packages = c("seqMeta", "GMMAT", "Matrix"),
  .export = c("prepScores", "run_one_region_test", "acat_combine")
) %dopar% {
  run_one_region_test(
    i = i,
    simdata = simdata,
    seqMeta_cov_base = seqMeta_cov_base,
    Z = region$Z,
    SNPInfo = region$SNPInfo,
    kin_m2 = kin_m2,
    kin_n = kin_n,
    family_type = family_type,
    prepScores_fun = prepScores
  )
}

doParallel::stopImplicitCluster()

stats_out <- file.path(power_dir, paste0("alternative_stats_", trait, "_", tag, ".csv"))
write.csv(res, stats_out, row.names = FALSE)

power_results <- vapply(stat_cols, function(test_name) {
  alt_stat <- res[, test_name]
  crit <- empirical_critical_values[[test_name]]
  mean(alt_stat > crit, na.rm = TRUE)
}, numeric(1))

acat_power_results <- vapply(acat_cols, function(test_name) {
  alt_p <- res[, test_name]
  crit_p <- empirical_acat_thresholds[[test_name]]
  mean(alt_p < crit_p, na.rm = TRUE)
}, numeric(1))

base_labels <- c(
  Burden_T = "Burden",
  SKAT_Q = "SKAT",
  famBurden_T = "famBurden-B",
  famSKAT_Q = "famSKAT-B",
  famBurdenM_T = "famBurden-M",
  famSKATM_Q = "famSKAT-M",
  famBurdenN_T = "famBurden-N",
  famSKATN_Q = "famSKAT-N"
)

acat_labels <- c(
  ACAT_p = "ACAT-O",
  famACATB_p = "famACAT-B-O",
  famACATM_p = "famACAT-M-O",
  famACATN_p = "famACAT-N-O"
)

power_df_base <- data.frame(
  Test = unname(base_labels[stat_cols]),
  InternalColumn = stat_cols,
  EmpiricalCriticalValue = as.numeric(empirical_critical_values[stat_cols]),
  CriticalValueType = "upper-tail test statistic",
  SizeCorrectedPower = as.numeric(power_results[stat_cols]),
  alpha = alpha,
  stringsAsFactors = FALSE
)

power_df_acat <- data.frame(
  Test = unname(acat_labels[acat_cols]),
  InternalColumn = acat_cols,
  EmpiricalCriticalValue = as.numeric(empirical_acat_thresholds[acat_cols]),
  CriticalValueType = "lower-tail ACAT p-value",
  SizeCorrectedPower = as.numeric(acat_power_results[acat_cols]),
  alpha = alpha,
  stringsAsFactors = FALSE
)

power_df <- rbind(power_df_base, power_df_acat)

power_out <- file.path(power_dir, paste0("power_", trait, "_", tag, ".csv"))
write.csv(power_df, power_out, row.names = FALSE)

cat("Alternative statistics saved to:\n", stats_out, "\n")
cat("Size-corrected power saved to:\n", power_out, "\n")
