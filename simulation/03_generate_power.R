#!/usr/bin/env Rscript

# Generate alternative phenotypes for empirical power estimation.
#
# Usage:
#   Rscript simulation/03_generate_power.R <trait> <architecture> <hcr> <causal_sign> [nrep] [causal_r2]
#
# Examples:
#   Rscript simulation/03_generate_power.R continuous balanced 0.3 1.00 10000 0.015
#   Rscript simulation/03_generate_power.R binary     balanced 0.3 1.00 10000 0.03
#
# causal_sign is the proportion of causal heteroplasmic variants with the same
# (positive) direction of effect. Manuscript values include:
#   1.00, 0.80, 0.70, 0.60, 0.55, 0.50

suppressPackageStartupMessages({
  library(MASS)
  library(foreach)
  library(doParallel)
})

source("R/simulation_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop(
    "Usage: Rscript simulation/03_generate_power.R ",
    "<trait> <architecture> <hcr> <causal_sign> [nrep] [causal_r2]"
  )
}

trait <- validate_trait(args[1])
architecture <- validate_architecture(args[2])
hcr <- validate_hcr(args[3])
causal_sign <- as.numeric(args[4])
nrep <- if (length(args) >= 5) as.integer(args[5]) else 10000L
causal_r2 <- if (length(args) >= 6) as.numeric(args[6]) else default_causal_r2(trait)

if (is.na(causal_sign) || causal_sign < 0.5 || causal_sign > 1) {
  stop("causal_sign must be between 0.5 and 1.")
}
if (is.na(nrep) || nrep < 1) stop("nrep must be a positive integer.")
if (is.na(causal_r2) || causal_r2 <= 0) stop("causal_r2 must be > 0.")

v <- architecture_variances(architecture, causal_r2)

repo_root <- get_repo_root()
input_dir <- file.path(repo_root, "results", "simulation_inputs")
output_dir <- file.path(repo_root, "results", "power", trait, architecture)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

aaf_file <- file.path(input_dir, paste0("heteroplasmy_hcr_", hcr, ".rds"))
required <- c(aaf_file, file.path(input_dir, "cov_sim.csv"),
              file.path(input_dir, "kin_n.rds"), file.path(input_dir, "kin_m.rds"))
if (any(!file.exists(required))) {
  stop("Missing simulation inputs. Run simulation/00_generate_common_inputs.R first.")
}

aaf <- readRDS(aaf_file)
cov_sim <- data.matrix(read.csv(file.path(input_dir, "cov_sim.csv")))
kin_n <- readRDS(file.path(input_dir, "kin_n.rds"))
kin_m <- readRDS(file.path(input_dir, "kin_m.rds"))

ncores <- get_ncores()
doParallel::registerDoParallel(cores = ncores)
cat("Using", ncores, "workers.\n")

simdata <- foreach::foreach(
  i = seq_len(nrep),
  .combine = "rbind",
  .packages = "MASS",
  .export = "simulate_one_phenotype"
) %dopar% {
  simulate_one_phenotype(
    i = i,
    aaf = aaf,
    cov_sim = cov_sim,
    cov_beta = c(0.057, 0.6),
    varM = v$varM,
    varN = v$varN,
    kin_m = kin_m,
    kin_n = kin_n,
    trait = trait,
    causal_fraction = 0.8,
    causal_sign = causal_sign,
    causal_r2 = causal_r2,
    intercept = 10,
    residual_var = 0.4,
    case_prob = 0.2
  )
}

doParallel::stopImplicitCluster()

colnames(simdata) <- colnames(aaf)

tag <- paste0(
  "hcr_", hcr,
  "_sign_", format(causal_sign, nsmall = 2),
  "_nrep_", nrep,
  "_r2_", format(causal_r2, scientific = FALSE)
)

outfile <- file.path(output_dir, paste0("alternative_", trait, "_", tag, ".rds"))
saveRDS(simdata, outfile)

settings <- data.frame(
  trait = trait,
  architecture = architecture,
  hcr = hcr,
  causal_fraction = 0.8,
  causal_sign = causal_sign,
  causal_r2 = causal_r2,
  varM_background = v$varM,
  varN = v$varN,
  total_mt_target = v$total_mt,
  prevalence = if (trait == "binary") 0.2 else NA_real_,
  stringsAsFactors = FALSE
)
write.csv(
  settings,
  file.path(output_dir, paste0("settings_", tag, ".csv")),
  row.names = FALSE
)

cat("Alternative phenotypes saved to:\n", outfile, "\n")
