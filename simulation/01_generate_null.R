#!/usr/bin/env Rscript

# Generate null phenotypes for Type I error assessment / empirical null calibration.
#
# Usage:
#   Rscript simulation/01_generate_null.R <trait> <architecture> [nrep] [causal_r2_reference]
#
# Examples:
#   Rscript simulation/01_generate_null.R continuous balanced 10000 0.015
#   Rscript simulation/01_generate_null.R binary     balanced 10000 0.03
#
# causal_r2_reference is used only to allocate the background mitochondrial variance
# consistently with the corresponding power scenario. No causal heteroplasmy effect is
# included in these null phenotypes.

suppressPackageStartupMessages({
  library(MASS)
  library(foreach)
  library(doParallel)
})

source("R/simulation_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop(
    "Usage: Rscript simulation/01_generate_null.R ",
    "<trait> <architecture> [nrep] [causal_r2_reference]"
  )
}

trait <- validate_trait(args[1])
architecture <- validate_architecture(args[2])
nrep <- if (length(args) >= 3) as.integer(args[3]) else 10000L
causal_r2 <- if (length(args) >= 4) as.numeric(args[4]) else default_causal_r2(trait)

if (is.na(nrep) || nrep < 1) stop("nrep must be a positive integer.")
if (is.na(causal_r2) || causal_r2 <= 0) stop("causal_r2_reference must be > 0.")

v <- architecture_variances(architecture, causal_r2)

repo_root <- get_repo_root()
input_dir <- file.path(repo_root, "results", "simulation_inputs")
output_dir <- file.path(repo_root, "results", "null", trait, architecture)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
  "cov_sim.csv", "kin_n.rds", "kin_m.rds", "heteroplasmy_hcr_0.1.rds"
)
missing <- required[!file.exists(file.path(input_dir, required))]
if (length(missing) > 0) {
  stop(
    "Missing common simulation inputs:\n  ",
    paste(missing, collapse = "\n  "),
    "\nRun simulation/00_generate_common_inputs.R first."
  )
}

cov_sim <- data.matrix(read.csv(file.path(input_dir, "cov_sim.csv")))
kin_n <- readRDS(file.path(input_dir, "kin_n.rds"))
kin_m <- readRDS(file.path(input_dir, "kin_m.rds"))
aaf_ref <- readRDS(file.path(input_dir, "heteroplasmy_hcr_0.1.rds"))

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
    aaf = aaf_ref,
    cov_sim = cov_sim,
    cov_beta = c(0.057, 0.6),
    varM = v$varM,
    varN = v$varN,
    kin_m = kin_m,
    kin_n = kin_n,
    trait = trait,
    causal_fraction = 0,
    causal_sign = 1,
    causal_r2 = 0,
    intercept = 10,
    residual_var = 0.4,
    case_prob = 0.2
  )
}

doParallel::stopImplicitCluster()

colnames(simdata) <- colnames(aaf_ref)

outfile <- file.path(
  output_dir,
  paste0(
    "null_", trait,
    "_nrep_", nrep,
    "_r2ref_", format(causal_r2, scientific = FALSE),
    ".rds"
  )
)
saveRDS(simdata, outfile)

settings <- data.frame(
  trait = trait,
  architecture = architecture,
  nrep = nrep,
  causal_r2_reference = causal_r2,
  varM_background = v$varM,
  varN = v$varN,
  total_mt_target = v$total_mt,
  stringsAsFactors = FALSE
)
write.csv(
  settings,
  file.path(output_dir, "null_settings.csv"),
  row.names = FALSE
)

cat("Null phenotypes saved to:\n", outfile, "\n")
