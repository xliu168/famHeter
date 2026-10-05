#!/usr/bin/env Rscript

# Run conventional and family-based Burden/SKAT/ACAT-O tests under the null.
#
# Usage:
#   Rscript simulation/02_run_null_tests.R <trait> <architecture> <hcr> [nrep] [causal_r2_reference]
#
# Examples:
#   Rscript simulation/02_run_null_tests.R continuous balanced 0.3 10000 0.015
#   Rscript simulation/02_run_null_tests.R binary     balanced 0.3 10000 0.03

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
if (length(args) < 3) {
  stop(
    "Usage: Rscript simulation/02_run_null_tests.R ",
    "<trait> <architecture> <hcr> [nrep] [causal_r2_reference]"
  )
}

trait <- validate_trait(args[1])
architecture <- validate_architecture(args[2])
hcr <- validate_hcr(args[3])
nrep <- if (length(args) >= 4) as.integer(args[4]) else 10000L
causal_r2 <- if (length(args) >= 5) as.numeric(args[5]) else default_causal_r2(trait)

repo_root <- get_repo_root()
input_dir <- file.path(repo_root, "results", "simulation_inputs")
null_dir <- file.path(repo_root, "results", "null", trait, architecture)

simfile <- file.path(
  null_dir,
  paste0(
    "null_", trait,
    "_nrep_", nrep,
    "_r2ref_", format(causal_r2, scientific = FALSE),
    ".rds"
  )
)
if (!file.exists(simfile)) {
  stop("Null phenotype file not found:\n", simfile)
}

aaf_file <- file.path(input_dir, paste0("heteroplasmy_hcr_", hcr, ".rds"))
required <- c(aaf_file, file.path(input_dir, "cov_sim.csv"),
              file.path(input_dir, "kin_n.rds"), file.path(input_dir, "kin_m2.rds"))
if (any(!file.exists(required))) {
  stop("Missing simulation inputs. Run simulation/00_generate_common_inputs.R first.")
}

simdata <- as.matrix(readRDS(simfile))
aaf <- readRDS(aaf_file)
cov_sim <- read.csv(file.path(input_dir, "cov_sim.csv"))
kin_n <- readRDS(file.path(input_dir, "kin_n.rds"))
kin_m2 <- readRDS(file.path(input_dir, "kin_m2.rds"))

IDs <- colnames(aaf)
dimnames(kin_n) <- list(IDs, IDs)
dimnames(kin_m2) <- list(IDs, IDs)

region <- prepare_region_inputs(aaf, gene = "MT-CYB")
seqMeta_cov_base <- as.data.frame(cov_sim)
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

outfile <- file.path(
  null_dir,
  paste0(
    "null_stats_", trait,
    "_hcr_", hcr,
    "_nrep_", nrep,
    "_r2ref_", format(causal_r2, scientific = FALSE),
    ".csv"
  )
)
write.csv(res, outfile, row.names = FALSE)

cat("Null test statistics saved to:\n", outfile, "\n")
