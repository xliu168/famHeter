#!/usr/bin/env Rscript

# Generate the pedigree, heteroplasmy matrices, relationship matrices, and covariates
# used by both the null and power simulations.
#
# Run from the repository root:
#   Rscript simulation/00_generate_common_inputs.R

suppressPackageStartupMessages({
  library(kinship2)
  library(bdsmatrix)
})

source("R/simulation_helpers.R")

repo_root <- get_repo_root()
input_dir <- file.path(repo_root, "results", "simulation_inputs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

n_families <- 750L
n_loci <- 1000L
hcr_levels <- VALID_HCR

family_data <- make_family_data(n_families)
write.csv(
  family_data,
  file.path(input_dir, "family_data.csv"),
  row.names = FALSE
)

for (hcr in hcr_levels) {
  cat("Generating heteroplasmy matrix for HCR =", hcr, "\n")
  aaf <- simulate_heteroplasmy_HCR_sparse(
    family_data = family_data,
    hcr_target = hcr,
    n_loci = n_loci,
    prop_indiv = 0.03,
    seed = 11L
  )
  saveRDS(
    aaf,
    file.path(input_dir, paste0("heteroplasmy_hcr_", hcr, ".rds"))
  )
}

rel <- make_relationship_matrices(
  family_data,
  maternal_corr = 1
)
saveRDS(rel$kin_n, file.path(input_dir, "kin_n.rds"))
saveRDS(rel$kin_m, file.path(input_dir, "kin_m.rds"))
saveRDS(rel$kin_m2, file.path(input_dir, "kin_m2.rds"))

cov_sim <- make_covariates(family_data, seed = 123L)
write.csv(
  cov_sim,
  file.path(input_dir, "cov_sim.csv"),
  row.names = FALSE
)

cat("Common simulation inputs saved to:\n", input_dir, "\n")
