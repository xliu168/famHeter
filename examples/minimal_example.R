# Minimal runnable example for famHeter
#
# Purpose:
#   Demonstrate the family-based heteroplasmy testing workflow on a small
#   synthetic data set. This is a smoke test / usage example, not a simulation
#   designed to estimate type I error accurately.
#
# Run from the repository root:
#   Rscript examples/minimal_example.R

set.seed(20260922)

required <- c("GMMAT", "seqMeta", "kinship2", "Matrix", "MASS")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop(
    "Missing required packages: ", paste(missing, collapse = ", "),
    "\nInstall them before running this example."
  )
}

source(file.path("R", "prepScores_glmm_git.R"))

# 1. Small pedigree ---------------------------------------------------------
n_families <- 25L
family_size <- 4L
n <- n_families * family_size
ids <- as.character(1000L + seq_len(n))

family_data <- data.frame(
  ID = ids,
  famid = rep(seq_len(n_families), each = family_size),
  role = rep(c("father", "mother", "boy", "girl"), n_families),
  sex = rep(c(1, 0, 1, 0), n_families),
  fatherID = "0",
  motherID = "0",
  stringsAsFactors = FALSE
)

for (f in seq_len(n_families)) {
  ii <- which(family_data$famid == f)
  family_data$fatherID[ii[3:4]] <- family_data$ID[ii[1]]
  family_data$motherID[ii[3:4]] <- family_data$ID[ii[2]]
}

# Pedigree-based nuclear relationship matrix.
N <- 2 * kinship2::makekinship(
  family_data$famid,
  family_data$ID,
  family_data$fatherID,
  family_data$motherID
)
N <- as.matrix(N)
dimnames(N) <- list(ids, ids)

# Maternal-lineage relationship matrix.
# Individuals in the same short-span maternal lineage have M_ij = 1;
# otherwise M_ij = 0. Fathers form separate maternal lineages here.
M <- diag(n)
dimnames(M) <- list(ids, ids)
for (f in seq_len(n_families)) {
  ii <- which(family_data$famid == f)
  maternal <- ii[c(2, 3, 4)]
  M[maternal, maternal] <- 1
}

# 2. Synthetic heteroplasmy -------------------------------------------------
n_loci <- 40L
hcr <- 0.30
Z_locus_by_person <- matrix(
  0, nrow = n_loci, ncol = n,
  dimnames = list(paste0("v", seq_len(n_loci)), ids)
)

carrier_families <- sample(seq_len(n_families), 8L)

for (f in carrier_families) {
  ii <- which(family_data$famid == f)
  father <- ii[1]
  mother <- ii[2]
  children <- ii[3:4]

  father_loci <- sample(seq_len(n_loci), 5L)
  mother_loci <- sample(seq_len(n_loci), 6L)
  Z_locus_by_person[father_loci, father] <- 1
  Z_locus_by_person[mother_loci, mother] <- 1

  for (child in children) {
    n_child <- round(0.8 * length(mother_loci))
    n_shared <- round(hcr * (length(mother_loci) + n_child) / (1 + hcr))
    n_shared <- min(n_shared, length(mother_loci), n_child)

    shared <- sample(mother_loci, n_shared)
    Z_locus_by_person[shared, child] <- 1

    n_extra <- n_child - n_shared
    if (n_extra > 0) {
      available <- setdiff(seq_len(n_loci), mother_loci)
      extra <- sample(available, n_extra)
      Z_locus_by_person[extra, child] <- 1
    }
  }
}

# prepScores expects individuals in rows and variants in columns.
Z <- t(Z_locus_by_person)
keep <- apply(Z, 2, var) > 0
Z <- Z[, keep, drop = FALSE]

SNPInfo <- data.frame(
  Name = colnames(Z),
  gene = "MT-demo",
  stringsAsFactors = FALSE
)

# 3. Null phenotype ---------------------------------------------------------
dat <- data.frame(
  age = rnorm(n, mean = 58, sd = 5.89),
  sex = family_data$sex,
  row.names = ids
)

alpha_age <- 0.057
alpha_sex <- 0.6
intercept <- 10

b_M <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = 0.235 * M)
b_N <- MASS::mvrnorm(1, mu = rep(0, n), Sigma = 0.25 * N)
eps <- rnorm(n, sd = sqrt(0.4))

dat$y <- intercept + alpha_age * dat$age + alpha_sex * dat$sex +
  b_M + b_N + eps

# 4. Family-based tests -----------------------------------------------------
run_tests <- function(scores) {
  burden <- seqMeta::burdenMeta(
    scores,
    SNPInfo = SNPInfo,
    wts = 1,
    snpNames = "Name",
    aggregateBy = "gene",
    verbose = FALSE
  )
  skat <- seqMeta::skatMeta(
    scores,
    SNPInfo = SNPInfo,
    wts = 1,
    snpNames = "Name",
    aggregateBy = "gene",
    verbose = FALSE
  )

  c(
    burden_Z = burden$beta / burden$se,
    burden_p = burden$p,
    SKAT_Q = skat$Q,
    SKAT_p = skat$p
  )
}

scores_B <- prepScores(
  Z = Z,
  formula = y ~ age + sex,
  family = gaussian(),
  SNPInfo = SNPInfo,
  snpNames = "Name",
  aggregateBy = "gene",
  kins = N,
  mito.kins = M,
  sparse = TRUE,
  data = dat
)

scores_M <- prepScores(
  Z = Z,
  formula = y ~ age + sex,
  family = gaussian(),
  SNPInfo = SNPInfo,
  snpNames = "Name",
  aggregateBy = "gene",
  mito.kins = M,
  sparse = TRUE,
  data = dat
)

scores_N <- prepScores(
  Z = Z,
  formula = y ~ age + sex,
  family = gaussian(),
  SNPInfo = SNPInfo,
  snpNames = "Name",
  aggregateBy = "gene",
  kins = N,
  sparse = TRUE,
  data = dat
)

out <- rbind(
  famBoth = run_tests(scores_B),
  famMaternal = run_tests(scores_M),
  famNuclear = run_tests(scores_N)
)

print(round(out, 6))

# Equal-weight ACAT combination of famBurden and famSKAT p-values.
acat_combine <- function(pmat, eps = 1e-15) {
  pmat <- as.matrix(pmat)

  pmat[pmat < eps] <- eps
  pmat[pmat > 1 - eps] <- 1 - eps

  Q_stat <- rowSums(tan((0.5 - pmat) * pi)) / ncol(pmat)
  p_acat <- 0.5 - atan(Q_stat) / pi

  return(p_acat)
}

cat("\nACAT-O p-values:\n")
print(c(
  famBoth = acat_combine(matrix(out["famBoth", c("burden_p", "SKAT_p")], nrow = 1)),
  famMaternal = acat_combine(matrix(out["famMaternal", c("burden_p", "SKAT_p")], nrow = 1)),
  famNuclear = acat_combine(matrix(out["famNuclear", c("burden_p", "SKAT_p")], nrow = 1))
))

cat("\nSession information:\n")
sessionInfo()
