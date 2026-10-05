# Shared helper functions for the famHeter simulation study.
# These functions were reorganized from the manuscript simulation scripts to
# remove SCC-specific paths and duplicated code while preserving the simulation logic.

suppressPackageStartupMessages({
  library(MASS)
  library(kinship2)
  library(bdsmatrix)
  library(seqMeta)
  library(GMMAT)
  library(Matrix)
})

VALID_HCR <- c(0.1, 0.2, 0.3, 0.4, 0.5)

get_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  farg <- grep("^--file=", args, value = TRUE)
  if (length(farg) > 0) {
    script_path <- sub("^--file=", "", farg[1])
    # Simulation scripts live one directory below repository root.
    return(normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE))
  }
  normalizePath(".", mustWork = TRUE)
}

get_ncores <- function() {
  n <- suppressWarnings(as.integer(Sys.getenv("SIM_NCORES")))
  if (!is.na(n) && n >= 1) return(n)

  n <- suppressWarnings(as.integer(Sys.getenv("NSLOTS")))
  if (!is.na(n) && n >= 1) return(n)

  detected <- parallel::detectCores()
  if (is.na(detected)) return(1L)
  max(1L, min(4L, detected - 1L))
}

validate_trait <- function(trait) {
  trait <- tolower(trait)
  if (!trait %in% c("continuous", "binary")) {
    stop("trait must be 'continuous' or 'binary'.")
  }
  trait
}

validate_hcr <- function(hcr) {
  hcr <- as.numeric(hcr)
  if (is.na(hcr) || !any(abs(hcr - VALID_HCR) < 1e-12)) {
    stop("HCR must be one of: ", paste(VALID_HCR, collapse = ", "))
  }
  hcr
}

validate_architecture <- function(architecture) {
  architecture <- tolower(architecture)
  aliases <- c(
    balanced = "balanced",
    nuclear = "nuclear_dominant",
    nuclear_dominant = "nuclear_dominant",
    mito = "mitochondrial_dominant",
    mitochondrial = "mitochondrial_dominant",
    mitochondrial_dominant = "mitochondrial_dominant"
  )
  if (!architecture %in% names(aliases)) {
    stop(
      "architecture must be one of: balanced, nuclear_dominant, ",
      "mitochondrial_dominant."
    )
  }
  unname(aliases[[architecture]])
}

default_causal_r2 <- function(trait) {
  trait <- validate_trait(trait)
  if (trait == "continuous") 0.015 else 0.03
}

architecture_variances <- function(architecture, causal_r2) {
  architecture <- validate_architecture(architecture)

  # The manuscript scenarios define total mitochondrial contribution and
  # nuclear contribution as follows. The causal heteroplasmy component is
  # part of the total mitochondrial contribution, so the remaining maternal-
  # lineage polygenic background is total_MT - causal_r2.
  vals <- switch(
    architecture,
    balanced = c(total_mt = 0.25, nuclear = 0.25),
    nuclear_dominant = c(total_mt = 0.05, nuclear = 0.45),
    mitochondrial_dominant = c(total_mt = 0.45, nuclear = 0.05)
  )

  maternal_background <- unname(vals["total_mt"]) - causal_r2
  if (maternal_background <= 0) {
    stop(
      "causal_r2 is too large for architecture '", architecture,
      "': total mitochondrial variance = ", vals["total_mt"], "."
    )
  }

  list(
    architecture = architecture,
    total_mt = unname(vals["total_mt"]),
    varN = unname(vals["nuclear"]),
    varM = maternal_background,
    causal_r2 = causal_r2
  )
}

make_family_data <- function(n_families = 750L) {
  n_indiv <- n_families * 4L
  IDs <- 1000 + seq_len(n_indiv)

  dat <- data.frame(
    ID = IDs,
    famid = rep(seq_len(n_families), each = 4),
    role = rep(c("father", "mother", "boy", "girl"), times = n_families),
    sex = rep(c(1, 0, 1, 0), times = n_families),
    fatherID = NA_character_,
    motherID = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_families)) {
    idx <- which(dat$famid == i)
    father_id <- dat$ID[idx[1]]
    mother_id <- dat$ID[idx[2]]

    dat$fatherID[idx[1:2]] <- "0"
    dat$motherID[idx[1:2]] <- "0"
    dat$fatherID[idx[3:4]] <- as.character(father_id)
    dat$motherID[idx[3:4]] <- as.character(mother_id)
  }

  dat
}

simulate_heteroplasmy_HCR_sparse <- function(
    family_data,
    hcr_target,
    n_loci = 1000L,
    prop_indiv = 0.03,
    seed = 11L) {

  hcr_target <- validate_hcr(hcr_target)
  n_indiv <- nrow(family_data)
  n_families <- length(unique(family_data$famid))
  id_index <- setNames(seq_len(nrow(family_data)), family_data$ID)

  mat <- matrix(0, nrow = n_loci, ncol = n_indiv)
  colnames(mat) <- family_data$ID

  n_fam_het <- ceiling((prop_indiv * n_indiv) / 4)
  set.seed(seed + round(hcr_target * 100))
  het_fam_ids <- sample(seq_len(n_families), n_fam_het, replace = FALSE)

  for (fam_id in het_fam_ids) {
    fam <- family_data[family_data$famid == fam_id, , drop = FALSE]

    idx_father <- id_index[[as.character(fam$ID[fam$role == "father"])]]
    idx_mother <- id_index[[as.character(fam$ID[fam$role == "mother"])]]
    idx_boy <- id_index[[as.character(fam$ID[fam$role == "boy"])]]
    idx_girl <- id_index[[as.character(fam$ID[fam$role == "girl"])]]

    het_counts <- sample(
      c(16, 17, 18),
      size = 2,
      replace = TRUE,
      prob = c(0.94, 0.04, 0.02)
    )

    father_loci <- sample(seq_len(n_loci), het_counts[1], replace = FALSE)
    mother_loci <- sample(seq_len(n_loci), het_counts[2], replace = FALSE)

    mat[father_loci, idx_father] <- 1
    mat[mother_loci, idx_mother] <- 1

    n_mother_het <- length(mother_loci)

    for (child_idx in c(idx_boy, idx_girl)) {
      child_het <- rep(0, n_loci)
      n_child_het <- round(0.8 * n_mother_het)

      n_concordant <- round(
        (hcr_target * (n_mother_het + n_child_het)) /
          (1 + hcr_target)
      )
      n_concordant <- min(n_concordant, n_mother_het, n_child_het)

      if (n_concordant > 0) {
        concordant_loci <- sample(mother_loci, n_concordant)
        child_het[concordant_loci] <- 1
      }

      n_nonconcordant <- n_child_het - n_concordant
      if (n_nonconcordant > 0) {
        candidate_loci <- setdiff(seq_len(n_loci), mother_loci)
        extra_loci <- sample(candidate_loci, n_nonconcordant)
        child_het[extra_loci] <- 1
      }

      mat[, child_idx] <- child_het
    }
  }

  mat
}

convert_nDNA_to_mtDNA <- function(nDNA_matrix, family_data, maternal_corr = 1) {
  # Individuals in the same short-span maternal lineage receive correlation 1,
  # matching the maternal-lineage matrix definition used in the manuscript.
  mtDNA_matrix <- as.matrix(nDNA_matrix)
  ids <- rownames(nDNA_matrix)

  for (i in seq_len(nrow(family_data))) {
    ind <- family_data$ID[i]
    father <- family_data$fatherID[i]
    mother <- family_data$motherID[i]

    ind_idx <- which(ids == ind)
    father_idx <- which(ids == father)
    mother_idx <- which(ids == mother)

    if (length(father_idx) > 0) {
      mtDNA_matrix[ind_idx, father_idx] <- 0
      mtDNA_matrix[father_idx, ind_idx] <- 0
    }

    if (length(mother_idx) > 0) {
      mtDNA_matrix[ind_idx, mother_idx] <- maternal_corr
      mtDNA_matrix[mother_idx, ind_idx] <- maternal_corr
    }
  }

  sibling_pairs <- which(
    mtDNA_matrix == 0.5 |
      mtDNA_matrix == 0.25 |
      mtDNA_matrix == 0.125 |
      mtDNA_matrix == 0.0625,
    arr.ind = TRUE
  )

  if (nrow(sibling_pairs) > 0) {
    for (pair in seq_len(nrow(sibling_pairs))) {
      row <- sibling_pairs[pair, 1]
      col <- sibling_pairs[pair, 2]
      if (row != col) {
        mtDNA_matrix[row, col] <- maternal_corr
        mtDNA_matrix[col, row] <- maternal_corr
      }
    }
  }

  mtDNA_matrix
}

make_relationship_matrices <- function(family_data, maternal_corr = 1) {
  kin_n <- 2 * kinship2::makekinship(
    family_data$famid,
    family_data$ID,
    family_data$fatherID,
    family_data$motherID
  )

  kin_m <- convert_nDNA_to_mtDNA(
    kin_n,
    family_data,
    maternal_corr = maternal_corr
  )

  n_families <- length(unique(family_data$famid))
  block_list <- lapply(seq_len(n_families), function(i) {
    idx <- ((i - 1) * 4 + 1):(i * 4)
    kin_m[idx, idx]
  })
  blocks_vector <- unlist(lapply(block_list, as.vector))
  kin_m2 <- bdsmatrix::bdsmatrix(
    blocksize = rep(4, n_families),
    blocks = blocks_vector
  )

  list(kin_n = kin_n, kin_m = kin_m, kin_m2 = kin_m2)
}

make_covariates <- function(family_data, seed = 123L) {
  set.seed(seed)
  cov_sim <- family_data[, c("ID", "sex")]
  cov_sim$age <- rnorm(nrow(cov_sim), 58, 5.89)
  rownames(cov_sim) <- cov_sim$ID
  cov_sim <- cov_sim[, c("age", "sex")]
  colnames(cov_sim) <- c("X1", "X2")
  data.matrix(cov_sim)
}

simulate_one_phenotype <- function(
    i,
    aaf,
    cov_sim,
    cov_beta = c(0.057, 0.6),
    varM,
    varN,
    kin_m,
    kin_n,
    trait = c("continuous", "binary"),
    causal_fraction = 0,
    causal_sign = 1,
    causal_r2 = 0,
    intercept = 10,
    residual_var = 0.4,
    case_prob = 0.2) {

  trait <- match.arg(trait)
  set.seed(10000 + i)

  IDs <- colnames(aaf)
  fixeff_cov <- as.vector(cov_sim %*% cov_beta)

  fixeff_mtgen <- rep(0, ncol(aaf))

  if (causal_fraction > 0) {
    rownames(aaf) <- seq_len(nrow(aaf))
    aaf_nonzero <- aaf[apply(aaf, 1, var, na.rm = TRUE) > 0, , drop = FALSE]
    eligible <- as.numeric(rownames(aaf_nonzero))

    n_causal <- round(length(eligible) * causal_fraction)
    if (n_causal < 1) stop("No causal variants were selected.")

    loci_causal <- sort(sample(eligible, size = n_causal, replace = FALSE))
    AAF_causal <- aaf_nonzero[as.character(loci_causal), , drop = FALSE]

    n_negative <- round(length(loci_causal) * (1 - causal_sign))
    sign_negative <- if (n_negative > 0) {
      sort(sample.int(length(loci_causal), size = n_negative))
    } else {
      integer(0)
    }

    V_sign <- rep(1, length(loci_causal))
    V_sign[sign_negative] <- -1

    D <- cor(t(AAF_causal))
    denom <- as.numeric(t(V_sign) %*% D %*% V_sign)
    if (!is.finite(denom) || denom <= 0) {
      stop("Invalid causal-variant correlation scaling denominator.")
    }
    C <- causal_r2 / denom

    loci_var <- apply(AAF_causal, 1, var)
    loci_var <- loci_var * (ncol(AAF_causal) - 1) / ncol(AAF_causal)
    if (any(!is.finite(loci_var)) || any(loci_var <= 0)) {
      stop("Invalid causal-variant variance encountered.")
    }

    beta_causal <- sqrt(C / loci_var)
    beta_causal[sign_negative] <- -beta_causal[sign_negative]

    fixeff_mtgen <- as.vector(t(AAF_causal) %*% beta_causal)
  }

  err <- rnorm(n = ncol(aaf), sd = sqrt(residual_var))
  M <- MASS::mvrnorm(
    1,
    mu = rep(0, length(IDs)),
    Sigma = varM * kin_m
  )
  N <- MASS::mvrnorm(
    1,
    mu = rep(0, length(IDs)),
    Sigma = varN * kin_n
  )

  latent <- fixeff_cov + fixeff_mtgen + M + N + intercept + err

  if (trait == "continuous") return(latent)

  cutoff <- quantile(latent, probs = 1 - case_prob)
  as.numeric(latent > cutoff)
}

build_snp_info <- function(aaf, gene = "MT-CYB") {
  data.frame(
    Name = rownames(aaf),
    gene = rep(gene, nrow(aaf)),
    stringsAsFactors = FALSE,
    row.names = rownames(aaf)
  )
}

prepare_region_inputs <- function(aaf, gene = "MT-CYB") {
  if (is.null(rownames(aaf))) rownames(aaf) <- as.character(seq_len(nrow(aaf)))

  keep <- apply(aaf, 1, var, na.rm = TRUE) > 0
  aaf2 <- aaf[keep, , drop = FALSE]
  snp_info <- build_snp_info(aaf2, gene = gene)

  list(
    Z = t(aaf2),
    SNPInfo = snp_info
  )
}


# Equal-weight ACAT combination used for ACAT-O and famACAT-O.
# This function is reproduced from the analysis code supplied by the study authors.
acat_combine <- function(pmat, eps = 1e-15) {
  pmat <- as.matrix(pmat)

  pmat[pmat < eps] <- eps
  pmat[pmat > 1 - eps] <- 1 - eps

  Q_stat <- rowSums(tan((0.5 - pmat) * pi)) / ncol(pmat)
  p_acat <- 0.5 - atan(Q_stat) / pi

  return(p_acat)
}

run_one_region_test <- function(
    i,
    simdata,
    seqMeta_cov_base,
    Z,
    SNPInfo,
    kin_m2,
    kin_n,
    family_type,
    prepScores_fun) {

  seqMeta_cov <- seqMeta_cov_base
  seqMeta_cov$y <- simdata[i, ]

  if (!family_type %in% c("gaussian", "binomial")) {
    stop("family_type must be 'gaussian' or 'binomial'.")
  }

  # seqMeta::prepScores2 expects a character string, whereas the custom
  # prepScores/GMMAT path expects an R family object.
  family_obj <- if (family_type == "gaussian") gaussian() else binomial()

  out <- setNames(rep(NA_real_, 20), c(
    "Burden_T", "SKAT_Q",
    "famBurden_T", "famSKAT_Q",
    "famBurdenM_T", "famSKATM_Q",
    "famBurdenN_T", "famSKATN_Q",
    "Burden_p", "SKAT_p",
    "famBurden_p", "famSKAT_p",
    "famBurdenM_p", "famSKATM_p",
    "famBurdenN_p", "famSKATN_p",
    "ACAT_p", "famACATB_p", "famACATM_p", "famACATN_p"
  ))

  tryCatch({
    scores0 <- seqMeta::prepScores2(
      Z = Z,
      formula = y ~ X1 + X2,
      family = family_type,
      SNPInfo = SNPInfo,
      snpNames = "Name",
      aggregateBy = "gene",
      sparse = TRUE,
      data = seqMeta_cov
    )

    burden0 <- seqMeta::burdenMeta(
      scores0,
      SNPInfo = SNPInfo,
      wts = 1,
      snpNames = "Name",
      aggregateBy = "gene",
      verbose = FALSE
    )
    skat0 <- seqMeta::skatMeta(
      scores0,
      SNPInfo = SNPInfo,
      wts = 1,
      snpNames = "Name",
      aggregateBy = "gene",
      verbose = FALSE
    )

    scoresB <- prepScores_fun(
      Z = Z,
      formula = y ~ X1 + X2,
      family = family_obj,
      SNPInfo = SNPInfo,
      snpNames = "Name",
      aggregateBy = "gene",
      kins = kin_n,
      mito.kins = kin_m2,
      sparse = TRUE,
      data = seqMeta_cov
    )
    burdenB <- seqMeta::burdenMeta(
      scoresB, SNPInfo = SNPInfo, wts = 1,
      snpNames = "Name", aggregateBy = "gene", verbose = FALSE
    )
    skatB <- seqMeta::skatMeta(
      scoresB, SNPInfo = SNPInfo, wts = 1,
      snpNames = "Name", aggregateBy = "gene", verbose = FALSE
    )

    scoresM <- prepScores_fun(
      Z = Z,
      formula = y ~ X1 + X2,
      family = family_obj,
      SNPInfo = SNPInfo,
      snpNames = "Name",
      aggregateBy = "gene",
      mito.kins = kin_m2,
      sparse = TRUE,
      data = seqMeta_cov
    )
    burdenM <- seqMeta::burdenMeta(
      scoresM, SNPInfo = SNPInfo, wts = 1,
      snpNames = "Name", aggregateBy = "gene", verbose = FALSE
    )
    skatM <- seqMeta::skatMeta(
      scoresM, SNPInfo = SNPInfo, wts = 1,
      snpNames = "Name", aggregateBy = "gene", verbose = FALSE
    )

    scoresN <- prepScores_fun(
      Z = Z,
      formula = y ~ X1 + X2,
      family = family_obj,
      SNPInfo = SNPInfo,
      snpNames = "Name",
      aggregateBy = "gene",
      kins = kin_n,
      sparse = TRUE,
      data = seqMeta_cov
    )
    burdenN <- seqMeta::burdenMeta(
      scoresN, SNPInfo = SNPInfo, wts = 1,
      snpNames = "Name", aggregateBy = "gene", verbose = FALSE
    )
    skatN <- seqMeta::skatMeta(
      scoresN, SNPInfo = SNPInfo, wts = 1,
      snpNames = "Name", aggregateBy = "gene", verbose = FALSE
    )

    p_acat <- acat_combine(cbind(burden0$p, skat0$p))
    p_acat_B <- acat_combine(cbind(burdenB$p, skatB$p))
    p_acat_M <- acat_combine(cbind(burdenM$p, skatM$p))
    p_acat_N <- acat_combine(cbind(burdenN$p, skatN$p))

    out[] <- c(
      burden0$beta / burden0$se, skat0$Q,
      burdenB$beta / burdenB$se, skatB$Q,
      burdenM$beta / burdenM$se, skatM$Q,
      burdenN$beta / burdenN$se, skatN$Q,
      burden0$p, skat0$p,
      burdenB$p, skatB$p,
      burdenM$p, skatM$p,
      burdenN$p, skatN$p,
      p_acat, p_acat_B, p_acat_M, p_acat_N
    )
  }, error = function(e) {
    message("Iteration ", i, " failed: ", conditionMessage(e))
  })

  out
}

statistic_columns <- function() {
  c(
    "Burden_T", "SKAT_Q",
    "famBurden_T", "famSKAT_Q",
    "famBurdenM_T", "famSKATM_Q",
    "famBurdenN_T", "famSKATN_Q"
  )
}

acat_p_columns <- function() {
  c("ACAT_p", "famACATB_p", "famACATM_p", "famACATN_p")
}
