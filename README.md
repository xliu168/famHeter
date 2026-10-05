# famHeter

Code accompanying the manuscript:

**Association Testing of Mitochondrial DNA Heteroplasmy Using Pedigree Data**

This repository implements family-based region-level association tests for mtDNA
heteroplasmy using generalized linear mixed models that can account for
pedigree-based nuclear relatedness, maternal-lineage correlation, or both.

## Repository contents

```text
famHeter/
├── R/
│   ├── prepScores_glmm_git.R
│   └── simulation_helpers.R
├── examples/
│   └── minimal_example.R
├── simulation/
│   ├── 00_generate_common_inputs.R
│   ├── 01_generate_null.R
│   ├── 02_run_null_tests.R
│   ├── 03_generate_power.R
│   ├── 04_run_power_tests.R
│   └── 99_save_session_info.R
├── REPRODUCIBILITY_NOTES.md
├── .gitignore
└── README.md
```

`prepScores_glmm_git.R` is the custom score-preparation function used by the
family-based Burden and SKAT analyses. `simulation_helpers.R` contains shared
simulation and testing functions, including the equal-weight ACAT-O combination.


## Example

From the repository root, run:

```bash
Rscript examples/minimal_example.R
```

The example uses a small synthetic pedigree and synthetic heteroplasmy data and
demonstrates:

- famBurden-B / famSKAT-B: adjustment for both nuclear and maternal-lineage correlation;
- famBurden-M / famSKAT-M: maternal-lineage correlation only;
- famBurden-N / famSKAT-N: pedigree-based nuclear correlation only;
- equal-weight ACAT-O combinations of the corresponding Burden and SKAT p-values.


## Main simulation design

The simulation code reproduces the manuscript framework using 750 four-member nuclear
families (3,000 individuals) and 1,000 potential heteroplasmic loci. Heteroplasmy
Concordance Rate (HCR) can be set to 0.1, 0.2, 0.3, 0.4, or 0.5.

The maternal-lineage matrix uses:

- `M_ij = 1` for individuals in the same short-span maternal lineage;
- `M_ij = 0` otherwise.

Three polygenic architectures are supported:

- `balanced`: total mitochondrial contribution 0.25, nuclear contribution 0.25;
- `nuclear_dominant`: total mitochondrial contribution 0.05, nuclear contribution 0.45;
- `mitochondrial_dominant`: total mitochondrial contribution 0.45, nuclear contribution 0.05.

For power simulations, 80% of observed heteroplasmic variants are selected as causal.
The causal heteroplasmy component explains 1.5% of variance for continuous traits and
3% of the underlying continuous liability for binary traits. Binary phenotypes are
created by dichotomizing the underlying continuous phenotype at the 80th percentile,
yielding approximately 20% prevalence.

## Simulation workflow

Run all commands from the repository root.

### 1. Generate common inputs

```bash
Rscript simulation/00_generate_common_inputs.R
```

### 2. Generate null phenotypes

```bash
Rscript simulation/01_generate_null.R continuous balanced 10000 0.015
Rscript simulation/01_generate_null.R binary balanced 10000 0.03
```

### 3. Run null tests

Example for HCR = 0.3:

```bash
Rscript simulation/02_run_null_tests.R continuous balanced 0.3 10000 0.015
Rscript simulation/02_run_null_tests.R binary balanced 0.3 10000 0.03
```

The output includes conventional Burden/SKAT, family-based Burden/SKAT under the
three correlation-adjustment strategies, and the four corresponding ACAT-O p-values.

### 4. Generate alternative phenotypes

The `causal_sign` argument gives the proportion of causal variants with the same
direction of effect. Manuscript settings include 1.00, 0.80, 0.70, 0.60, 0.55, and 0.50.

```bash
Rscript simulation/03_generate_power.R continuous balanced 0.3 1.00 10000 0.015
Rscript simulation/03_generate_power.R binary balanced 0.3 1.00 10000 0.03
```

### 5. Run empirical-size-corrected power analyses

```bash
Rscript simulation/04_run_power_tests.R continuous balanced 0.3 1.00 10000 0.015 0.001
Rscript simulation/04_run_power_tests.R binary balanced 0.3 1.00 10000 0.03 0.001
```

For Burden and SKAT, the method-specific critical value is the 99.9th percentile of
the corresponding empirical null test-statistic distribution at `alpha = 0.001`.
For ACAT-O, significance is represented by smaller p-values, so the corresponding
empirical threshold is the 0.1th percentile of the null ACAT-O p-value distribution.


## Parallel execution

On SGE/SCC, the scripts use `NSLOTS` when available. For local testing, they default
to at most four cores. The number of workers can also be set explicitly:

```bash
export SIM_NCORES=2
```

## Main R dependencies

The workflow uses `GMMAT`, `seqMeta`, `kinship2`, `Matrix`, `MASS`, `bdsmatrix`,
`foreach`, and `doParallel`.


