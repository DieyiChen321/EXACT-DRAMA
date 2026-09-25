# EXACT analysis code

This repository contains the analysis code accompanying the EXACT manuscript. It
provides two complementary methods that operate on different molecular and disease
association data:

1. **EXACT** (Elucidate seX-specific And Cell Type-specific regulatory effects)
   extends the JOBS framework to refine female- and male-specific single-cell eQTL
   (sc-eQTL) effects across cell types by borrowing information from large
   sex-combined bulk and single-cell eQTL datasets.
2. **DRAMA** jointly models sex-stratified and sex-combined GWAS summary statistics,
   accounting for the sex proportion of each study, to refine female- and
   male-specific GWAS effects and increase power for downstream colocalization and
   TWAS analyses.

EXACT and DRAMA are separate methods. Their refined results are brought together in
downstream sex-stratified analyses, such as GWAS-sc-eQTL colocalization.

## Repository contents

- `R/drama.R`: reusable DRAMA estimator.
- `scripts/run_drama.R`: command-line DRAMA analysis in 10 computational chunks.
- `scripts/prepare_drama_inputs.R`: constructs standardized DRAMA inputs from GWAS files.
- `scripts/estimate_jobs_weights.R`: EXACT Step 1, estimating chromosome-specific
  NNLS weights with JOBS.
- `scripts/run_jobs_eqtl.R`: EXACT Step 2, refining sex-specific sc-eQTL effects.
- `scripts/install_dependencies.R`: installs the R dependencies and JOBS package.
- `tests/smoke_test.R`: synthetic test that requires no controlled data.
- `examples/`: runnable EXACT and DRAMA simulations with inputs and expected outputs.
- `example_data/README.md`: input schemas and data-access boundary.

## Installation

The code was prepared for R 4.5. Run from the repository root:

```sh
Rscript scripts/install_dependencies.R
```

JOBS is installed from its public source repository. No login or personal
information is required.

## EXACT: refining sex-specific sc-eQTL effects

EXACT builds on JOBS, which models bulk eQTL effects as weighted combinations of
cell-type-specific sc-eQTL effects. EXACT additionally models the sex-combined
sc-eQTL effect for each cell type using its female- and male-specific effects. This
joint analysis improves estimation and detection of sex-stratified and
sex-difference sc-eQTLs.

### Step 1: estimate weights

Estimate weights separately for chromosomes 1 through 22. For chromosome 1, run:

```sh
Rscript scripts/estimate_jobs_weights.R 1 \
  data/jobs_input/sex_comb/beta \
  data/jobs_input/sex_comb/se \
  data/weights_est
```

Step 1 calls `JOBS::jobs.nnls.weights()` after reversing `beta.est` to align the
effect direction. It labels the returned weights using beta columns 3 through 16,
with the `_beta` suffix removed.

### Step 2: refine sex-specific sc-eQTL effects

After all 22 weight files have been generated, apply their across-chromosome mean
to the matched eQTL summary statistics for each chromosome:

```sh
Rscript scripts/run_jobs_eqtl.R \
  data/matched_eqtl_beta_chr1.txt \
  data/matched_eqtl_se_chr1.txt \
  data/weights_est \
  results/exact_refined_sceqtl_chr1.txt
```

Step 2 averages the weight estimates across chromosomes, reverses `beta.est` to
match the effect direction used in the manuscript, and calls
`JOBS::jobs.eqtls(..., COR = TRUE)` to obtain the EXACT-refined effects.

## DRAMA: refining sex-stratified GWAS effects

Construct aligned standardized beta and standard-error matrices and the female-study
proportion table for disease ID 1 through 12:

```sh
Rscript scripts/prepare_drama_inputs.R 1 \
  data/sex_stratified data/combined data/combined/N_eff_sex_ratio.txt data/prepared
```

The disease order is SLE, CD, PBC, RA, UC, VIT, AS, CELIAC, MS, PSOA, SJOGREN,
and T1D. The example above uses disease ID 1 and therefore creates
`sle_beta.tsv`, `sle_se.tsv`, and `sle_female_proportions.tsv` in
`data/prepared`. Then run one of the 10 computational chunks:

```sh
Rscript scripts/run_drama.R \
  data/prepared/sle_beta.tsv \
  data/prepared/sle_se.tsv \
  data/prepared/sle_female_proportions.tsv 1 \
  results/drama_chunk_1.tsv
```

Run chunk IDs 1 through 10 to reproduce all rows. Each study is represented by one
row of the design matrix, `(female proportion, 1 - female proportion)`. DRAMA uses
inverse-variance weighted generalized least squares and a Moore-Penrose inverse to
estimate female and male effects and their standard errors.

## Verification

The synthetic DRAMA test does not use manuscript data:

```sh
Rscript tests/smoke_test.R
```

Users can reproduce deterministic EXACT and DRAMA simulations, including
before-versus-after accuracy summaries, without manuscript data:

```sh
Rscript examples/run_all_examples.R
```

See `examples/README.md`. These examples demonstrate execution and improved
estimation accuracy under their stated simulations; they do not recreate the
manuscript's empirical power estimates.

## Data availability

This repository contains code only. It intentionally excludes controlled-access,
third-party, and large summary-statistic files. Source datasets, accession numbers,
and access conditions should be listed in the manuscript's key resources table and
data availability statement. Users can inspect and run the code without a
GitHub account.

## License and citation

Code in this repository is released under the MIT License. See `CITATION.cff` and
the associated manuscript for citation information.
