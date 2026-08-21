# Input data dictionary

Individual-level genotype and phenotype data are not distributed in this repository.
The scripts consume tab-delimited summary-statistic files with the schemas below.

## EXACT inputs

EXACT refines sex-specific, cell-type-specific sc-eQTL effects using information
from sex-stratified and sex-combined sc-eQTL summary statistics and bulk eQTL
summary statistics.

- eQTL beta and standard-error files are tab-delimited matrices with row names.
- The beta file contains the column `beta.est`.
- For Step 1 weight estimation, beta columns 3 through 16 must be the 14 cell-type
  beta columns. Their names must end in `_beta`; this suffix is removed in the
  output.
- Step 1 generates `weight_est_13par_chr1.txt` through
  `weight_est_13par_chr22.txt`, with identical columns in every file.
- Step 2 uses the across-chromosome mean of these weights to refine the matched
  sex-specific sc-eQTL effects.

## DRAMA inputs

- `beta.tsv`: first column `snp`; remaining columns are standardized effect estimates,
  with one column per sex-specific or sex-combined study.
- `se.tsv`: same rows and columns as `beta.tsv`, containing standard errors.
- `female_proportions.tsv`: columns `study` and `female_proportion`. Study names and
  order must correspond to the columns in `beta.tsv`. Use 1 and 0 for the directly
  observed female-only and male-only studies, respectively.

The analysis standardizes a study's effect and standard error as `z / sqrt(N_eff)`
and `1 / sqrt(N_eff)` before constructing these inputs.

Controlled or third-party data must be obtained from the repositories and under the
access terms listed in the manuscript's resource table.
