# Reviewer examples

These deterministic simulations let reviewers run both methods without controlled manuscript data. Inputs, known true effects, and expected outputs are committed.

```sh
Rscript scripts/install_dependencies.R
Rscript examples/run_all_examples.R
```

Reproduced files are written to `examples/reproduced_output/` and checked numerically against committed expected outputs.

The EXACT example contains 240 gene-SNP pairs across 14 cell types. Higher-precision bulk eQTL effects are weighted combinations of noisy female and male cell-type effects. The DRAMA example contains 300 variants and combines noisy sex-only GWAS with three more precise sex-combined studies having different female proportions.

The EXACT precision summary compares original and refined model-based standard errors. The DRAMA accuracy summary compares RMSE against known simulated truth. They illustrate the methods under the stated simulations; they do not recreate manuscript power estimates.
