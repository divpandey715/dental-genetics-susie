# Dental Genetics SuSiE Fine-Mapping

This repository contains an R workflow for SuSiE-RSS and mvSuSiE-RSS fine-mapping of dental GWAS traits using summary statistics and locus-specific LD matrices.

Results are published as an interactive [workflowr](https://workflowr.github.io/workflowr/) site: https://divpandey715.github.io/dental-genetics-susie/

GWAS summary statistics come from Shungin et al., ["Genome-wide analysis of dental caries and periodontitis combining clinical and self-reported data,"](https://www.nature.com/articles/s41467-019-10630-1) *Nature Communications*, 2019.

## Traits

- **DMFS** (Decayed, Missing, and Filled Surfaces)
- **Nteeth** (number of teeth)
- **DFSS** (Decayed and Filled Surfaces)

## Pipeline

Scripts are numbered and should be run from the repository root, in order. Steps 1-3 (single-trait) can be run in any order relative to each other; the joint mvSuSiE steps (4-7) each depend only on their own trait pair's single-trait outputs, not on each other; step 8 must run last.

**1. Single-trait SuSiE-RSS**

| Trait | Fine-mapping | Summarize | Position mapping | eQTL colocalization |
|---|---|---|---|---|
| DMFS | `02_run_susie.R` | `03_summarize_results.R` | `04_position_mapping.R` | `05_eqtl_colocalization.R` |
| Nteeth | `02_run_susie_nteeth.R` | `03_summarize_results_nteeth.R` | `04_position_mapping_nteeth.R` | `05_eqtl_colocalization_nteeth.R` |
| DFSS | `08_run_susie_dfss.R` | `08b_summarize_results_dfss.R` | `09_position_mapping_dfss.R` | `10_dfss_eqtl_analysis.R` |

`01_prepare_input_data.R` converts the raw GWAS summary statistic `.txt` files (see [Input files](#input-files)) into `.RDS` files used by all of the above. `06_display_results.R` is a legacy script for previewing result tables outside the workflowr site; the published site (`analysis/run_susier.Rmd`) renders these tables directly and is the canonical output.

**2. Joint mvSuSiE-RSS** (all trait-pair combinations)

- `07_mvsusie_dmfs_nteeth.R` — DMFS + Nteeth
- `11_mvsusie_dmfs_dfss.R` — DMFS + DFSS
- `12_mvsusie_nteeth_dfss.R` — Nteeth + DFSS
- `13_mvsusie_dmfs_nteeth_dfss.R` — DMFS + Nteeth + DFSS

Each of these depends on `01` plus the single-trait outputs (steps in the table above) for the traits it combines.

**3. Convergence pass**

- `14_batch_refit_nonconverged.R` — re-fits, with a higher iteration budget (`max_iter = 3000` vs. `1000`), any locus from step 2 that did not converge on the first pass. Patches the results/summary files from step 2 in place. Safely re-runnable: it automatically skips loci that already show `converged = TRUE`.

### Key parameters

- `L = 10` single effects, prior variance matrix `matrix(c(0.2, 0.1, 0.1, 0.2), nrow = 2)` for all joint analyses
- `max_iter = 1000`, `tol = 1e-04` for the initial pass (steps 1-2); `max_iter = 3000` for the convergence pass (step 3)
- Per-locus LD matrices are checked for positive semi-definiteness and regularized (shrunk toward the identity matrix by the minimum amount necessary) before every fit — real reference-panel LD blocks are not always PSD, and fitting a non-PSD matrix causes the model's ELBO to diverge instead of converge. See any of the scripts above for the exact regularization code.
- `set.seed(20260523)` is used wherever randomness is involved.

Model convergence (`fit$converged`) is verified for every locus in every analysis; see `analysis/run_susier.Rmd` for the current convergence status and the reasoning behind the two-pass (steps 2 + 3) approach.

## Input files

Raw GWAS summary statistics should be placed in:

```text
data_samples/
```

Locus-specific LD matrices (`.ld`, `.snplist`, `.id_map` files per locus) are expected under:

```text
LDMatrix/
```

Both directories are gitignored (too large for version control) — regenerate or re-download them before running the pipeline from scratch.

## Reproducibility

- R session info (R version, package versions) and the exact git commit each page was built from are recorded automatically at the bottom of every page on the published site, via workflowr.
- Intermediate results (`results/*.RDS`, most `results/*.csv`) are gitignored, since they're large and fully regenerable from the pipeline above — only the final published HTML and a handful of small summary CSVs are tracked in git.
- To rebuild a single page locally: `workflowr::wflow_build("analysis/run_susier.Rmd")` (requires pandoc; if not on `PATH`, set the `RSTUDIO_PANDOC` environment variable to a pandoc binary, e.g. one bundled with RStudio).
