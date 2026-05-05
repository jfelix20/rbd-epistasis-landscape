# Processed data tables

This folder contains processed W/O genotype count and landscape analysis tables used in the dissertation-associated SARS-CoV-2 RBD yeast display analysis.

These files are derived from upstream position-level sequencing tables and are small enough to be tracked directly in GitHub. Large upstream position tables and raw sequencing files are not stored in this repository.

## Files

### `wo_counts.csv`

Replicate-level W/O genotype count table.

This is the most upstream processed table currently tracked in the repository. It contains W/O genotype calls and UMI counts for each included sequencing sample before pooling across replicates.

Expected columns include:

- `wo`
- `umi_count`
- `library_id`
- `condition`
- `replicate`
- `sample`

This table is the preferred starting point for cleaned recalculation workflows and sample-level read-depth auditing.

### `WO_genotype_counts_per_library_condition.csv`

Pooled W/O genotype count table.

This file is derived from `wo_counts.csv` by pooling replicate counts by library, condition, and W/O genotype. It contains the pooled UMI count, total condition count, and normalized frequency for each genotype within each library-condition group.

Expected columns include:

- `library_id`
- `condition`
- `wo`
- `umi_count`
- `total_cond`
- `freq`

### `WO_landscape_per_library.csv`

Primary W/O landscape table.

This file contains genotype-level normalized frequencies and enrichment values for each library. It includes positive-selection, negative-selection, and pre-selection frequency values, along with calculated enrichment metrics such as Pos/Pre, Neg/Pre, mutational load, and delta-delta-style comparisons.

This table is a primary source table for dissertation-associated landscape analyses and figures.

### `WO_landscape_with_pvalues.csv`

Landscape table with statistical comparisons.

This file extends the primary W/O landscape table with count columns and statistical comparisons, including p-values and multiple-testing-adjusted q-values for relevant condition contrasts.

This table is useful for reproducing volcano plots and significance-based summaries.

### `WO_landscape_ddelta_posneg_volcano_source.csv`

Volcano-plot source table.

This file is a downstream figure-source table focused on Pos-vs-Neg delta-delta-style enrichment comparisons and volcano plot categorization.

It is included as a checkpoint from the final analysis run and as a source for reproducing dissertation-associated volcano visualizations.

### `amino_acid_frequency_matrix_prefilter_dissertation_heatmap.csv`

Processed amino-acid frequency matrix used to reproduce the dissertation-facing amino-acid replicate correlation heatmap.

Rows correspond to sequencing samples, and columns correspond to site/amino-acid frequency features. This matrix reflects the pre-filtered amino-acid frequency representation used for the original dissertation figure.

### `amino_acid_frequency_matrix_design_sites_filtered_heatmap.csv`

Filtered/design-site amino-acid frequency matrix retained for comparison.

Rows correspond to sequencing samples, and columns correspond to amino-acid frequency features across the designed mutation sites. This matrix produces a different replicate correlation structure from the dissertation-facing pre-filtered matrix and is retained to document the effect of filtering/feature selection on the amino-acid frequency heatmap.

## Processing hierarchy

The tracked processed tables follow this approximate analysis hierarchy:

wo_counts.csv
    -> pool replicates by library + condition + genotype
WO_genotype_counts_per_library_condition.csv
    -> normalize frequencies and calculate enrichment values
WO_landscape_per_library.csv
    -> add statistical comparisons and multiple-testing correction
WO_landscape_with_pvalues.csv
    -> format Pos-vs-Neg delta-delta volcano source table
WO_landscape_ddelta_posneg_volcano_source.csv

## Notes

The cleaned analysis workflow should start from `wo_counts.csv` when possible, because it preserves replicate-level sample information needed for read-depth auditing and low-read-count sensitivity analysis.

The downstream landscape tables are retained as checkpoints from the final dissertation-associated analysis run and will be used to verify that cleaned notebooks reproduce the original calculations.
