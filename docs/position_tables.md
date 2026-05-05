# Position tables

The original position-level sequencing tables are upstream intermediate files generated from the deduplicated sequencing outputs. These files contain position-level amino acid calls and UMI counts for each sequencing library.

The full position tables are not tracked directly in this GitHub repository because several files are very large and are better suited for deposition in a dedicated data archive.

The processed W/O genotype count table used for the main dissertation-associated analysis is included at:

`data/processed/wo_counts.csv`

The expected position table files are listed in:

`data/metadata/position_table_manifest.csv`

## Sequencing library coverage notes

The final position table inventory contains 24 sequencing library files. Three expected sequencing libraries are absent from the position table set because they returned zero reads:

- `Library-2-Pos-1`
- `Library-3-Neg-1`
- `Library-3-Pre-1`

Two additional sequencing libraries returned only one read each. These very low-read samples are retained in the archived upstream data inventory but will be evaluated during cleanup and sensitivity analysis to determine whether they were included in the original processed workflow and whether excluding them affects the downstream landscape calculations.

A specific goal of the cleanup analysis is to verify whether the original pipeline excluded zero-read and near-zero-read samples before pooling counts by library and condition. If near-zero-read samples were included, the cleaned workflow will document where they enter the analysis and evaluate whether removing them changes the reported enrichment patterns.

## Data availability

Full upstream position tables will be deposited in an appropriate external archive or made available upon reasonable request. GitHub is used here for code, processed analysis tables, metadata, and documentation rather than large sequencing-derived intermediate files.
