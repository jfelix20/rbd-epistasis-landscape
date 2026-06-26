# rbd-epistasis-landscape

Reproducible analysis workflows and processed data for mapping ACE2-associated enrichment landscapes in the SARS-CoV-2 receptor-binding domain (RBD).

This repository accompanies the doctoral dissertation:

**Jonathan Felix. *Experimental Mapping of the ACE2-Binding Fitness Landscape in the SARS-CoV-2 Receptor-Binding Domain.***
Doctoral dissertation, Keck Graduate Institute.
Available through ProQuest: https://www.proquest.com/docview/3333498950

## Project overview

This project examines how combinations of mutations in the SARS-CoV-2 spike receptor-binding domain affect ACE2-associated enrichment in a yeast display selection system. The analysis is based on a defined binary Wuhan/Omicron combinatorial RBD variant library, enabling genotype-level comparison across a complete 8-site mutational landscape.

Sequencing-derived genotype counts are used to calculate normalized frequencies, enrichment values, statistical comparisons, mutational-load trends, path/network summaries, epistasis decompositions, and dissertation-associated landscape visualizations.

The repository is intended to support reproducibility of the dissertation-associated analyses and provide a transparent record of the computational workflow used to interpret the experimental selection data.

## Repository purpose

This is an archival dissertation reproducibility repository. It contains processed data, analysis notebooks, checkpoint tables, and figure-generation workflows associated with the dissertation.

The repository is not intended as a general-purpose software package. It is organized as a numbered notebook workflow that reconstructs the dissertation analyses from processed genotype-level count and landscape tables.

## Analysis scope

The workflows in this repository support analyses including:

* project setup and data inventory
* recalculation of genotype frequencies and enrichment metrics
* low-read sample sensitivity analysis
* regeneration of dissertation-associated landscape figures and summary tables
* replicate and genotype-frequency correlation analyses
* network and path analysis across the binary genotype landscape
* pairwise and third-order epistasis calculations
* Walsh-style variance decomposition by interaction order
* residual landscape reconstruction from archived human-model comparison files
* final repository-level reproducibility inventory

## Notebook workflow

The notebooks are intended to be read and run in numerical order:

| Notebook                                                 | Purpose                                                                             |
| -------------------------------------------------------- | ----------------------------------------------------------------------------------- |
| `00_project_setup_and_data_inventory.ipynb`              | Project setup, path checks, and data inventory                                      |
| `01_recalculate_landscape_tables.ipynb`                  | Recalculate pooled frequencies, enrichment metrics, and landscape tables            |
| `02_low_read_filter_sensitivity.ipynb`                   | Evaluate sample-level read-depth cutoff sensitivity                                 |
| `03_dissertation_figures_from_landscape_tables.ipynb`    | Regenerate dissertation figures and summary tables from processed landscape outputs |
| `04_network_and_path_analysis.ipynb`                     | Reproduce network/path analyses and weighted optimal path calculations              |
| `05_epistasis_analysis.ipynb`                            | Reproduce pairwise, third-order, and variance-decomposition epistasis analyses      |
| `06_chris_model_comparison_and_residual_landscape.ipynb` | Reconstruct archived human-model comparison and residual landscape workflow         |
| `07_repository_index_and_reproducibility_summary.ipynb`  | Final repository index and reproducibility summary                                  |

Notebook `07` performs no new scientific analysis. It provides a repository-level audit and public-release index.

## Repository structure

```text
archive/      Archived legacy scripts and source files used for reconstruction
config/       Analysis configuration files
data/         Processed and analysis-ready data tables
docs/         Supporting documentation
figures/      Generated analysis figures
notebooks/    Numbered reproducibility notebooks
outputs/      Intermediate checkpoint tables and inventories
src/          Reusable project code, if applicable
```

## Data availability

Processed genotype count tables and analysis-ready files required to reproduce the dissertation-associated analyses are included in this repository.

Raw sequencing files are not hosted directly in this GitHub repository due to file size and archival considerations. Raw sequencing data will be deposited in an appropriate public archive or made available upon reasonable request.

## Environment

The repository includes both Conda and pip environment specifications.

To create the Conda environment:

```bash
conda env create -f environment.yml
conda activate rbd-epistasis
```

To install dependencies with pip instead:

```bash
pip install -r requirements.txt
```

The notebooks were reconstructed using the project environment described in `environment.yml` and `requirements.txt`.

## Reproducibility notes

The numbered notebooks write intermediate tables to `outputs/checkpoints/` and figures to `figures/`. Several notebooks also produce inventory tables documenting their generated outputs.

Some archived files are retained to document the original dissertation-analysis context. The final numbered notebooks are intended to stand alone and should be treated as the primary reproducibility workflow.

## Citation

If using this repository or associated analyses, please cite the dissertation:

Felix, Jonathan. *Experimental Mapping of the ACE2-Binding Fitness Landscape in the SARS-CoV-2 Receptor-Binding Domain.* Doctoral dissertation, Keck Graduate Institute. ProQuest, 2026. https://www.proquest.com/docview/3333498950
