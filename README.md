# rbd-epistasis-landscape

Analysis workflows and processed data for mapping ACE2-binding enrichment landscapes in the SARS-CoV-2 receptor-binding domain (RBD).

This repository accompanies the doctoral dissertation:

**Jonathan Felix. _Experimental Mapping of the ACE2-Binding Fitness Landscape in the SARS-CoV-2 Receptor-Binding Domain._**  
Doctoral dissertation, Keck Graduate Institute.  
Available through ProQuest: https://www.proquest.com/docview/3333498950

## Project overview

This project examines how combinations of mutations in the SARS-CoV-2 spike receptor-binding domain affect ACE2-associated enrichment in a yeast display selection system. The analysis is based on a combinatorial RBD variant library spanning a binary Wuhan/Omicron sequence design, enabling genotype-level comparison across a defined mutational landscape.

Sequencing-derived genotype counts are used to calculate normalized frequencies, enrichment values, statistical comparisons, mutational-load trends, and landscape visualizations. The repository is intended to support reproducibility of the dissertation-associated analyses and provide a transparent record of the computational workflow used to interpret the experimental selection data.

## Repository purpose

This repository is maintained as an archival analysis repository linked to the dissertation record. It contains code, processed analysis tables, and figure-generation workflows associated with the dissertation.

The repository is not currently intended as a general-purpose software package. Instead, it is organized to document and reproduce the specific analyses performed for this SARS-CoV-2 RBD yeast display landscape study.

## Analysis scope

The workflows in this repository support analyses including:

- parsing and organizing genotype-level sequencing count data
- calculating normalized genotype frequencies across library conditions
- estimating enrichment values for ACE2-associated positive selection and control conditions
- comparing enrichment patterns across replicate libraries
- evaluating mutational-load trends across the binary RBD landscape
- visualizing genotype-level enrichment and epistatic structure
- performing sensitivity analyses, including low-read-count filtering checks

## Repository status

This repository is linked in the dissertation record and is being organized for long-term reproducibility. Code and documentation may be cleaned, clarified, and reorganized after dissertation publication, but the repository is intended to preserve the analysis logic associated with the dissertation.

Additional manuscript-specific analyses, revised quality-control filters, or publication-ready figure workflows may be added in future releases or maintained separately.

## Data availability

Processed genotype count tables and analysis-ready files required to reproduce the dissertation-associated analyses will be included in this repository when available.

Raw sequencing files are not hosted directly in this GitHub repository due to file size and archival considerations. Raw sequencing data will be deposited in an appropriate public archive or made available upon reasonable request.

## Citation

If using this repository or associated analyses, please cite the dissertation:

Felix, Jonathan. _Experimental Mapping of the ACE2-Binding Fitness Landscape in the SARS-CoV-2 Receptor-Binding Domain._ Doctoral dissertation, Keck Graduate Institute. ProQuest, 2026. https://www.proquest.com/docview/3333498950
