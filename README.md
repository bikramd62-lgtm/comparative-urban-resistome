# Comparative Urban Resistome

This repository contains code and processed analysis outputs for a Master's thesis comparing Antimicrobial Resistance patterns in urban sewage and public-transit metagenomes across matched global cities.

## Project overview and Research questions

This thesis investigates urban antimicrobial resistance patterns across two contrasting but potentially connected urban environments: wastewater systems and public-transit surfaces. The analysis is structured around three research questions:
- RQ1. Do cities share a common urban resistome across sewage and public-transit environments?
This question examines the overlap in antimicrobial resistance gene profiles between matched sewage and transit datasets across the final city panel.
- RQ2. Which ARG features are more specific to sewage or public-transit environments?
This question identifies environment-specific and environment-associated ARG groups and resistance classes using prevalence-based comparisons.
- RQ3. Are city-level characteristics associated with sewage–transit resistome overlap?
This question explores whether climatic, demographic, geographic, socioeconomic, and sanitation-related city-level variables are associated with ARG overlap between sewage and transit environments.

## Workflow overview

The workflow combines bioinformatic preprocessing, ARG profiling, presence/absence matrix construction, statistical analysis, and figure generation.

The main workflow steps are:

Metadata parsing and matched-city selection
Sewage and transit metadata were parsed to identify the final matched city panel used for cross-environment comparison.
Read quality control and trimming
Raw sequencing reads were quality-checked and trimmed using standard read preprocessing workflows.
Human-read depletion
Reads were screened against a human reference genome to reduce host-associated sequence content before downstream ARG profiling.
City-wise read pooling
Reads were pooled at the city level within each environment to create one sewage and one transit profile per matched city.
Subsampling
City-level read pools were standardized to a fixed sequencing depth to reduce technical imbalance between samples.
ARG mapping and profiling
Reads were mapped against MEGARes v3.00, and ARG profiles were processed using an 80% gene-fraction criterion.
Presence/absence matrix construction
ARG count tables were converted into binary presence/absence matrices. A 10% prevalence filter was applied to retain recurrent features across city-environment profiles.
RQ1 shared resistome analysis
Jaccard-based overlap, ordination, PERMANOVA, within-vs-between city similarity, and shared feature analyses were performed.
RQ2 environment-specificity analysis
Strict environment-exclusive features and prevalence-based sewage/transit-associated features were identified and visualized.
RQ3 urban covariate analysis
City-level climatic, demographic, geographic, socioeconomic, and sanitation-related variables were compared against sewage–transit ARG overlap metrics.
Sensitivity and robustness analysis
Additional analyses tested the influence of sequencing depth, ARG count imbalance, richness effects, sample number, and turnover/nestedness structure.

## Repository structure

- `scripts/`: Bioinformatics and downstream analysis scripts
- `scripts/01_metadata/`: Metadata parsing script for matched city selection across urban sewage and public-transit systems
- `scripts/02_download/`: Script for downloading selected urban sewage raw reads across matched cities with RUN IDs after sample selection, assigned to each sample of a specific city.
- `scripts/Sewage/`: Script for bioinformatic preprocessing steps for the retrieved sewage sample raw reads, including read-quality assessment and trimming, human-read removal, city-wise pooling, subsampling to 5 million reads, ARG mapping, and prevalence filtering using 10% filtering criteria.
- `scripts/Transit/`: Script for bioinformatic preprocessing steps for the retrieved transit sample raw reads, including read-quality assessment and trimming, human-read removal, city-wise pooling, subsampling to 5 million reads, ARG mapping, and prevalence filtering using 10% filtering criteria.
- `scripts/appendix_figures/`: Script for all the figures provided in the appendix section of the thesis report.
- `scripts/comparative_resistome_RQ1/`: Scripts for the whole RQ1 analysis used for the thesis with relative paths.
- `scripts/comparative_resistome_RQ1/06_RQ1_figures_used_in_thesis/`: Scripts for the figures of RQ1 anaylsis used in the thesis report with relative paths.
- `scripts/environment_specific_RQ2/`: Scripts for the whole RQ2 analysis used for the thesis with relative paths.
- `scripts/environment_specific_RQ2/03_RQ2_figures_used_in_thesis/`: Scripts for the figures of RQ2 analysis used in the thesis report with relative paths.
- `scripts/urban_covariate_analysis_RQ3/`: Scripts for the whole RQ3 analysis used for the thesis with relative paths.
- `scripts/urban_covariate_analysis_RQ3/03_RQ3_figures_used_in_thesis/`: Scripts for the figures of RQ3 analysis used in the thesis report with relative paths.
- `scripts/sensitivity_analysis/`: Scripts for the whole sensitivity analysis used for the thesis with relative paths.
- `scripts/sensitivity_analysis/sensitivity_analysis_figures_used_in_the_thesis/`: Scripts for the figures of sensitivity analysis used in the thesis report with relative paths.
- `data/metadata/`: Metadata and city-level covariate files.
- `data/processed/`: Processed matrices and derived analysis tables.
- `results/`: Statistical outputs and summary tables.
- `results/RQ1_figures/`: Figures generated during the RQ1 analysis.
- `results/RQ1_output_tables/`: Tables generated during the RQ1 analysis.
- `results/RQ2_figures/`: Figures generated during the RQ2 analysis.
- `results/RQ2_output_tables/`: Tables generated during the RQ2 analysis.
- `results/RQ3_figures/`: Figures generated during the RQ3 analysis.
- `results/RQ3_output_tables/`: Tables generated during the RQ3 analysis.
- `results/Sensitivity_analysis_figures/`: Figures generated during the sensitivity analysis.
- `results/Sensitivity_analysis_output_tables/`: Tables generated during the sensitivity analysis.
- `docs/`: whole Workflow documentation

## Data and privacy note

Raw sequencing reads are not included in this repository due to file size and privacy considerations. The repository contains scripts, metadata, processed matrices, and derived outputs required to inspect and reproduce the downstream analysis.

## Software requirements

All the tools, software, and databases used in this thesis are listed below with their corresponding versions wherever applicable.

- Ubuntu/Linux - 24.04.3 LTS, kernel 6.14.0-37-generic
- micromamba - 2.6.2-1
- GNU Wget - 1.21.4
- Unix cat - 9.4
- FastQC - 0.12.1
- MultiQC - 1.32
- fastp - 0.23.4
- KneadData - 0.12.4
- Bowtie2 - 2.5.4
- SAMtools - 1.22.1
- seqtk - 1.5-r133
- hg39 Bowtie2 index - Index files dated 17 June 2021
- MEGARes - 3.00
- AMR++/ResistomeAnalyzer - Not applicable
- R - 4.5.2
- readr - 2.1.5
- dplyr - 1.1.4
- stringr - 1.6.0
- tidyr - 1.3.1
- ggplot2 - 4.0.0
- vegan - 2.7.2
- ape - 5.8.1
- pheatmap - 1.0.13
- ggrepel - 0.9.6
- tibble - 3.3.0
- purrr - 1.2.0

## Analysis reproducibility


## Main outputs


## Contact

