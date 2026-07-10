# Comparative Urban Resistome

This repository contains code and processed analysis outputs for a Master's thesis comparing Antimicrobial Resistance gene profiles in urban sewage and public-transit metagenomes across matched global cities.

## Project overview

The thesis investigates whether urban sewage and public-transit environments share a common resistome across matched cities, which antimicrobial resistance features distinguish these environments, and whether city-level characteristics are associated with sewage-transit resistome overlap.

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

## Data availability

Raw sequencing reads are not included in this repository due to file size and privacy considerations. The repository contains scripts, metadata, processed matrices, and derived outputs required to inspect and reproduce the downstream analysis.

## Software versions

All the tools, software, and databases used in this thesis are listed below with their corresponding versions wherever applicable.

Ubuntu/Linux - 24.04.3 LTS, kernel 6.14.0-37-generic
micromamba - 2.6.2-1
GNU Wget - 1.21.4
Unix cat - 9.4
FastQC - 0.12.1
MultiQC - 1.32
fastp - 0.23.4
KneadData - 0.12.4
Bowtie2 - 2.5.4
SAMtools - 1.22.1
seqtk - 1.5-r133
hg39 Bowtie2 index - Index files dated 17 June 2021
MEGARes - 3.00
AMR++/ResistomeAnalyzer - Not applicable
R - 4.5.2
readr - 2.1.5
dplyr - 1.1.4
stringr - 1.6.0
tidyr - 1.3.1
ggplot2 - 4.0.0
vegan - 2.7.2
ape - 5.8.1
pheatmap - 1.0.13
ggrepel - 0.9.6
tibble - 3.3.0
purrr - 1.2.0
