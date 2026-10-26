# (MPXV) - Côte d'Ivoire (2025)

This repository contains the scriptsR under a *One Health* approach to analyze Mpox virus (MPXV) strains circulating in Côte d'Ivoire in 2025.

## 🔬 Objectives of the Analysis

The project focuses on two main components:
1. **Phylogenetic Reconstruction and Visualization:** Generation of robust phylogenetic trees (`Phylogenetique_Tree_Mpox.R`) integrating epidemiological metadata (`annotations.txt`, `Mpox_Tree_file.treefile`).
2. **Transmission Chain Characterization:** Identification and differentiation of human-to-human transmission events versus zoonotic spillover introductions using APOBEC3-mediated mutation profile analysis (`APOBEC3.R`, `apobec3_summary_Mpox2025.csv`).

## Repository Structure

- `Phylogenetique_Tree_Mpox.R` : R script for tree formatting and annotated figure generation.
- `APOBEC3.R` : Script for mutation profile analysis and transmission route identification.
- `Mpox_Tree_file.treefile` : Input phylogenetic tree file.
- `annotations.txt` : Epidemiological metadata table.
- `Mpox_Date_2025.csv` & `apobec3_summary_Mpox2025.csv` : Processed data tables.

## Prerequisites and Dependencies (R)

The scripts require the following R packages:

```R
install.packages(c("ape", "phytools", "ggplot2", "dplyr", "ggnewscale"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ggtree")