# Fecal Microbiome and Metabolome Signatures of Severe and Attenuated Familial Adenomatous Polyposis
This repository contains the analysis code for the paper **"Fecal Microbiome and Metabolome Signatures ofSevere and Attenuated Familial Adenomatous Polyposis"**.

## Overview

Familial Adenomatous Polyposis (FAP) is an early-onset colorectal cancer syndrome. This study integrates gut microbiome profiles (16S rRNA), fecal metabolomics (LC-MS), and clinical blood biochemistry to identify robust multi-omic signatures that track with disease severity (Classical vs. Attenuated phenotypes).

The pipeline performs quality control, batch correction, taxonomic/metabolomic profiling, and integrated multi-omic analysis combining taxonomy, metabolomics and blood biochemistry test results to predict and characterize FAP.

## Repository Structure

The analysis is organized into a sequential series of R Markdown notebooks and helper scripts:

| Directory/File | Description |
| :--- | :--- |
| `python/metabolomics_data_preprocessing.ipynb` | **Pre-processing:** Preprocessing and ID mapping of raw metabolomics data. |
| `R/00_process_metadata.Rmd` | **Cohort Definition:** Cleaning clinical metadata and defining study groups (FAP vs. Control). |
| `R/01_QC.Rmd` | **Microbiome QC:** Rarefaction curves, prevalence filtering, and sequencing depth evaluation. |
| `R/02_compare_batches.Rmd` | **Batch Correction:** Evaluation and correction of batch effects (Batch 22 vs. 24). |
| `R/03_taxonomic_analysis.Rmd` | **Microbiome:** Figure 1: Alpha/Beta diversity analysis and differential abundance testing. |
| `R/04_metabolomic_analysis.Rmd` | **Metabolome:** Figure 2: PCA, univariate enrichment analysis, and KEGG pathway over-representation analysis (ORA). |
| `R/05_multiomic_analysis.Rmd` | **Integration:** Figure 3: Logistic regression classification (AUC analysis) and multi-omic module discovery using MintTea. |
| `R/statistics.R` | Helper functions for statistical tests (Wilcoxon, Permutation ANCOVA). |
| `R/visualization.R` | Custom `ggplot2` themes and plotting wrappers. |
| `R/utils.R` | Small helper functions.  |
