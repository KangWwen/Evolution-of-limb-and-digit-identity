Evolution of limb and digit identity genes since the tetrapod ancestor

## Overview

This repository contains the analysis code and figure scripts associated 
with the manuscript:
Kang, W., Fei, JF., Liang, C. et al. Evolution of limb and digit identity genes since the tetrapod ancestor. Nat Commun (2026). https://doi.org/10.1038/s41467-026-73821-7

Comparative transcriptomics of digits from 6 species was used to identify 
gene expression signatures associated with the evolution of digit number 
and morphology, and to characterize morphological fingerprint genes (MFGs).

---

## Data Availability

| Dataset | Source | Accession |
|---|---|---|
| axolotl, turtle, crocodile, chicken, emu, ostrich (this study) | NCBI SRA | PRJNA1258050 |
| mouse and Anolis digits (published) | NCBI GEO | GSE108337 |


## Repository Structure
```
├── Alignment.sh                      # Read alignment pipeline (raw FASTQ → BAM/counts)
├── CellRangerPipelines.sh            # Cell Ranger pipeline for scRNA-seq data
├── DEGs.r                            # Differential expression analysis (DESeq2/edgeR)
├── TauValueCaculationFunction.r      # Tau tissue-specificity index calculation
├── PhyDGET/                          # Branch-specific forelimb and hindlimb expression changes
└── FigureScripts/
├──── Fig1b Clusters Between Species.R
├──── Fig1c Overlapped Number of MFGs.R
├──── Fig1d Heatmap of MADness Genes.R
├──── Fig2c&2d Expression Heatmap for DEGs.R
├──── Fig3a Number of D1 MFGs in Forelimb and Hindlimb.R
├──── Fig3a Piechart of MFGs number.R
├──── Fig3b Density of Vertebrate Genes.R
├──── Fig3c Filtered GO Enrichment of MFGs.R
├──── Fig3d&3e scRNA-seq Plot Markers.R
└── Supplementary Figures/
```
