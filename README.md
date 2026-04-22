Evolution of limb and digit identity genes since the tetrapod ancestor

## Overview

This repository contains the analysis code and figure scripts associated 
with the manuscript:
"Evolution of limb and digit identity genes since the tetrapod ancestor"

Comparative transcriptomics of digits from 6 species was used to identify 
gene expression signatures associated with the evolution of digit number 
and morphology, and to characterize morphological fingerprint genes (MFGs).

---

## Data Availability

| Dataset | Source | Accession |
|---|---|---|
| axolotl, turtle, crocodile, chicken, emu, ostrich (this study) | NCBI SRA | PRJNA1258050 |
| mouse and Anolis digits (published) | NCBI GEO | GSE108337 |

Processed input files used directly by figure scripts are provided in 
the `/InputFiles` folder.

## Repository Structure
├── Alignment.sh                      # Read alignment pipeline (raw FASTQ → BAM/counts)
├── CellRangerPipelines.sh            # Cell Ranger pipeline for scRNA-seq data
├── DEGs.r                            # Differential expression analysis (DESeq2/edgeR)
├── TauValueCaculationFunction.r      # Tau tissue-specificity index calculation
├── PhyDGET/                          # Branch-specific forelimb and hindlimb expression changes
├── InputFiles/
│   ├── TPM_withSqrt/                 # Square-root-transformed TPM matrices; used for
│   │                                 #   clustering analysis (Fig1b, FigS2, FigS5)
│   ├── Average_TPM/                  # Per-species average TPM matrices; used for
│   │                                 #   expression heatmaps (Fig1d, Fig2c&2d,
│   │                                 #   FigS12&14&15)
│   └── Tau_inEachSpecies/            # Tau scores and MFG lists per species; used for
│                                     #   all MFG-related analyses (Fig1c, Fig3a–3e,
│                                     #   FigS3, FigS6, FigS8, FigS16, FigS25)
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
---

## How to Run

###  Raw data processing
raw sequencing reads would be found through Data Availability Accession
bash Alignment.sh             # Align bulk RNA-seq reads; generates count matrices

##	Pre-computed processed files are already available in `/InputFiles` 
Differential expression and Tau calculation
DEGs.r
TauValueCaculationFunction.r
Output files from these scripts are pre-computed and provided in 
`/InputFiles/Tau_inEachSpecies`.

##	Generate figures
Each script in `/FigureScripts` reads input files directly from 
`/InputFiles`. Run scripts individually in R to reproduce the 
corresponding figure panels.