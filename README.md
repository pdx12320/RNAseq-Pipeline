# RNA-seq Analysis Pipeline

A professional, end-to-end bioinformatics pipeline for bulk RNA-seq data analysis. This project covers the entire workflow from raw FASTQ files to advanced systems biology analysis (WGCNA).

## Description
This repository provides a standardized framework for:
- Raw data quality control and trimming.
- Efficient read alignment and post-processing.
- Gene-level quantification.
- Differential expression analysis (DEA) with high-quality visualizations.
- Functional interpretation via GO, KEGG, and GSEA.
- Gene co-expression network construction (WGCNA).

## Repository Structure
- `scripts/`: Contains all modular Bash and R scripts.
- `README.md`: Documentation and usage guide.
- `LICENSE`: MIT License.

---

## Workflow Overview

### 1. Upstream Processing (Shell Scripts)
All upstream steps are automated in shell scripts located in `/scripts`.

* [cite_start]**Environment Setup**: Configuration using Conda and SRA-Toolkit[cite: 5, 9].
* [cite_start]**Preprocessing**: Quality control with FastQC and adapter trimming with Cutadapt[cite: 18, 23]. 
    * *Script*: `scripts/01_preprocessing.sh`
* [cite_start]**Alignment**: Supports Hisat2 and STAR for mapping reads to a reference genome[cite: 60, 62].
    * *Script*: `scripts/02_alignment.sh`
* [cite_start]**Quantification**: Counting reads per gene using featureCounts[cite: 112, 113].
    * *Script*: `scripts/03_quantification.sh`

### 2. Downstream Analysis (R Scripts)
Downstream statistical analysis is performed using R (v4.0+).

* [cite_start]**Differential Expression**: Utilizing `DESeq2` for normalization and identifying DEGs[cite: 126, 151].
    * *Script*: `scripts/04_deseq2_analysis.R`
* [cite_start]**Functional Enrichment**: GO/KEGG Over-representation Analysis (ORA) and Gene Set Enrichment Analysis (GSEA)[cite: 239, 330, 365].
    * *Script*: `scripts/05_enrichment.R`
* [cite_start]**Network Analysis**: WGCNA for identifying highly correlated gene modules[cite: 445, 469].
    * *Script*: `scripts/06_wgcna.R`

---

## Quick Start

### Prerequisites
Ensure you have `Conda` installed. Create the environment:
```bash
conda create -n rnaseq_env bioconda::sra-tools bioconda::fastqc bioconda::cutadapt bioconda::hisat2 bioconda::subread
conda activate rnaseq_env
