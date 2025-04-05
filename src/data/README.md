# NASH-HCC vs Normal Liver miRNA Expression Analysis

This directory contains scripts to acquire and analyze miRNA expression data for comparing NASH-HCC (Non-alcoholic steatohepatitis-associated hepatocellular carcinoma) with Normal Liver tissue. The goal is to identify differentially expressed miRNAs that can serve as potential therapeutic targets.

## Data Acquisition Options

We provide multiple scripts to acquire data from different sources:

### 1. UCSC Xena Browser (TCGA Data)

The `fetch_xena_data.py` script downloads miRNA expression data from the TCGA Liver Hepatocellular Carcinoma (LIHC) cohort via the UCSC Xena Browser.

```bash
python src/data/fetch_xena_data.py --output-dir hackathonbro/src/data/processed
```

This script will:
- Download miRNA expression data from TCGA LIHC dataset
- Download clinical metadata
- Process and merge the data to create a format suitable for differential expression analysis
- Automatically classify samples as "Normal" or "NASH-HCC" based on sample type in clinical data

### 2. Gene Expression Omnibus (GEO)

The `fetch_geo_data.py` script downloads miRNA expression data from GEO for studies specifically focused on NASH-HCC or related liver conditions.

```bash
python src/data/fetch_geo_data.py --geo-id GSE123456 --output-dir hackathonbro/src/data/processed
```

This script will:
- Download the specified GEO dataset
- Extract sample metadata and expression data
- Classify samples as "Normal" or "NASH-HCC" based on keyword matching in sample descriptions
- Process the data into a format suitable for differential expression analysis

## Analysis Options

We provide two analysis pipelines:

### 1. Python-based Statistical Analysis (t-test/Wilcoxon)

The existing `target_prioritizer.py` script performs differential expression analysis using t-tests or Wilcoxon rank-sum tests.

```bash
python src/target_prioritizer.py --input hackathonbro/src/data/processed/TCGA_LIHC_miRNA_processed.csv --id-col miRNA
```

This is suitable for normalized data (e.g., log2(RPM+1) from Xena Browser).

### 2. R-based DESeq2 Analysis

For count-based data and more sophisticated differential expression analysis, we provide the `deseq2_analysis.R` script.

```bash
Rscript src/data/deseq2_analysis.R --input hackathonbro/src/data/processed/TCGA_LIHC_counts.csv --output hackathonbro/results/deseq2_results.csv
```

DESeq2 is better suited for raw count data and provides more robust statistical modeling for RNA-seq data.

## Installation

### Python Dependencies

Install Python dependencies with:

```bash
pip install -r requirements.txt
```

### R Dependencies (for DESeq2 analysis)

Install R dependencies within R:

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "dplyr", "readr", "argparse", "ggplot2", "pheatmap"))
```

## Workflow Example

Here's a typical workflow:

1. Download data from UCSC Xena Browser:
   ```bash
   python src/data/fetch_xena_data.py --output-dir hackathonbro/src/data/processed
   ```

2. Run differential expression analysis:
   - For normalized data (e.g., log2(RPM+1)):
     ```bash
     python src/target_prioritizer.py --input hackathonbro/src/data/processed/TCGA_LIHC_miRNA_processed.csv --id-col miRNA
     ```
   - For raw count data:
     ```bash
     Rscript src/data/deseq2_analysis.R --input hackathonbro/src/data/processed/TCGA_LIHC_counts.csv --output hackathonbro/results/deseq2_results.csv
     ```

3. Review results:
   - Differential expression results are saved to CSV files
   - Visualization plots are generated (volcano plots, heatmaps, etc.)

## Notes

- For TCGA data, normal samples are labeled as "Solid Tissue Normal" in clinical data
- NASH-HCC samples are a subset of HCC samples; additional filtering based on clinical variables (e.g., steatosis, cirrhosis) may be required to specifically identify NASH-HCC
- DESeq2 analysis requires raw count data; if only normalized data is available, use the Python-based analysis 