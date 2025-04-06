# AI-Powered miRNA Drug Discovery Platform

A comprehensive computational platform for identifying, prioritizing, and evaluating miRNA-based therapeutics, demonstrated for Hepatocellular Carcinoma (HCC).

## Project Overview

This platform addresses the critical challenge of targeted drug delivery and efficient drug discovery in the rapidly growing field of RNA therapeutics. By leveraging disease-specific microRNA (miRNA) expression patterns and employing advanced AI and computational biology techniques, we identify and optimize potential miRNA drug candidates.

Our system integrates multi-omics data (genomics, transcriptomics, protein interactions) to:
1.  Identify dysregulated miRNAs in disease states (initially focused on NASH-HCC vs. Normal Liver).
2.  Predict and prioritize miRNA target genes based on network centrality and biological relevance.
3.  Assess the druggability of identified targets using resources like Open Targets.
4.  Analyze miRNA binding site accessibility using RNA structure prediction (ViennaRNA).
5.  Perform *in silico* safety checks for potential off-target effects.
6.  Generate comprehensive reports and interactive visualizations to guide therapeutic design.

## Hackathon Focus

Due to time constraints, our hackathon implementation focuses on one of two core components:

### Option A: Target Prioritizer
Identifies and ranks miRNAs differentially expressed between NASH-HCC and Normal Liver tissue. This helps select the most promising miRNA targets for therapeutic design.

#### New Features Added:
- Automated data acquisition from multiple sources (TCGA/Xena Browser, GEO)
- Flexible data processing pipelines
- Statistical analysis using both Python (t-test/Wilcoxon) and R (DESeq2)
- End-to-end workflow for differential expression analysis

### Option B: Accessibility Predictor
Optimizes the placement of miRNA binding sites (particularly miR-122 for liver targeting) within mRNA constructs by analyzing RNA secondary structure and accessibility.

## Key Features

-   **End-to-End Workflow:** Integrated pipeline from data acquisition to prioritized candidate list.
-   **Multi-Omics Integration:** Leverages TCGA, miRTarBase, STRING DB, Open Targets, and miRBase.
-   **Differential Expression Analysis:** Robust statistical analysis (Wilcoxon, t-test) with FDR correction.
-   **Network-Based Prioritization:** Identifies high-impact hub genes using graph theory (NetworkX).
-   **Druggability Assessment:** Filters targets based on known druggability characteristics.
-   **Structural Accessibility Prediction:** Optimizes binding site design using RNA secondary structure (`RNAfold`).
-   **Off-Target Safety Analysis:** Assesses potential binding to unintended sequences.
-   **Advanced Visualization:** Includes volcano plots, heatmaps (Seaborn, Plotly), UMAP, and interactive network graphs (PyVis).
-   **Pathway Enrichment:** Identifies biological pathways associated with target genes (g:Profiler).

## Project Structure

```
hackathonbro/
├── data/                     # Data directory
│   ├── annotations/          # Annotation files (miRBase, mappings)
│   ├── mirnet/               # miRNA target databases (miRTarBase, etc.)
│   ├── open_targets/         # Open Targets data (Parquet files)
│   ├── string_db/            # STRING protein interaction data
│   ├── tcga/                 # TCGA data (expression, clinical)
│   └── utr_sequences/        # Reference UTR sequences (e.g., HBB)
├── results/                  # Analysis results (CSVs, plots, reports)
│   ├── accessibility_analysis/ # Accessibility prediction outputs
│   └── networks/             # Network visualizations (HTML, GraphML)
├── scripts/                  # Main analysis and workflow scripts
│   ├── advanced_visualizations.py
│   ├── convert_mimat_to_mirna.py
│   ├── mirna_differential_expression.py
│   ├── off_target_check.py
│   ├── pathway_enrichment.py
│   └── target_prioritization_analysis.py
├── src/                      # Source code for core modules
│   ├── data/                 # Data fetching/processing modules (e.g., fetch_xena_data_api.py)
│   ├── models/               # Core computational models (e.g., rna_structure.py)
│   ├── utils/                # Utility functions (e.g., mirna_utils.py)
│   └── visualization/        # Visualization modules (e.g., structure_viz.py)
├── run_mirna_switch_workflow.py # Main orchestrator script for the platform
├── requirements.txt          # Project dependencies
└── README.md                 # This file
```
*(Note: Structure based on recent commit; some subdirectories might vary based on exact data organization)*

## Setup

1.  **Prerequisites:**
    *   Python 3.8+
    *   ViennaRNA package (See [ViennaRNA Installation](https://www.tbi.univie.ac.at/RNA/documentation.html#install))
    *   R environment (Optional, if using R-based analysis components like DESeq2)

2.  **Python Environment & Dependencies:**
    ```bash
    # Create and activate virtual environment (recommended)
    python3 -m venv venv
    source venv/bin/activate

    # Install Python packages
    pip install --upgrade pip
    pip install -r requirements.txt
    ```
    *Key Dependencies include:* `pandas`, `numpy`, `scipy`, `statsmodels`, `matplotlib`, `seaborn`, `plotly`, `networkx`, `pyvis`, `requests`, `biopython`, `umap-learn`, `gprofiler-official`, `pyarrow`, `xenaPython`.

3.  **R Dependencies (Optional):**
    ```R
    # Run in R console
    if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(c("DESeq2", "dplyr", "readr", "argparse", "ggplot2", "pheatmap"))
    ```

4.  **Data:**
    *   Download necessary data files (e.g., from STRING, miRBase, Open Targets) and place them in the corresponding `data/` subdirectories.
    *   Alternatively, use provided scripts (like `fetch_xena_data_api.py`) to acquire data where possible.

## Usage

The primary way to run the full analysis pipeline is using the main workflow script:

```bash
python run_mirna_switch_workflow.py
```

Individual analysis steps can also be run using the scripts in the `scripts/` directory, typically requiring specific input files generated by previous steps. Check the individual script arguments (`--help`) for details.

Example (Differential Expression):
```bash
# Ensure data/tcga/TCGA_LIHC_miRNA_processed.csv exists first
python scripts/mirna_differential_expression.py
```

Example (Target Prioritization - requires DE results and target data):
```bash
# Ensure results/mirna_de_significant_with_names.csv and data/mirnet/ files exist
python scripts/target_prioritization_analysis.py
```

## Scientific Background

MicroRNAs (miRNAs) are small non-coding RNAs crucial for post-transcriptional gene regulation. They bind primarily to the 3' UTR of target mRNAs, often leading to mRNA degradation or translational repression. Differential expression of miRNAs between healthy and diseased states provides opportunities for developing highly specific therapeutics. Effective miRNA targeting depends on factors like seed sequence complementarity, binding site accessibility (influenced by local RNA secondary structure), and target abundance.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

This platform utilizes data and tools from several sources, including:
- The Cancer Genome Atlas (TCGA)
- UCSC Xena Browser
- miRBase
- miRTarBase
- STRING Database
- Open Targets
- ViennaRNA Package
- g:Profiler
