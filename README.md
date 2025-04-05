# miRNA Target Optimization Platform

A computational platform for designing mRNA therapeutics with enhanced tissue specificity through microRNA-based targeting.

## Project Overview

This platform addresses a critical challenge in drug design: targeted delivery. By leveraging microRNA expression patterns that differentiate between diseased and healthy tissues, we can create mRNA therapeutics that function selectively where needed while sparing healthy cells.

The platform focuses on NASH-HCC (Non-Alcoholic Steatohepatitis-driven Hepatocellular Carcinoma) as a proof-of-concept. With approximately 25% of the world's population suffering from Non-Alcoholic Fatty Liver Disease, this represents an urgent medical need and significant commercial opportunity.

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

- **Differential miRNA Analysis**: Identify disease-specific miRNA signatures
- **Structure-based Site Optimization**: Leverage ViennaRNA for RNA structure prediction
- **Visualization Tools**: Interpret complex RNA interactions visually

## Project Structure

```
hackathonbro/
├── src/                 # Source code
│   ├── data/            # Data handling modules
│   │   ├── fetch_xena_data.py     # TCGA data acquisition
│   │   ├── fetch_geo_data.py      # GEO data acquisition
│   │   ├── deseq2_analysis.R      # R-based DESeq2 analysis
│   │   └── run_nash_hcc_analysis.py # End-to-end workflow
│   ├── models/          # Core computational models
│   ├── utils/           # Utility functions
│   └── visualization/   # Visualization modules
├── data/                # Data directory
│   ├── raw/             # Raw data
│   └── processed/       # Processed data
├── results/             # Analysis results
│   └── plots/           # Visualization plots
├── tools/               # Helper scripts and tools
├── requirements.txt     # Project dependencies
└── README.md            # This file
```

## Setup

1. Ensure you have Python 3.8+ installed
2. Install ViennaRNA package (required for RNA structure prediction):
   ```bash
   # macOS (using Homebrew)
   brew tap brewsci/bio
   brew install brewsci/bio/viennarna
   
   # Ubuntu/Debian
   sudo apt-get install viennarna
   ```

3. Install required Python packages:
   ```bash
   # Create and activate virtual environment
   python3 -m venv venv
   source venv/bin/activate
   
   # Install dependencies
   pip install -r requirements.txt
   ```

4. Install R dependencies (optional, for DESeq2 analysis):
   ```R
   if (!require("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   
   BiocManager::install(c("DESeq2", "dplyr", "readr", "argparse", "ggplot2", "pheatmap"))
   ```

## Usage

### Target Prioritizer

Standard analysis:
```bash
python src/target_prioritizer.py --input path/to/data.csv --output results/
```

New end-to-end workflow:
```bash
python src/data/run_nash_hcc_analysis.py --source xena --analysis python
```

### Accessibility Predictor
```bash
python src/accessibility_predictor.py --utr sample_utr.fa --mirna miR-122
```

## Scientific Background

MicroRNAs (miRNAs) are small non-coding RNAs that regulate gene expression through binding to messenger RNAs (mRNAs). Their expression patterns differ between diseased and healthy tissues, offering a natural mechanism for targeted therapeutics.

Key aspects:
- miRNA targeting relies on seed region (nucleotides 2-8) complementarity
- Secondary structure and accessibility significantly impact binding efficiency
- Canonical seed matches include: 8mer, 7mer-m8, and 7mer-A1 sites

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- TCGA (The Cancer Genome Atlas) for providing data
- UCSC Xena Browser for data access
- Gene Expression Omnibus (GEO) for additional datasets
