#!/usr/bin/env Rscript
# DESeq2 Analysis for miRNA Differential Expression in NASH-HCC vs Normal Liver
#
# This script performs differential expression analysis using DESeq2 on miRNA data
# for NASH-HCC vs Normal Liver comparison. It requires raw count data.
#
# Usage:
#   Rscript deseq2_analysis.R --input counts.csv --output results/deseq2_results.csv
#
# Author: Hackathon-Bio Team

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(argparse)
  library(ggplot2)
  library(pheatmap)
})

# Parse command line arguments
parser <- ArgumentParser(description="Perform DESeq2 analysis on miRNA expression data")

parser$add_argument("--input", dest="input_file", required=TRUE,
                   help="Input file with miRNA count data (CSV format)")
parser$add_argument("--id-col", dest="id_col", default="miRNA",
                   help="Column name for miRNA IDs (default: miRNA)")
parser$add_argument("--group-col", dest="group_col", default="Group",
                   help="Column name for sample groups (default: Group)")
parser$add_argument("--normal-label", dest="normal_label", default="Normal",
                   help="Label for normal liver samples (default: Normal)")
parser$add_argument("--disease-label", dest="disease_label", default="NASH-HCC",
                   help="Label for NASH-HCC samples (default: NASH-HCC)")
parser$add_argument("--output", dest="output_file", default="deseq2_results.csv",
                   help="Output file for DESeq2 results (default: deseq2_results.csv)")
parser$add_argument("--plots-dir", dest="plots_dir", default="plots",
                   help="Directory for output plots (default: plots)")
parser$add_argument("--padj-cutoff", dest="padj_cutoff", type="double", default=0.05,
                   help="Adjusted p-value cutoff for significance (default: 0.05)")
parser$add_argument("--log2fc-cutoff", dest="log2fc_cutoff", type="double", default=1.0,
                   help="Log2 fold change cutoff for significance (default: 1.0)")

args <- parser$parse_args()

# Ensure plots directory exists
dir.create(args$plots_dir, showWarnings = FALSE, recursive = TRUE)

# Function to load data
load_data <- function(input_file, id_col, group_col) {
  cat("Loading data from", input_file, "...\n")
  
  # Read the input file
  data <- read_csv(input_file, show_col_types = FALSE)
  
  # Verify required columns exist
  if(!id_col %in% colnames(data)) {
    stop(paste("ID column", id_col, "not found in the input file"))
  }
  
  if(!group_col %in% colnames(data)) {
    stop(paste("Group column", group_col, "not found in the input file"))
  }
  
  return(data)
}

# Function to prepare count matrix and sample data
prepare_deseq_input <- function(data, id_col, group_col, normal_label, disease_label) {
  # Extract sample data
  sample_data <- data %>%
    select(all_of(c(id_col, group_col))) %>%
    filter(!!sym(group_col) %in% c(normal_label, disease_label))
  
  # Convert to factor and relevel to make normal the reference
  sample_data[[group_col]] <- factor(sample_data[[group_col]], 
                                    levels = c(normal_label, disease_label))
  
  # Extract count data
  # In real data, would need to determine sample columns, for now assume all columns except id_col and group_col
  count_cols <- setdiff(colnames(data), c(id_col, group_col))
  
  if(length(count_cols) == 0) {
    stop("No count columns found in the data")
  }
  
  count_data <- data %>%
    select(all_of(c(id_col, count_cols))) %>%
    column_to_rownames(id_col)
  
  # Ensure counts are integers (DESeq2 requirement)
  count_data <- round(count_data)
  
  return(list(
    count_data = count_data,
    sample_data = sample_data
  ))
}

# Main function to run DESeq2 analysis
run_deseq2 <- function(count_data, sample_data, group_col) {
  cat("Running DESeq2 analysis...\n")
  
  # Create DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = sample_data,
    design = as.formula(paste("~", group_col))
  )
  
  # Run DESeq2 analysis
  dds <- DESeq(dds)
  
  # Extract results, with disease vs normal contrast
  res <- results(dds, contrast = c(group_col, levels(sample_data[[group_col]])[2], levels(sample_data[[group_col]])[1]))
  
  # Return the DESeqDataSet and results
  return(list(
    dds = dds,
    res = res
  ))
}

# Function to create plots
create_plots <- function(dds, res, plots_dir, padj_cutoff, log2fc_cutoff) {
  cat("Creating plots...\n")
  
  # 1. MA Plot
  png(file.path(plots_dir, "ma_plot.png"), width = 800, height = 600)
  DESeq2::plotMA(res, ylim = c(-5, 5), main = "MA Plot: NASH-HCC vs Normal Liver")
  abline(h = c(-log2fc_cutoff, log2fc_cutoff), col = "red", lty = 2)
  dev.off()
  
  # 2. Volcano Plot
  res_df <- as.data.frame(res)
  res_df$miRNA <- rownames(res_df)
  res_df$significant <- res_df$padj < padj_cutoff & abs(res_df$log2FoldChange) > log2fc_cutoff
  
  # Top 20 miRNAs to label
  top_mirnas <- res_df %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    head(20) %>%
    pull(miRNA)
  
  res_df$label <- ifelse(res_df$miRNA %in% top_mirnas, res_df$miRNA, NA)
  
  p <- ggplot(res_df %>% filter(!is.na(padj)), 
              aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red")) +
    geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    geom_text(aes(label = label), hjust = -0.3, size = 3, check_overlap = TRUE) +
    labs(
      title = "Volcano Plot: NASH-HCC vs Normal Liver",
      x = "log2 Fold Change",
      y = "-log10(adjusted p-value)"
    ) +
    theme_minimal()
  
  ggsave(file.path(plots_dir, "volcano_plot.png"), p, width = 10, height = 8)
  
  # 3. Sample Clustering Heatmap
  # Transform counts for visualization
  vst <- vst(dds, blind = FALSE)
  
  # Get top variable miRNAs
  topVarGenes <- head(order(rowVars(assay(vst)), decreasing = TRUE), 50)
  
  # Create heatmap
  mat <- assay(vst)[topVarGenes, ]
  rownames(mat) <- rownames(assay(vst))[topVarGenes]
  
  # Prepare annotation
  annotation_col <- data.frame(
    Group = colData(vst)[[group_col]],
    row.names = colnames(mat)
  )
  
  pheatmap(mat, 
           annotation_col = annotation_col,
           show_rownames = TRUE,
           show_colnames = FALSE,
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation",
           filename = file.path(plots_dir, "heatmap_top50_variable.png"),
           width = 10,
           height = 12)
  
  # 4. PCA Plot
  p <- plotPCA(vst, intgroup = group_col) +
    theme_minimal() +
    ggtitle("PCA Plot: NASH-HCC vs Normal Liver")
  
  ggsave(file.path(plots_dir, "pca_plot.png"), p, width = 8, height = 6)
  
  cat("Plots saved to", plots_dir, "\n")
}

# Main execution
main <- function() {
  # Load data
  data <- load_data(args$input_file, args$id_col, args$group_col)
  
  # Prepare DESeq2 input
  cat("Preparing data for DESeq2...\n")
  input <- prepare_deseq_input(data, args$id_col, args$group_col, args$normal_label, args$disease_label)
  
  # Print sample counts
  sample_counts <- table(input$sample_data[[args$group_col]])
  cat("\nSample counts:\n")
  print(sample_counts)
  
  if(any(sample_counts < 3)) {
    warning("At least one group has fewer than 3 samples, which may affect statistical power")
  }
  
  # Run DESeq2
  deseq_out <- run_deseq2(input$count_data, input$sample_data, args$group_col)
  
  # Create plots
  create_plots(deseq_out$dds, deseq_out$res, args$plots_dir, args$padj_cutoff, args$log2fc_cutoff)
  
  # Process and save results
  cat("Processing results...\n")
  res_df <- as.data.frame(deseq_out$res)
  res_df$miRNA <- rownames(res_df)
  
  # Reorder columns
  res_df <- res_df %>%
    select(miRNA, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
  
  # Save results
  write_csv(res_df, args$output_file)
  cat("Results saved to", args$output_file, "\n")
  
  # Summary
  sig_up <- sum(res_df$padj < args$padj_cutoff & res_df$log2FoldChange > args$log2fc_cutoff, na.rm = TRUE)
  sig_down <- sum(res_df$padj < args$padj_cutoff & res_df$log2FoldChange < -args$log2fc_cutoff, na.rm = TRUE)
  
  cat("\nAnalysis summary:\n")
  cat("Total miRNAs analyzed:", nrow(res_df), "\n")
  cat("Significant miRNAs (adjusted p <", args$padj_cutoff, "and |log2FC| >", args$log2fc_cutoff, "):", sig_up + sig_down, "\n")
  cat("  - Upregulated in NASH-HCC:", sig_up, "\n")
  cat("  - Downregulated in NASH-HCC:", sig_down, "\n")
  
  # Show top significant miRNAs
  cat("\nTop 10 differentially expressed miRNAs:\n")
  top10 <- res_df %>%
    filter(!is.na(padj)) %>%
    filter(padj < args$padj_cutoff & abs(log2FoldChange) > args$log2fc_cutoff) %>%
    arrange(padj) %>%
    head(10)
  
  if(nrow(top10) > 0) {
    for(i in 1:nrow(top10)) {
      direction <- ifelse(top10$log2FoldChange[i] > 0, "↑", "↓")
      cat(sprintf("%s %s (log2FC: %.2f, FDR: %.2e)\n", 
                 top10$miRNA[i], direction, top10$log2FoldChange[i], top10$padj[i]))
    }
  } else {
    cat("No significant differentially expressed miRNAs found.\n")
  }
}

# Run the main function
main() 