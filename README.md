# R-Studio

Project Overview
This project involves the analysis of RNA sequencing (RNA-seq) data to explore gene expression differences across various conditions. The primary objectives are to identify differentially expressed genes, perform Gene Ontology (GO) enrichment analysis, and visualize the results using plots such as PCA, MA plots, and heatmaps.

Project Structure
The project is organized into several directories and files:

raw_data: Contains raw RNA-seq count data and annotation files.

scripts: Includes R scripts for data analysis, including DESeq2 analysis, PCA, MA plots, heatmaps, and GO enrichment analysis.

results: Holds output files from the analysis, such as DESeq2 results, PCA plots, and GO enrichment plots.

logs: Stores log files from the analysis process.

README.md: This file provides an overview of the project.

Key Files
em.csv: Raw count data.

annotation.csv: Gene annotation data.

col_data.csv: Sample metadata.

DESeq2_analysis.R: Script for performing DESeq2 analysis.

GO_enrichment.R: Script for GO enrichment analysis.

Analysis Workflow
Data Preparation: Merge count data with annotation data.

DESeq2 Analysis: Perform differential expression analysis using DESeq2.

Visualization:

PCA plot to visualize sample clustering.

MA plot and volcano plot to identify differentially expressed genes.

Heatmap to display gene expression patterns.

GO Enrichment Analysis: Identify enriched GO terms using clusterProfiler.

Software and Libraries
R: Primary programming language.

DESeq2: Package for differential expression analysis.

ggplot2: Package for data visualization.

clusterProfiler: Package for GO enrichment analysis.

Running the Analysis
Clone the repository.

Install necessary R packages.

Run the DESeq2_analysis.R script for differential expression analysis.

Run the GO_enrichment.R script for GO enrichment analysis.
