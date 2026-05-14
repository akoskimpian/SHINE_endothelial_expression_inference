# SHINE_endothelial_expression_inference
SNP-based inference of endothelial gene expression for pathway analysis in shock-induced endotheliopathy

The endothelium is a dynamic organ central to vascular homeostasis, and its dysfunction is implicated in a wide range of diseases. In acute critical illness, shock-induced endotheliopathy (SHINE) has been proposed as a unifying mechanism linking distinct shock states, yet direct clinical assessment of the endothelium remains infeasible in critically ill patients. To overcome this, we leverage patient genomic variation  to infer endothelial gene expression from SNP data, reconstructing patient-specific molecular states. These inferred expression profiles are then used for pathway-level analysis to explore biological processes potentially underlying endothelial dysfunction in SHINE, with the broader goal of identifying molecular signatures that could support risk stratification and guide future therapies.

This repository currently contains the preprocessing pipeline for this work; the project is ongoing and the codebase is updated as development continues.


# Script descriptions
The input of the scripts originates from the accession GSE30169, and can be accessed through the dbGaP database through the accession phs002057.v1.p1.
The scripts are numbered in order of execution where relevant

## initial_preprocessing
The majority of these scripts are data exploratory and data wrangling in nature.

### 1_dealing_with_metadata.R
Loads and parses the dbGaP metadata files. Extracts subject-level sex information and formats it for use downstream.

### 2_microarray_processing.R
Reads raw Affymetrix CEL files and applies RMA normalization. Filters to control samples only. Some subjects have multiple microarray replicates, these are resolved by averaging their expression values into a single profile per subject. The result is a probe-level expression matrix (subjects × probes) saved as a CSV.

### 3_snp_processing.R
Matches subjects in the expression data to their corresponding samples in the genotype VCF file. Some subjects have multiple genotype samples, analysis confirms these are technical replicates, so one sample per subject is kept. The final list of sample IDs to retain is saved as keep_samples.txt.

### 4_probe_to_gene_to_peer.R
Converts the probe-level expression matrix to gene-level by mapping probes to Entrez IDs using the hthgu133a.db annotation package. Probes mapping to multiple genes are discarded. Probes that map to the same gene are averaged. The final gene-level matrix is saved both as a standard CSV and as a transposed, header-free version ready for PEER input.

### 5_create_peer_covariates.R
Assembles the covariate matrix to be passed to PEER alongside the expression data. Combines the top 5 genotype PCs (computed via smartpca) with sex, aligned to match the sample order in the expression matrix. Saved as a header-free CSV ready for PEER input.

### smartpca_visualisation.R
Parses the .evec output file from smartpca and visualizes the genotype PCA results. Adds demographic annotations from metadata. Plots PC1 vs PC2 and PC1 vs PC3, colored by annotated race, to check for population stratification before using the PCs as covariates in PEER.

### peer_visualisation.R
Visualizes the output of PEER by plotting factor relevance (1/Alpha) for each factor. Produces both a raw relevance plot and a log-scaled bar chart.

## Predixcan_training

The guide I have followed: https://github.com/hakyimlab/PredictDB-Tutorial/blob/master/README.md

### 1_bed_to_dosages.sh
Converts PLINK binary format genotype files to per-chromosome dosage and SNP annotation files.

### 2_fix_sample_ids.R
Maps expression matrix column names from subject IDs to sample IDs using the dbGaP metadata. The updated expression matrix is saved as a CSV.

### 3_gene_annotation_creation.R
Builds a gene annotation file from the Entrez IDs in the expression matrix using the `hthgu133a.db` package. Retrieves genomic coordinates, gene symbol, and chromosome for each gene, filtering to autosomal entries with complete coordinate information and one row per  gene. Also filters the expression matrix to annotated genes only and transposes it to a samples × genes format ready for PrediXcan training input.

### 4_covariates_file_creation.R
Reads the PEER factor matrix output and assigns sample IDs. Saves the labelled covariate matrix file ready for PrediXcan training.

### 5_training.R / gtex_v7_nested_cv_elnet.R
Trains cis-genetic expression prediction models across all 22 autosomes using the
PrediXcan elastic net framework (`gtex_v7_nested_cv_elnet.R`). For each chromosome, palindromic SNPs (A/T and C/G) are removed and a MAF
filter of 1% is applied. Each gene is modelled using SNPs within a 1 Mb *cis*-window,
with expression values first residualised against the PEER covariate matrix. Model
performance is estimated via nested cross-validation (5 outer folds, 10 inner folds) and
p-values are combined across folds using Fisher's and Stouffer's methods. Per-chromosome
outputs are written to `summary/`, `weights/`, and `covariances/` directories.
