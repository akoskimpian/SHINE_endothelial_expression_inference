# SHINE_endothelial_expression_inference
SNP-based inference of endothelial gene expression for pathway analysis in shock-induced endotheliopathy

The endothelium is a dynamic organ central to vascular homeostasis, and its dysfunction is implicated in a wide range of diseases. In acute critical illness, shock-induced endotheliopathy (SHINE) has been proposed as a unifying mechanism linking distinct shock states, yet direct clinical assessment of the endothelium remains infeasible in critically ill patients. To overcome this, we leverage patient genomic variation  to infer endothelial gene expression from SNP data, reconstructing patient-specific molecular states. These inferred expression profiles are then used for pathway-level analysis to explore biological processes potentially underlying endothelial dysfunction in SHINE, with the broader goal of identifying molecular signatures that could support risk stratification and guide future therapies.

This repository currently contains the preprocessing pipeline for this work; the project is ongoing and the codebase is updated as development continues.


# Script descriptions
The input of the scripts originates from the accession GSE30169, and can be accessed through the dbGaP database through the accession phs002057.v1.p1.
The majority of these scripts are data exploratory and data wrangling in nature.
The scripts are numbered in order of execution.

1_dealing_with_metadata.R
Loads and parses the dbGaP metadata files. Extracts subject-level sex information and formats it for use downstream.

2_microarray_processing.R
Reads raw Affymetrix CEL files and applies RMA normalization. Filters to control samples only. Some subjects have multiple microarray replicates, these are resolved by averaging their expression values into a single profile per subject. The result is a probe-level expression matrix (subjects × probes) saved as a CSV.

3_snp_processing.R
Matches subjects in the expression data to their corresponding samples in the genotype VCF file. Some subjects have multiple genotype samples, analysis confirms these are technical replicates, so one sample per subject is kept. The final list of sample IDs to retain is saved as keep_samples.txt.

4_probe_to_gene_to_peer.R
Converts the probe-level expression matrix to gene-level by mapping probes to Entrez IDs using the hthgu133a.db annotation package. Probes mapping to multiple genes are discarded. Probes that map to the same gene are averaged. The final gene-level matrix is saved both as a standard CSV and as a transposed, header-free version ready for PEER input.

5_create_peer_covariates.R
Assembles the covariate matrix to be passed to PEER alongside the expression data. Combines the top 5 genotype PCs (computed via smartpca) with sex, aligned to match the sample order in the expression matrix. Saved as a header-free CSV ready for PEER input.

smartpca_visualisation.R
Parses the .evec output file from smartpca and visualizes the genotype PCA results. Adds demographic annotations from metadata. Plots PC1 vs PC2 and PC1 vs PC3, colored by annotated race, to check for population stratification before using the PCs as covariates in PEER.

peer_visualisation.R
Visualizes the output of PEER by plotting factor relevance (1/Alpha) for each factor. Produces both a raw relevance plot and a log-scaled bar chart.
