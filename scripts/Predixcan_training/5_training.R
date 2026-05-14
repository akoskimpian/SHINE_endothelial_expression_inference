source("gtex_v7_nested_cv_elnet.R")
"%&%" <- function(a,b) paste(a,b, sep='')

for (chrom in 1:22) {
  cat("Processing chromosome", chrom, "\n")


snp_annot_file  <- "processed/snp_annot.chr" %&% chrom %&% ".txt"
gene_annot_file <- "gene_annot.parsed.txt"
genotype_file   <- "processed/genotype.chr" %&% chrom %&% ".txt"
expression_file <- "transformed_expression.txt"
covariates_file <- "covariates.txt"
prefix          <- "Model_training"

main(snp_annot_file, gene_annot_file, genotype_file, expression_file, covariates_file, as.numeric(chrom), prefix, null_testing=FALSE)

