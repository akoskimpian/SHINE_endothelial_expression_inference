library(hthgu133a.db)
library(AnnotationDbi)
library(tidyverse)

expr <- read.csv("GSE30169_Expression_RMA_genes.csv",
                 row.names = 1)
# Get all the columns needed for gene_annot in one query
gene_annot_map <- AnnotationDbi::select(hthgu133a.db,
                                        keys = rownames(expr),  # these are already Entrez IDs
                                        columns = c("SYMBOL", "CHRLOC", "CHRLOCEND", "CHR", "GENENAME"),
                                        keytype = "ENTREZID")

gene_annot <- gene_annot_map %>%
  filter(!is.na(CHR), !is.na(CHRLOC), !is.na(CHRLOCEND)) %>%
  mutate(
    start     = abs(CHRLOC),     # negative = minus strand
    end       = abs(CHRLOCEND),
    chr       = CHR,
    gene_id   = ENTREZID,
    gene_name = SYMBOL,
    gene_type = "protein_coding"
  ) %>%
  group_by(gene_id) %>%
  slice(1) %>%                     # one row per gene
  ungroup() %>%
  select(gene_name, gene_id, start, end, chr, gene_type) %>%
  filter(chr %in% as.character(1:22)) # autosomes only

cat("Genes in expression matrix:", nrow(expr), "\n")
cat("Genes with annotation:", nrow(gene_annot), "\n")
cat("Overlap:", length(intersect(rownames(expr), gene_annot$gene_id)), "\n")

# Check chr format — should be "1" not "chr1"
head(gene_annot$chr)

# Save gene annotation
write.table(gene_annot,
            "gene_annot.parsed.txt",
            sep="\t", row.names=FALSE, quote=FALSE)

# Format expression for training
# Filter to only genes that have annotation
expr_matched <- expr[rownames(expr) %in% gene_annot$gene_id, ]

# training script needs rows=samples, cols=genes
expr_t <- as.data.frame(t(expr_matched))


write.table(expr_t,
            "transformed_expression.txt",
            sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)








