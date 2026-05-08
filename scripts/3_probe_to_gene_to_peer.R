library(hthgu133a.db)
library(AnnotationDbi)
library(tidyverse)
library(limma)

# 1. Load Data
expr_matrix <- read.csv("GSE30169_Subject_Expression_RMA.csv", 
                        row.names = 1, 
                        check.names = FALSE)

# 2. Get Mapping
probe_map <- AnnotationDbi::select(hthgu133a.db, 
                                   keys = rownames(expr_matrix), 
                                   columns = "ENTREZID", 
                                   keytype = "PROBEID")

# 3. Clean Mapping
# We remove NAs and probes that map to multiple genes (ambiguous)
clean_map <- probe_map %>%
  filter(!is.na(ENTREZID)) %>%
  group_by(PROBEID) %>%
  filter(n_distinct(ENTREZID) == 1) %>%
  ungroup()

# 4. Align Matrix and Map
# Filter the matrix to only include the "clean" probes
expr_filtered <- expr_matrix[clean_map$PROBEID, ]

# 5. Handle Many-to-One (Many probes -> One gene)
# Use avereps on the filtered matrix, using the matching ENTREZIDs
# This returns a matrix where rownames are the unique Entrez IDs
expr_final <- avereps(expr_filtered, ID = clean_map$ENTREZID)

# Result Check
print(paste("Final number of unique genes:", nrow(expr_final)))

# Save the standard CSV (includes Entrez IDs and Sample Names)
write.csv(expr_final, "GSE30169_Expression_RMA_genes.csv", row.names = TRUE)

# Save the purely numerical table for peer
# We use write.table to easily remove headers and use a comma separator
# Transpose the matrix so genes become columns
expr_flipped <- t(expr_final)
write.table(expr_flipped, "GSE30169_Expr_RMA_genes_PEER_input.csv", 
            sep = ",", 
            col.names = FALSE, 
            row.names = FALSE)

identical(rownames(expr_flipped), as.character(peer_dbgap_order))





# Discovery below
# how many genes would be lost if i deleted probes that map to multiple genes, and thats its only probe
# Identify which probes are ambiguous (map to > 1 gene)
probe_counts <- table(probe_map$PROBEID)
ambiguous_probes <- names(probe_counts[probe_counts > 1])
unique_probes <- names(probe_counts[probe_counts == 1])

# Get the sets of genes
genes_on_ambiguous <- unique(probe_map$ENTREZID[probe_map$PROBEID %in% ambiguous_probes])
genes_on_unique <- unique(probe_map$ENTREZID[probe_map$PROBEID %in% unique_probes])

# Calculate "Lost" Genes: 
# Genes that appear in the ambiguous list but NEVER in the unique list
lost_genes <- setdiff(genes_on_ambiguous, genes_on_unique)

# Results
length(lost_genes)
