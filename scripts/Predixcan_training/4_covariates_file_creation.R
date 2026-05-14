peer_factors = read.csv(file = "peer_out/X.csv", header = FALSE)

gene_exp_transpose <- read.table("transformed_expression.txt", 
                                 header=TRUE,
                                 row.names=1,
                                 sep="\t")
#Set the column names for the PEER factors as the subject IDs
colnames(peer_factors) = rownames(gene_exp_transpose)

# Write out covariates matrix
write.table(peer_factors, file = "covariates.txt", sep = "\t",
            row.names = TRUE)
