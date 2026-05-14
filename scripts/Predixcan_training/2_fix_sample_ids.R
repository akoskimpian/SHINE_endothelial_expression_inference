library(tidyverse)

read_dbgap_v3 <- function(path) {
  # Read the whole file into a temporary object
  raw_data <- read.delim(path, comment.char = "#", header = FALSE, stringsAsFactors = FALSE)
  
  # Find the row that contains 'dbGaP_Sample_ID' or 'dbGaP_Subject_ID'
  # We search the first 5 rows for these keywords
  header_row_index <- which(apply(raw_data[1:5, ], 1, function(x) any(grepl("dbGaP_Sample_ID|dbGaP_Subject_ID", x))))[1]
  
  # If we found it, set it as the header
  if(!is.na(header_row_index)) {
    colnames(raw_data) <- as.character(raw_data[header_row_index, ])
    # Remove the header row and any variable ID rows above it
    raw_data <- raw_data[(header_row_index + 1):nrow(raw_data), ]
  }
  
  # Clean up whitespace
  colnames(raw_data) <- trimws(colnames(raw_data))
  return(raw_data)
}

# Load the table
sample_data <- read_dbgap_v3("path/to/dataset/phs002057.v1.p1/phs002057.v1.pht010400.v1.p1.ECs_Sample.MULTI.txt.gz")

expr <- read.csv("C:/Users/AKIM0004/Documents/Predixcan/expression_data/GSE30169_Expression_RMA_genes.csv",
                 header=TRUE)

#To check how the genotypeheader looks like
geno_header <- scan("processed/genotype.chr1.txt",
                    what="", nlines=1, sep="\t")

geno_samples <- geno_header[-1]

expr_samples <- colnames(expr)




# Strip the "X" from expression column names
# Skip the first column (the gene names column)
colnames(expr)[-1] <- gsub("^X", "", colnames(expr)[-1])

# Create a named vector for mapping (Subject_ID -> Sample_ID)
id_map <- setNames(as.character(sample_data$SAMPLE_ID), 
                   as.character(sample_data$dbGaP_Subject_ID))

# Overwrite expression headers with Sample IDs
# Keep the first column name and map the rest
current_expr_ids <- colnames(expr)[-1]
new_sample_ids <- id_map[current_expr_ids]

# Update the column names
colnames(expr)[-1] <- new_sample_ids

# Check for duplicates and count
num_unique <- length(unique(colnames(expr)[-1]))
cat("Number of unique Sample IDs:", num_unique, "\n")


# Define file path
output_path <- "GSE30169_Expression_RMA_genes.csv"

# Write the updated data frame to CSV
write.csv(expr, file = output_path, row.names = FALSE, quote = FALSE)

