# Function that reads in metadata files
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
#  Load the tables
geo_mapping   <- read_dbgap_v3("path/to/dbgap/dataset/data/phs002057.v1.p1/phs002057.v1.pht010401.v1.p1.ECs_GEO.MULTI.txt.gz")
sample_attr   <- read_dbgap_v3("path/to/dbgap/dataset/data/phs002057.v1.p1/phs002057.v1.pht010403.v1.p1.c1.ECs_Sample_Attributes.GRU-PUB.txt.gz")
subject_pheno <- read_dbgap_v3("path/to/dbgap/dataset/data/phs002057.v1.p1/phs002057.v1.pht010402.v1.p1.c1.ECs_Subject_Phenotypes.GRU-PUB.txt.gz")
subject_multi <- read_dbgap_v3("path/to/dbgap/dataset/data/phs002057.v1.p1/phs002057.v1.pht010399.v1.p1.ECs_Subject.MULTI.txt.gz")
sample_data <- read_dbgap_v3("path/to/dbgap/dataset/data/phs002057.v1.p1/phs002057.v1.pht010400.v1.p1.ECs_Sample.MULTI.txt.gz")

####Load in ancestry PCs
file_path <- "C:/Users/AKIM0004/Documents/azure_docker_pca_peer/smartpca_input/pca_output.evec"

# Read the file (skip eigenvalue line)

pca_data <- read.table(file_path, header = FALSE, skip = 1, stringsAsFactors = FALSE)

# Assign column names

# First column = ID, last column = label (e.g., NOTHING)

colnames(pca_data) <- c("ID", paste0("PC", 1:(ncol(pca_data)-2)), "Label")

# Remove the label column if not needed

pca_data <- pca_data[, -ncol(pca_data)]

# Convert PCs to numeric (just in case)

pca_data[, -1] <- lapply(pca_data[, -1], as.numeric)

# Change "Sample87:Sample87" to just "Sample87"
pca_data$ID <- sub(":.*", "", pca_data$ID)


#Read in expression table
expr_rma <- read.csv("GSE30169_Subject_Expression_RMA.csv",
                     check.names = FALSE,
                     header = TRUE,
                     row.names = 1)

# Get the expr_rma column order (these are dbGaP Sample IDs) 
peer_dbgap_order <- colnames(expr_rma)

# Map dbGaP Sample IDs → Sample IDs using sample_data 
sample_data_filtered <- sample_data[sample_data$SAMPLE_ID %in% pca_data$ID, ]
dbgap_to_sample <- setNames(sample_data_filtered$SAMPLE_ID,
                            sample_data_filtered$dbGaP_Subject_ID)

peer_sample_order <- dbgap_to_sample[peer_dbgap_order]

# Check for failed conversions 
if (any(is.na(peer_sample_order))) {
  cat("WARNING:", sum(is.na(peer_sample_order)), "dbGaP Sample IDs could not be converted to Sample IDs\n")
} else {
  cat("All", length(peer_sample_order), "dbGaP Sample IDs successfully converted to Subject IDs")
}

# Check alignment with pca_data 
same_set   <- setequal(peer_sample_order, pca_data$ID)
same_order <- identical(as.character(peer_sample_order), as.character(pca_data$ID))
cat("Same set:  ", same_set, "\n")
cat("Same order:", same_order, "\n")

if (same_set & !same_order) {
  mismatched_positions <- which(as.character(peer_sample_order) != as.character(pca_data$ID))
  cat("Order differs at", length(mismatched_positions), "positions\n")
}

# Reorder pca_data to match expr_rma column order 
pca_data_reordered <- pca_data[match(peer_sample_order, pca_data$ID), ]

same_set   <- setequal(pca_data_reordered$ID, peer_sample_order)
same_order <- identical(as.character(pca_data_reordered$ID), as.character(peer_sample_order))
cat("Same set:  ", same_set, "\n")
cat("Same order:", same_order, "\n")

# Load sex data (first column is dbGaP Subject ID) 
sex_data <- read.table("sex_data.ind",
                       header = FALSE,
                       col.names = c("dbGaP_Subject_ID", "Sex", "Extra"))

# Convert sex to numeric (M=0, F=1) 
sex_data$Sex_numeric <- ifelse(sex_data$Sex == "M", 0, 1)


# Reorder sex_data to match peer_dbgap_order using Subject IDs directly 
sex_data_reordered <- sex_data[match(as.character(peer_dbgap_order), sex_data$dbGaP_Subject_ID), ]

# Verify 
cat("Any NAs after reordering?", any(is.na(sex_data_reordered$dbGaP_Subject_ID)), "\n")
cat("Order matches peer_dbgap_order:",
    identical(as.character(sex_data_reordered$dbGaP_Subject_ID),
              as.character(peer_dbgap_order)), "\n")

# Combine PC1-5 and sex into covariate matrix 
peer_covariates <- cbind(pca_data_reordered[, c("PC1", "PC2", "PC3", "PC4", "PC5")],
                         Sex = sex_data_reordered$Sex_numeric)
head(peer_covariates)
dim(peer_covariates)

write.table(peer_covariates,
            file = "GSE30169_PEER_Covariates.csv",
            sep = ",",
            row.names = FALSE,
            col.names = FALSE)
