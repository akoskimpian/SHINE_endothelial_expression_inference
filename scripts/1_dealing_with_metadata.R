library(dplyr)
library(tidyverse)

#  Set Working Directory
setwd("path/to/dbgap/dataset/data")

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
geo_mapping   <- read_dbgap_v3("phs002057.v1.p1/phs002057.v1.pht010401.v1.p1.ECs_GEO.MULTI.txt.gz")
sample_attr   <- read_dbgap_v3("phs002057.v1.p1/phs002057.v1.pht010403.v1.p1.c1.ECs_Sample_Attributes.GRU-PUB.txt.gz")
subject_pheno <- read_dbgap_v3("phs002057.v1.p1/phs002057.v1.pht010402.v1.p1.c1.ECs_Subject_Phenotypes.GRU-PUB.txt.gz")
subject_multi <- read_dbgap_v3("phs002057.v1.p1/phs002057.v1.pht010399.v1.p1.ECs_Subject.MULTI.txt.gz")
sample_data <- read_dbgap_v3("phs002057.v1.p1/phs002057.v1.pht010400.v1.p1.ECs_Sample.MULTI.txt.gz")



# Getting sex data since it is needed for PCA
# sex data is in subject_multi
ind_data <- subject_multi %>%
  # 1. Select the ID and Sex columns
  select(dbGaP_Subject_ID, SEX) %>%
  # 2. Format columns for smartpca
  mutate(
    # smartpca expects M/F/U or 1/2. 
    # If your data is "Male"/"Female", we should map it to M/F.
    SEX = case_when(
      SEX %in% c("Male", "M", "1") ~ "M",
      SEX %in% c("Female", "F", "2") ~ "F",
      TRUE ~ "U"
    ),
    # 3. Add a placeholder population column (everyone gets the same label)
    POP = "NOTHING"
  )

# Write to a tab-delimited file without headers
# The .ind format is: ID  SEX  POP
write_delim(ind_data, "sex_data.ind", delim = "\t", col_names = FALSE)
