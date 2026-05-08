library(ggplot2)

# Path to your file

file_path <- "path/to/smartpca_output/pca_output.evec"

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

# Step A: Link PCA to Sample Metadata
# Using ID from PCA and SAMPLE_ID from sample_data
pca_with_samples <- merge(pca_data, sample_data, by.x = "ID", by.y = "SAMPLE_ID")

# Step B: Link the result to Subject Pheno
# Now we use SUBJECT_ID (which exists in both tables)
final_merged_data <- merge(pca_with_samples, subject_pheno, by = "SUBJECT_ID")
# --- Plot PC1 vs PC2 colored by RACE ---
ggplot(final_merged_data, aes(x = PC1, y = PC2, color = RACE)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(title = "PCA: PC1 vs PC2", 
       subtitle = "Colored by Race (self report)",
       color = "Race")

# --- Plot PC1 vs PC3 colored by RACE ---
ggplot(final_merged_data, aes(x = PC1, y = PC3, color = RACE)) +
  geom_point(alpha = 0.6, size = 2) +
  theme_minimal() +
  labs(title = "PCA: PC1 vs PC3", color = "Race")
