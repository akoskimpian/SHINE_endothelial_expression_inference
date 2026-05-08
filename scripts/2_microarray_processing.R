library(tidyverse)
library(affy)
library(limma)

# used in tandem with data_exploration! some variables and functions there

# affy documentation:
# https://bioconductor.posit.co/packages/3.19/bioc/html/affy.html


#####
# GSE30169 (the only one)
####

#raw cel data
raw_30169 <- ReadAffy(celfile.path = "path/to/GSE30169_control_cels/")

meta_30169 <- read.csv("path/to/GSE30169_metadata.csv")

eset <- rma(raw_30169)


### DEALING with replicates - inspection
current_ids_raw <- colnames(exprs(eset))
current_ids_clean <- sub("\\.cel.*", "", current_ids_raw, ignore.case = TRUE)

gsm_subset <- gsm_titles[current_ids_clean]

# keep only samples whose gsm_titles contain "control" (case-insensitive)
ctrl_mask <- grepl("control", gsm_subset, ignore.case = TRUE)


# subset the gsm/title vectors and the ExpressionSet to controls only
gsm_subset <- gsm_subset[ctrl_mask]
current_ids_clean <- current_ids_clean[ctrl_mask]

# subset the ExpressionSet columns to match the same order/selection
eset <- eset[, colnames(exprs(eset))[ctrl_mask]]



groups <- sub("_rep.*", "", gsm_subset)
expr_collapsed <- sapply(unique(groups), function(g){
  cols <- which(groups == g)
  rowMeans(exprs(eset)[, cols, drop = FALSE])
})
dim(expr_collapsed)

##checking replicate and sample counts
sample_counts_full <- table(meta_30169$dbGaP_Subject_ID)
table(sample_counts_full)
current_gsm_ids <- sub("\\.cel.*", "", colnames(exprs(eset)), ignore.case = TRUE)
meta_subset <- meta_30169[meta_30169$GEO_ACCESSION %in% current_gsm_ids, ]
sample_counts_subset <- table(meta_subset$dbGaP_Subject_ID)
summary_table <- table(sample_counts_subset)
print(summary_table)

#Averaging expression replicates 
#USED in the end
subject_ids <- as.character(meta_aligned$dbGaP_Subject_ID)

expr_averaged <- sapply(unique(subject_ids), function(s) {
  cols <- which(subject_ids == s)
  rowMeans(ex[, cols, drop = FALSE])
})

colnames(expr_averaged) <- unique(subject_ids)

output_path <- "GSE30169_Subject_Expression_RMA.csv"
write.csv(expr_averaged, file = output_path, row.names = TRUE)



# Basic dimensions
dim(old)
dim(new)


# column names
head(colnames(old))
head(colnames(new))


##### Dropped one replicate where there was 2 and getting expression table 
# NOT USED
# during SNP processing, we chose one sample to go on with.
# however for expression data, the best practice is to take the average expression, as done above

# 1. Get the expression matrix
ex <- exprs(eset)

# 2. Clean the column names (GSM IDs)
# This strips .CEL, .cel, .gz, etc.
current_gsm_ids <- gsub("\\.cel.*|\\.gz.*", "", colnames(ex), ignore.case = TRUE)

# 3. Match metadata to the expression columns
# We use match to make sure meta_aligned is in the EXACT same order as colnames(ex)
matching_indices <- match(current_gsm_ids, meta_30169$GEO_ACCESSION)

# --- IMPORTANT CHECK: Did every column find a match? ---
if (any(is.na(matching_indices))) {
  missing_count <- sum(is.na(matching_indices))
  stop(paste("Error:", missing_count, "columns in your data don't match the GEO_ACCESSIONs in your metadata! Check for typos or extra suffixes."))
}

meta_aligned <- meta_30169[matching_indices, ]
subject_ids <- as.character(meta_aligned$dbGaP_Subject_ID)

# 4. Select One Sample per Subject (Keep First Instance) 
# but use Sample ID for the header
subject_ids <- as.character(meta_aligned$dbGaP_Subject_ID)
sample_ids  <- as.character(meta_aligned$dbGaP_Sample_ID)

keep_idx    <- !duplicated(subject_ids)
expr_final  <- ex[, keep_idx]

# Name the columns using the Sample IDs of the kept samples
colnames(expr_final) <- sample_ids[keep_idx]

# 5. Final Check
print(paste("Original columns:", ncol(ex)))
print(paste("Final subjects:", ncol(expr_final)))

# Define the output path
output_path <- "C:/Users/AKIM0004/Documents/GSE30169_Subject_Expression_RMA.csv"

# Save as a CSV
# Using row.names=TRUE ensures your Probe/Gene IDs are saved in the first column
write.csv(expr_final, file = output_path, row.names = TRUE)
