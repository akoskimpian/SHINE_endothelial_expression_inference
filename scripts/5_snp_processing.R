library(vcfR)
library(VariantAnnotation)
library(Rsamtools)
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
sample_data <- read_dbgap_v3("phs002057.v1.p1/phs002057.v1.pht010400.v1.p1.ECs_Sample.MULTI.txt.gz")




vcf_path <- "path/to/dbgap/dataset/phg001489.v1.EndothelialCells.genotype-calls-vcf.Genome-Wide_Human_SNP_Array_6_0.c1.GRU-PUB/dbGAP-vari-gene-submission2.vcf.gz"
vcf_bgz <- gsub(".gz$", ".bgz", vcf_path)

headers <- scanVcfHeader(vcf_bgz)

samples <- samples(headers)
head(samples)

#### checking if we have snps for our expression - THEY are here, skip this block

# 1. Read only the header of your CSV
# row.names=1 treats the first column (Gene IDs) as row names, 
# so colnames() will give you just the Subject IDs.
csv_headers <- colnames(read.csv("GSE30169_Subject_Expression_RMA.csv", 
                                 check.names = FALSE, 
                                 row.names = 1, 
                                 nrows = 1))

# Map those dbGaP_Subject_IDs to SAMPLE_IDs
# We look for matches in the 'dbGaP_Subject_ID' column of the mapping table
matched_metadata <- sample_data[sample_data$dbGaP_Sample_ID %in% csv_headers, ]

# Get  VCF Sample List 
vcf_samples <- samples(headers)

# Final Comparison
# Check if the 'SAMPLE_ID' from the mapping exists in the VCF
# (Note: If VCF says "Sample87" but mapping says "87", we use paste0 below)
present_in_vcf <- matched_metadata$SAMPLE_ID %in% vcf_samples

# The Results 
found_ids <- matched_metadata$SAMPLE_ID[present_in_vcf]
missing_ids <- matched_metadata$SAMPLE_ID[!present_in_vcf]

cat("Unique Subjects in CSV: ", length(unique(csv_headers)), "\n")
cat("Samples found in VCF: ", length(unique(found_ids)), "\n")

#END OF SKIP

#### FILTERING SNP SAMPLES


target_samples <- sample_data$SAMPLE_ID[sample_data$dbGaP_Sample_ID %in% csv_headers]
keep_samples <- intersect(target_samples, samples)

# This looks at length of chromosomes to check ref genome: verifies its hg19
seq_info <- seqinfo(headers)
print(seq_info)

# we save keep_samples to continue working with it in PLINK
write.table(
  data.frame(FID = keep_samples, IID = keep_samples),
  file = "keep_samples.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)







# Checking how similar multiple samples for 1 subject are
# Results: Same subject samples are about 95% similar, a random set of 100 pairs
# of unrelated samples is around 57% similar.
# tehrefore the technical replicates are similar enough to just take one

# 1. Identify subjects with multiple samples
subject_counts <- table(sample_data$dbGaP_Subject_ID[sample_data$dbGaP_Subject_ID %in% csv_headers])
multi_sample_subjects <- names(subject_counts[subject_counts > 1])

cat("Number of subjects with multiple samples:", length(multi_sample_subjects), "\n")

# 2. Check concordance for one or more of these duplicates
# We'll pull a small chunk of genotypes (e.g., first 5000 SNPs) to verify similarity

# prepare results container
concordance_results <- data.frame(
  subject = character(),
  sample1 = character(),
  sample2 = character(),
  concordance = numeric(),
  stringsAsFactors = FALSE
)

for (test_subject in multi_sample_subjects) {
  test_samples <- sample_data$SAMPLE_ID[sample_data$dbGaP_Subject_ID == test_subject]
  
  # read genotypes for all samples of this subject
  param_check <- ScanVcfParam(samples = test_samples, which = GRanges("1", IRanges(1, 10e6))) # adjust range as needed
  vcf_check <- readVcf(vcf_bgz, genome = "hg19", param = param_check)
  genotypes <- geno(vcf_check)$GT
  
  if (is.null(genotypes) || ncol(genotypes) < 2) {
    cat("No or insufficient genotype columns for subject", test_subject, "\n")
    next
  }
  
  sample_names <- colnames(genotypes)
  # compute pairwise concordances
  pairs_mat <- combn(seq_along(sample_names), 2)
  for (k in seq_len(ncol(pairs_mat))) {
    i <- pairs_mat[1, k]
    j <- pairs_mat[2, k]
    conc <- sum(genotypes[, i] == genotypes[, j], na.rm = TRUE) / nrow(genotypes)
    concordance_results <- rbind(
      concordance_results,
      data.frame(
        subject = test_subject,
        sample1 = sample_names[i],
        sample2 = sample_names[j],
        concordance = conc,
        stringsAsFactors = FALSE
      )
    )
  }
}

# print results (formatted as percentages)
if (nrow(concordance_results) == 0) {
  cat("No pairwise concordance results computed.\n")
} else {
  concordance_results$percent <- round(concordance_results$concordance * 100, 2)
  print(concordance_results[, c("subject", "sample1", "sample2", "percent")])
}

# Strategy: Keep only the first sample ID for each subject to match the 1:1 CSV ratio
final_mapping <- sample_data[sample_data$dbGaP_Subject_ID %in% csv_headers, ]
# Keep only the first occurrence of each Subject ID
final_mapping_unique <- final_mapping[!duplicated(final_mapping$dbGaP_Subject_ID), ]

keep_samples <- intersect(final_mapping_unique$SAMPLE_ID, samples)
cat("Final sample count for VCF:", length(keep_samples), "\n")

# Estimate concordance among unrelated samples via random sampling

set.seed(42)

# define how many random pairs
n_pairs <- 100

# map sample -> subject for filtering
sample_to_subject <- setNames(sample_data$dbGaP_Subject_ID, sample_data$SAMPLE_ID)

# keep only samples that made it into VCF
valid_samples <- intersect(samples, sample_data$SAMPLE_ID)

unrelated_results <- data.frame(
  sample1 = character(),
  sample2 = character(),
  concordance = numeric(),
  stringsAsFactors = FALSE
)

pair_count <- 0

while (pair_count < n_pairs) {
  pair <- sample(valid_samples, 2, replace = FALSE)
  
  # skip if same subject
  if (sample_to_subject[pair[1]] == sample_to_subject[pair[2]]) next
  
  # read genotypes for this pair
  param_check <- ScanVcfParam(samples = pair, which = GRanges("1", IRanges(1, 10e6)))
  vcf_check <- readVcf(vcf_bgz, genome = "hg19", param = param_check)
  genotypes <- geno(vcf_check)$GT
  
  if (is.null(genotypes) || ncol(genotypes) < 2) next
  
  conc <- sum(genotypes[,1] == genotypes[,2], na.rm = TRUE) / nrow(genotypes)
  
  unrelated_results <- rbind(
    unrelated_results,
    data.frame(
      sample1 = pair[1],
      sample2 = pair[2],
      concordance = conc,
      stringsAsFactors = FALSE
    )
  )
  
  pair_count <- pair_count + 1
}

# summarize
unrelated_results$percent <- round(unrelated_results$concordance * 100, 2)

cat("Unrelated sample concordance summary:\n")
print(summary(unrelated_results$percent))






# other thing to check

library(vcfR)
vcf <- read.vcfR(vcf_bgz)          # read VCF into vcfR object
gt_mat <- extract.gt(vcf, element = "GT")  # matrix: variants x samples

# sample -> subject mapping (assumes 'sample_data' has SAMPLE_ID and dbGaP_Subject_ID)
map <- setNames(sample_data$dbGaP_SubJECT_ID, sample_data$SAMPLE_ID) # adjust column names if needed

# keep only samples present in the VCF
vcf_samples <- colnames(gt_mat)
map <- map[names(map) %in% vcf_samples]

# identify duplicated subjects
subject_to_samples <- split(names(map), map)         # list: subject -> sample IDs
dups <- subject_to_samples[sapply(subject_to_samples, length) > 1]

# helper: normalize missing and separators
is_missing <- function(x) grepl("^\\.|\\./\\.|\\.\\|\\.$", x) # crude but practical
normalize_gt <- function(x) gsub("\\|", "/", x)

discord_summary <- lapply(names(dups), function(subj) {
  sids <- dups[[subj]]
  g1 <- normalize_gt(gt_mat[, sids[1]])
  g2 <- normalize_gt(gt_mat[, sids[2]])
  nonmiss <- !(is_missing(g1) | is_missing(g2))
  total_comp <- sum(nonmiss)
  if(total_comp == 0) {
    discord_rate <- NA
  } else {
    discord_rate <- mean(g1[nonmiss] != g2[nonmiss])
  }
  callrate1 <- mean(!is_missing(g1))
  callrate2 <- mean(!is_missing(g2))
  list(subject = subj,
       sample1 = sids[1], sample2 = sids[2],
       discord_rate = discord_rate,
       total_comp = total_comp,
       callrate1 = callrate1, callrate2 = callrate2)
})
discord_df <- do.call(rbind, lapply(discord_summary, as.data.frame))
print(discord_df)

