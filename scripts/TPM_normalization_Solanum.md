TPM normalization sheet

2025-03-25

You have expression data in the form of raw read counts. What should you
do first? Would you use these raw counts directly for downstream
analyses like clustering, PCA, or differential expression? Why or why
not? What factors influence the number of reads mapped to a gene? What
is the purpose of normalization in RNA-seq analysis?

Basically, TPM (Transcripts Per Million) adjusts raw read counts by
accounting for gene length and sequencing depth for a fair comparison of
expression levels across genes and across samples.

[]{#anchor}Load required libraries
==================================

*\# Install dplyr if needed*\
**install.packages**(\"dplyr\")

*\# Load libraries*\
**library**(dplyr)\
\
*\# Check your current directory*\
**getwd**()

*\# Set your working directory*\
**setwd**(\"J:/Dropbox/Prague/Teaching/Course genomics of
speciation/2024-2025/RNA-Seq practical/TPM normalization\")

[]{#anchor}Import raw read counts and gene lengths
==================================================

*\# Read count matrix: rows = genes, columns = samples*\
read\_counts \<- **read.csv**(\"your\_read\_counts.csv\", row.names = 1)

*\# Gene lengths in kilobases (must contain columns: gene\_id, length)*\
gene\_length\_df \<- **read.csv**(\"your\_gene\_lengths.csv\")

[]{#anchor}Filter and match gene IDs
====================================

*\# Find common gene IDs between read counts and gene length file*\
common\_ids \<- **intersect**(**rownames**(read\_counts),
gene\_length\_df**\$**gene\_id)

*\# Subset and align*\
read\_counts \<- read\_counts\[common\_ids, \]

gene\_length\_df \<- gene\_length\_df\[**match**(common\_ids,
gene\_length\_df**\$**gene\_id), \]

*\# Extract lengths as numeric vector (already in kb)*\
lengths\_kb \<- gene\_length\_df**\$**length

[]{#anchor}TPM normalization
============================

Now that you have read counts and gene length, how would you calculate
transcripts per million ? Remember, TPM adjusts for gene length and
sequencing depth. Hint: calculate reads per kilobase for each gene and
then normalize to per million.

*\# Define TPM normalization function*\
calculate\_tpm \<- **function**(counts, lengths\_kb) {\
rpk \<- **sweep**(counts, 1, lengths\_kb, FUN = \"/\") *\# Reads Per
Kilobase*\
tpm \<- **sweep**(rpk, 2, **colSums**(rpk), FUN = \"/\") **\*** 1e6 *\#
Normalize to per million*\
**return**(tpm)\
}\
\
*\# Run TPM normalization*\
tpm\_matrix \<- **calculate\_tpm**(**as.matrix**(read\_counts),
lengths\_kb)

[]{#anchor}Export TPM matrix
============================

**write.csv**(tpm\_matrix, \"your\_output.csv\", quote = FALSE)

You now have a normalized expression matrix[]{#anchor}.
