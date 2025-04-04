[]{#anchor}PCA visualization sheet

2025-03-25

Why are we using PCA in RNA-seq analysis? What can it tell us about the
samples? What kind of variation is PCA capturing --- technical or
biological?

Which input data would you use for PCA and Why? Should you use all genes
for PCA? Why or why not? These questions will help you use to understand
script below.

[]{#anchor-1}Get started
========================

*\# Check your current directory*\
**getwd**()

*\# Set your working directory*\
**setwd**(\"J:/Dropbox/Prague/Teaching/Course genomics of
speciation/2024-2025/RNA-Seq practical/PCA visualization\")

[]{#anchor-1}Import data
========================

*\# Load expression data and sample metadata (contains sample
information)*\
expr\_matrix \<- **read.csv**(\"your\_read\_counts.csv\", row.names = 1)

metadata \<- **read.csv**(\"your\_sample\_metadata.csv\", row.names = 1)

*\# Check the structure of the imported files*\
**head**(expr\_matrix)

**head**(metadata)

*\# Check that the samples in both files match*\
**all**(**colnames**(expr\_matrix) **==** **rownames**(metadata)) *\#
should return TRUE*

What kind of information about the experimental design can you find in
the metadata file?

[]{#anchor-1}Prepare files before running PCA
=============================================

*\# Add 1 to ensure that the smallest value is log2(1) = 0*\
log\_expr \<- **log2**(expr\_matrix **+** 1)

*\# Calculate the variance in expression of each gene across all
samples*\
var\_genes \<- **apply**(log\_expr, 1, var)

*\# Filter to keep only genes with non-zero variance*\
log\_expr\_filtered \<- log\_expr\[var\_genes **\>** 0, \]

Why did we apply a log transformation to expression values before
running PCA? Why did we remove genes with zero variance across all
samples?

[]{#anchor-1}Running PCA
========================

To perform PCA in R, use the prcomp() function with a matrix where rows
are samples and columns are variables (e.g., genes).

*\# Transpose the matrix so rows = samples, columns = genes*\
*\# Center and scale genes so each contributes equally to PCA*\
pca \<- **prcomp**(**t**(log\_expr\_filtered), scale. = TRUE)

*\# View a summary of the PCA results*\
**summary**(pca)

scale. = TRUE standardizes each gene (z-score), ensuring that PCA isn't
dominated by highly expressed genes.

[]{#anchor-1}Visualizing PCA
============================

*\# Basic PCA scatter plot using base R*\
**plot**(pca**\$**x\[,1\], pca**\$**x\[,2\],\
xlab = \"PC1\", ylab = \"PC2\",\
main = \"PCA Plot\",\
pch = 19, col = \"blue\")

*\# Calculate % variance explained*\
percent\_var \<- **round**(**summary**(pca)**\$**importance\[2,
1**:**2\] **\*** 100, 1)

*\# Extract PC scores for samples*\
pca\_scores \<- pca**\$**x\[, 1**:**2\]

*\# Get grouping from metadata (make sure order matches rownames of
pca\$x)*\
grouping \<- **as.factor**(metadata\[**rownames**(pca\_scores),
\"grouping\"\]) *\# \'grouping\' column in the metadata represents the
experimental condition you want to use to check if samples cluster
together or differ in the PCA*

*\# Assign color values based on grouping factor levels*\
colors \<- **as.factor**(grouping)

*\# Plot PCA*\
**plot**(pca\_scores,\
col = colors,\
pch = 19,\
xlab = **paste0**(\"PC1 (\", percent\_var\[1\], \"%)\"),\
ylab = **paste0**(\"PC2 (\", percent\_var\[2\], \"%)\"),\
main = \"PCA Plot Colored by Grouping\")

*\# Add legend*\
**legend**(\"topright\",\
legend = **levels**(grouping),\
col = 1**:length**(**levels**(grouping)),\
pch = 19)

Do your PCA results support your experimental design or hypothesis? Why
or why not?
