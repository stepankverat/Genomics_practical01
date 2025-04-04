Heatmap Visualization sheet

2025-03-25

Would you like to do a heatmap? Why?

Which input data would you use for heatmap and Why? Should you use all
genes or only a subset?

[]{#anchor}Load required libraries
==================================

*\# Install required packages if not already installed*\
**if** (**!requireNamespace**(\"BiocManager\", quietly = TRUE))
**install.packages**(\"BiocManager\")\
**if** (**!requireNamespace**(\"pheatmap\", quietly = TRUE))
**install.packages**(\"pheatmap\")\
**if** (**!requireNamespace**(\"RColorBrewer\", quietly = TRUE))
**install.packages**(\"RColorBrewer\")\
\
*\# Load libraries*\
**library**(pheatmap) *\# for creating heatmaps*\
**library**(RColorBrewer) *\# for color palettes*\
\
*\# Check your current directory*\
**getwd**()

*\# Set your working directory*\
**setwd**(\"J:/Dropbox/Prague/Teaching/Course genomics of
speciation/2024-2025/RNA-Seq practical/Heatmap visualization\")

[]{#anchor}Import data
======================

*\# Load e*[]{#anchor}*xpression data*\
tpm\_matrix \<- **read.csv**(\"your\_read\_counts.csv\", row.names = 1)

*\# (Optional) Filter genes for better visualization --- e.g., top
variable genes*\
top\_var\_genes \<- **apply**(tpm\_matrix, 1, var)

top\_genes \<- **names**(**sort**(top\_var\_genes, decreasing =
TRUE))\[1**:**50\] *\# top 50 variable genes*

mat\_log \<- **log2**(tpm\_matrix\[top\_genes, \] **+** 1) *\#
log-transform for better contrast*

[]{#anchor-1}Create a basic Heatmap
===================================

**pheatmap**(mat\_log,\
scale = \"row\", *\# normalize expression across each gene*\
clustering\_method = \"ward.D2\", *\# method for hierarchical
clustering*\
show\_rownames = TRUE, *\# show gene names*\
show\_colnames = TRUE, *\# show sample names*\
color = **colorRampPalette**(**rev**(**brewer.pal**(11,
\"RdBu\")))(50),\
main = \"Heatmap of Top Variable Genes\")

[]{#anchor-1}(Optional) Add column annotations by grouping onto the Heatmap
===========================================================================

Again, 'grouping' column in the metadata represents the experimental
condition you want to use to check if samples show similar expression
profile or not.

*\# Load sample metadata*\
metadata\_subset \<- **read.csv**(\"your\_sample\_metadata.csv\",
row.names = 1)

*\# Create annotation data frame for column coloring*\
annotation\_df \<- **data.frame**(grouping =
metadata\_subset**\$**grouping)

**rownames**(annotation\_df) \<- **rownames**(metadata\_subset) *\#
sample IDs should match column names in mat\_log*

*\# Heatmap with sample annotations*\
**pheatmap**(mat\_log,\
scale = \"row\",\
annotation\_col = annotation\_df,\
clustering\_method = \"ward.D2\",\
show\_rownames = TRUE,\
show\_colnames = TRUE,\
color = **colorRampPalette**(**rev**(**brewer.pal**(11,
\"RdBu\")))(50),\
main = \"Heatmap with grouping annotation\")

Do you observe any noticeable differences in gene expression patterns
between sample groupings?
