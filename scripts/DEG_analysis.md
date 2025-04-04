DEG analysis sheet

2025-03-25

[]{#anchor}Load required libraries
==================================

*\# Install required packages if not already installed*\
**if** (**!requireNamespace**(\"DESeq2\", quietly = TRUE)) {\
**install.packages**(\"BiocManager\")\
BiocManager**::install**(\"DESeq2\")\
}\
\
*\# Load libraries*\
**library**(DESeq2)\
\
*\# Check your current directory*\
**getwd**()

*\# Set your working directory*\
**setwd**(\"J:/Dropbox/Prague/Teaching/Course genomics of
speciation/2024-2025/RNA-Seq practical/DEG analysis\")

[]{#anchor}Import data
======================

Load the raw gene expression data (genes as rows, samples as columns) as
well as the sample information (experimental design, e.g.,
control/treatment).

*\# Load the counts data and metadata*\
counts \<- **read.csv**(\"your\_read\_counts.csv\", row.names = 1)

metadata \<- **read.csv**(\"your\_sample\_metadata.csv\", row.names = 1)

*\# Check the structure of each input files*\
**head**(counts)

**head**(metadata)

Think about your experimental design. Which experimental groups would
you like to compare? Based on what we learnt about contrasts in DESeq2,
would you apply the approach outlined in Chapter I (compare only two
groups) or Chapter II (compare two or more groups using a special
function)?

### []{#anchor}*Chapter I: Compare only two groups*

Use this approach when you're comparing exactly **two groups**, such as:

\- Control vs Treatment

\- Wildtype vs Mutant

\- Genotype A vs Genotype B

Here, you use the formula design = \~grouping, which tells DESeq2 to
model gene expression based on the **levels in the 'grouping' column**
of your metadata. So, 'grouping' refers to the column in your metadata
file that defines experimental groups. \# \~grouping includes an
intercept (baseline group)

By setting the levels manually (e.g., c(\"control\", \"treatment\")),
you're defining which is the **baseline (control)** and which is the
**comparison (treatment)** group.

[]{#anchor}Create a DESeq2 dataset and run the analysis
=======================================================

To perform differential expression analysis, DESeq2 needs both the **raw
gene counts** and the **sample information (metadata)** in a specific
format --- this is called a **DESeqDataSet**.

The DESeq2 dataset combines:

\- The **counts matrix** (genes x samples) --- tells DESeq2 how many
reads mapped to each gene.

\- The **metadata** (samples x variables) --- tells DESeq2 which samples
belong to which experimental groups (e.g., treatment vs control).

\- A **design formula** --- defines which variable(s) DESeq2 should use
to compare groups.

*\# Create a DESeq2 object*\
dds \<- **DESeqDataSetFromMatrix**(countData = counts,\
colData = metadata,\
design = **\~**grouping) *\# Make sure to replace \'grouping\' with the
column you\'re using for comparison in your metadata*

*\# Make sure grouping is a factor*\
dds**\$**grouping \<- **factor**(dds**\$**grouping, levels=**c**(\"level
1\", \"level 2\"))

*\# (OPTIONAL) Remove genes with very low counts across all samples*\
*\# Keep genes where the total count across all samples is at least 10*\
dds \<- dds\[**rowSums**(**counts**(dds)) **\>=** 10, \]

*\# Run the differential expression analysis*\
dds \<- **DESeq**(dds)

[]{#anchor-1}View and filter results
====================================

*\# Perform your desired contrast *\
results\_table \<- **results**(dds, contrast = **c**(\"grouping\",
\"level 1\", \"level 2\"))

*\# In a DESeq2 model, R compares two groups:*\
*\# - One group is used as the \"baseline\" (this is the intercept)*\
*\# - The difference between the other group and the baseline is the
\"coefficient\"*\
*\# This tells us how much more or less a gene is expressed in one group
compared to the other.*\
\
*\# Remove genes with missing values*\
results\_table \<- **na.omit**(results\_table)

*\# Filter for significant genes (adjusted p-value \< 0.05)*\
significant\_genes \<- results\_table\[results\_table**\$**padj **\<**
0.05, \]

*\# Check how many genes are significant*\
**nrow**(significant\_genes)

How many genes passed the significance threshold?

[]{#anchor-1}Identify Up- and Down-regulated genes
==================================================

*\# Up-regulated in inter (log2 Fold Change \> 1)*\
up\_genes \<- significant\_genes\[significant\_genes**\$**log2FoldChange
**\>** 1, \]

*\# Down-regulated in inter (log2 Fold Change \< -1)*\
down\_genes \<-
significant\_genes\[significant\_genes**\$**log2FoldChange **\<**
**-**1, \]

Out of all the significantly differentially expressed genes, how many
were upregulated and how many were downregulated?

[]{#anchor-1}Save the results
=============================

**write.csv**(**as.data.frame**(up\_genes), \"your\_output\_up.csv\")

**write.csv**(**as.data.frame**(down\_genes),
\"your\_output\_down.csv\")

### []{#anchor-1}*Chapter II: Compare two or more groups using contraster()*

Use this approach when:

\- You want to compare **a group vs. the average** of multiple other
groups

\- You want to customize contrasts beyond simple pairwise comparisons

\- Your experimental design involves more than two groups[]{#anchor-2}

Here, the formula is design = \~0 + grouping, which removes the default
intercept and lets you manually build contrasts using the contraster()
function.

The grouping column in your metadata should still define the biological
groups (e.g., genotypes, treatments), but now you have more flexibility.

For example:

group1 = list(c("grouping", "treatment A", "treatment B")) group2 =
list(c("grouping", "control"))

This compares the average of treatment A and treatment B to control.

[]{#anchor-3}Create a DESeq2 dataset and run the analysis
=========================================================

*\# Create a DESeq2 object*\
dds \<- **DESeqDataSetFromMatrix**(countData = counts,\
colData = metadata,\
design = **\~**0**+**grouping)\
*\# \~0 + grouping removes intercept and estimates each group
separately*

*\# Make sure to replace \'grouping\' with the column you\'re using for
comparison in your metadata*\
\
*\# Make sure grouping is a factor*\
dds**\$**grouping \<- **factor**(dds**\$**grouping, levels=**c**(\"level
1\", \"level 2\", \"level 3\"\...))\
\
*\# (OPTIONAL) Remove genes with very low counts across all samples*\
*\# Keep genes where the total count across all samples is at least 10*\
dds \<- dds\[**rowSums**(**counts**(dds)) **\>=** 10, \]\
\
*\# Run the differential expression analysis*\
dds \<- **DESeq**(dds)

[]{#anchor-4}Define the contraster function
===========================================

contraster \<- **function**(dds, *\# should contain colData and design*\
group1, *\# list of character vectors each with 2 or more items *\
group2, *\# list of character vectors each with 2 or more items*\
weighted = F\
){\
\
\
mod\_mat \<- **model.matrix**(**design**(dds), **colData**(dds))\
\
grp1\_rows \<- **list**()\
grp2\_rows \<- **list**()\
\
\
**for**(i **in** 1**:length**(group1)){\
\
grp1\_rows\[\[i\]\] \<- **colData**(dds)\[\[group1\[\[i\]\]\[1\]\]\]
**%in%** group1\[\[i\]\]\[2**:length**(group1\[\[i\]\])\]\
\
}\
\
\
**for**(i **in** 1**:length**(group2)){\
\
grp2\_rows\[\[i\]\] \<- **colData**(dds)\[\[group2\[\[i\]\]\[1\]\]\]
**%in%** group2\[\[i\]\]\[2**:length**(group2\[\[i\]\])\]\
\
}\
\
grp1\_rows \<- **Reduce**(**function**(x, y) x **&** y, grp1\_rows)\
grp2\_rows \<- **Reduce**(**function**(x, y) x **&** y, grp2\_rows)\
\
mod\_mat1 \<- mod\_mat\[grp1\_rows, ,drop=F\]\
mod\_mat2 \<- mod\_mat\[grp2\_rows, ,drop=F\]\
\
**if**(**!**weighted){\
\
mod\_mat1 \<- mod\_mat1\[**!duplicated**(mod\_mat1),,drop=F\]\
mod\_mat2 \<- mod\_mat2\[**!duplicated**(mod\_mat2),,drop=F\]\
\
}\
\
**return**(**colMeans**(mod\_mat1)**-colMeans**(mod\_mat2))\
\
\
}

[]{#anchor-4}View and filter results
====================================

*\# Perform your desired contrast *\
results\_table1 \<- **results**(dds,\
contrast = **contraster**(dds,\
group1 = **list**(**c**(\"grouping\", \"level 1\")),\
group2 = **list**(**c**(\"grouping\", \"level 3\", \"level 4\"))))

*\# Check the intercept and coefficients\**\
**contraster**(dds,\
group1 = **list**(**c**(\"grouping\", \"level 1\")),\
group2 = **list**(**c**(\"grouping\", \"level 3\", \"level 4\")))

*\# Remove genes with missing values*\
results\_table1 \<- **na.omit**(results\_table1)

*\# Filter for significant genes (adjusted p-value \< 0.05)*\
significant\_genes1 \<- results\_table1\[results\_table1**\$**padj
**\<** 0.05, \]

*\# Check how many genes are significant*\
**nrow**(significant\_genes1)

Simple Analogy for Intercepts and Coefficients:

Imagine you're comparing how tall people are in two different
classrooms: Class A and Class B.

-   The intercept is like the average height of the first group you're
    using as a baseline --- for example, Class A.
-   The coefficient is how much taller or shorter people in Class B are
    compared to Class A.

Here, we're doing the same thing --- just not with height, but with gene
expression.

[]{#anchor-4}Identify Up- and Down-regulated genes
==================================================

*\# Up-regulated in inter (log2 Fold Change \> 1)*\
up\_genes1 \<-
significant\_genes1\[significant\_genes1**\$**log2FoldChange **\>** 1,
\]

*\# Down-regulated in inter (log2 Fold Change \< -1)*\
down\_genes1 \<-
significant\_genes1\[significant\_genes1**\$**log2FoldChange **\<**
**-**1, \]

[]{#anchor-4}Save the results
=============================

*\# Save results to CSV*\
**write.csv**(**as.data.frame**(up\_genes1), \"your\_output\_up.csv\")

**write.csv**(**as.data.frame**(down\_genes1),
\"your\_output\_down.csv\")
