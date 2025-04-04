GO enrichment sheet

2025-03-27

Before starting, a few questions for you: what does a GO enrichment do?
What can it tell us? By definition, an enrichment means an
over-representation in a group compared to another? What are those two
groups? These questions are important to answer so that you use the
script below properly.

[]{#anchor}Getting started - special installation
=================================================

*\#install clusterProfiler (special way; for the other packages, the
normal way should work)*\
**if** (**!require**(\"BiocManager\", quietly = TRUE))\
**install.packages**(\"BiocManager\")\
\
BiocManager**::install**(\"clusterProfiler\")\
\
*\#load packages, set up the working directory*\
**library**(clusterProfiler)\
**library**(ggplot2)\
**library**(tidyr)\
**library**(dplyr)\
**getwd**()

**setwd**(\"J:/Dropbox/Prague/Teaching/Course genomics of
speciation/2024-2025/RNA-Seq practical/GO enrichment\")

[]{#anchor}Importing the reference GO files
===========================================

*\#The files should look like this: *\
*\#GO:0055114 Solyc02g080590 Biological Process:oxidation-reduction
process*\
*\#GO:0055114 Solyc08g007400 Biological Process:oxidation-reduction
process*\
*\#GO:0055114 Solyc03g031650 Biological Process:oxidation-reduction
process*\
\
SlycBP \<- **read.csv**(\"SlycBP.csv\")\
term2gene\_BP \<- SlycBP\[,**c**(1,2)\] *\#extract first 2 columns (GO
term and gene ID)*

term2name\_BP \<- SlycBP\[,**c**(1,3)\] *\#extract Go term and function*

SlycCC \<- **read.csv**(\"SlycCC.csv\")\
term2gene\_CC \<- SlycCC\[,**c**(1,2)\] *\#extract first 2 columns (GO
term and gene ID)*

term2name\_CC \<- SlycCC\[,**c**(1,3)\] *\#extract Go term and function*

SlycMF \<- **read.csv**(\"SlycMF.csv\")

term2gene\_MF \<- SlycMF\[,**c**(1,2)\] *\#extract first 2 columns (GO
term and gene ID)*

term2name\_MF \<- SlycMF\[,**c**(1,3)\] *\#extract Go term and function*

[]{#anchor}Reading your list of genes
=====================================

*\#the file should look like this (only 1 column with gene names):*\
*\#Solyc03g031650*\
*\#Solyc03g031650*\
*\#Solyc03g031650*\
*\#Solyc03g031650*\
\
candidate\_genes \<- **read.csv**(\"XXXXX.csv\")

candidate\_genes = **as.vector**(candidate\_genes\[,1\])

*\#Check if they intersect*\
**intersect**(term2gene\_BP**\$**locusName, candidate\_genes)

[]{#anchor}Testing for enrichment (PS: between what and what?)
==============================================================

*\#BP*\
candidate\_genes\_BP \<- **enricher**(gene= candidate\_genes,
pvalueCutoff = 0.05, pAdjustMethod = \"fdr\", TERM2GENE= term2gene\_BP,
TERM2NAME = term2name\_BP)

**head**(candidate\_genes\_BP)

candidate\_genes\_BP \<- candidate\_genes\_BP\[,**c**(2,3,4,6,8,9)\]
*\#function, p-value, geneIDs, count*

candidate\_genes\_BP**\$**GeneRatio \<- **paste0**(\"\'\",
candidate\_genes\_BP**\$**GeneRatio) *\# Add an apostrophe before each
GeneRatio value, otherwise it will be read as a date from excel*

**write.csv**(candidate\_genes\_BP, \'candidate\_genes\_GO\_BP.csv\')

This was for Biological Processes. Now you're gonna need to do the same
for Celllular Component and Molecular Function.

[]{#anchor-1}Preparing the file with all enriched terms for plotting
====================================================================

This puts together the enriched GO terms in a csv file that looks like
this:

Category Description GeneRatio p.adjust\
upregulated Molecular Function:DNA-binding transcription factor activity
1 0.023282605

*\#Give a name to your group of candidate genes (as an example here is
\"upregulated\"), add it as column to each of the enriched dfs; if you
have several groups of genes, do that for every group*\
candidate\_genes\_BP**\$**Category \<- \"upregulated\"

candidate\_genes\_CC**\$**Category \<- \"upregulated\"

candidate\_genes\_MF**\$**Category \<- \"upregulated\"

All\_GO \<- **rbind**(candidate\_genes\_BP,\
candidate\_genes\_CC,candidate\_genes\_MF) *\#if you have several groups
of genes, you can put them all together here if you want them plotted
together*

[]{#anchor-2}Plotting to visualize the enriched GO terms
========================================================

This will generate dotplots showing gene ratio and adj pvalue.

GO \<- All\_GO

GO**\$**Description \<- **as.factor**(GO**\$**Description)

*\# Edit our GeneRatio column*\
*\# Remove the leading \' from the strings*\
GO**\$**GeneRatio \<- **gsub**(\"\^\'\", \"\", GO**\$**GeneRatio)

*\# Use separate() to split the column into two new columns*\
\
GO \<- GO **%\>%** **separate**(GeneRatio, into = **c**(\"Numerator\",
\"Denominator\"), sep = \"/\")

*\# Convert the new columns to numeric*\
GO**\$**Numerator \<- **as.numeric**(GO**\$**Numerator)

GO**\$**Denominator \<- **as.numeric**(GO**\$**Denominator)

GO**\$**GeneRatio \<- GO**\$**Numerator**/**GO**\$**Denominator

*\# Remove extra columns*\
GO**\$**Numerator \<- NULL

GO**\$**Denominator \<- NULL

*\#plot all*\
**ggplot**(data = GO, **aes**(x = Category, y = Description,\
color = \`p.adjust\`, size = GeneRatio)) **+**\
**geom\_point**() **+**\
**scale\_color\_gradient**(low = \"red\", high = \"blue\") **+**\
**theme\_bw**() **+**\
**ylab**(\"\") **+**\
**xlab**(\"\") **+**\
**theme\_bw**() **+**\
*\#theme(axis.title.x = element\_text(size = 24), axis.text.x =
element\_text(size = 24)) +*\
**theme**(axis.text.y = **element\_text**(size = 12)) **+**\
*\#theme(legend.title = element\_blank(), legend.text =
element\_text(size = 24)) + *\
*\#theme(axis.text.x = element\_text(angle = 25, vjust = 0.8, hjust =1,
size = 12)) +*\
**theme**(axis.text.x = **element\_text**(vjust = 0.8, hjust =1, size =
12)) **+**\
**scale\_y\_discrete**(limits=rev)

**dev.off**() *\#this closes the graph, so if you want to see the graph,
don\'t run that :)*
