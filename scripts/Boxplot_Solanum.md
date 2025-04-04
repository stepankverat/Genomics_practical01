Boxplot visualization

2025-03-25

[]{#anchor}Load required libraries
==================================

*\# Load libraries*\
**library**(dplyr)\
**library**(tidyr)\
**library**(ggplot2)

What questions can boxplots help you answer in gene expression analysis?
What kind of data can you use as input[]{#anchor}?

*\# Load and reshape expression data*\
expression\_matrix \<- **read.csv**(\"your\_data.csv\", row.names = 1,
check.names = FALSE)

metadata \<- **read.csv**(\"your\_metadata.csv\") *\# must contain
\'sample\_id\' and \'grouping\'*

expression\_matrix**\$**gene\_id \<- **rownames**(expression\_matrix)

expression\_long \<- **pivot\_longer**(expression\_matrix,
**-**gene\_id, names\_to =

expression\_annotated \<- expression\_long **%\>%**
**left\_join**(metadata, by = \"sample\_id\")

expression\_annotated \<- expression\_annotated **%\>%**
**mutate**(log\_expression = **log2**(expression **+** 1))

*\# Boxplot to visualize gene expression per sample group*\
**ggplot**(expression\_annotated, **aes**(x = grouping, y =
log\_expression, fill = grouping)) **+**\
**geom\_boxplot**(outlier.size = 0.5, alpha = 0.7) **+**\
**theme\_bw**() **+**\
**theme**(axis.text.x = **element\_text**(angle = 45, hjust = 1)) **+**\
**labs**(title = \"Gene Expression by Grouping\",\
x = \"Grouping\",\
y = \"log2(expression + 1)\")
