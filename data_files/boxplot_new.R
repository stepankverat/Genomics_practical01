knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  error = TRUE
)
# Load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# Load and reshape expression data
expression_matrix <- read.csv("your_data.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("your_metadata.csv")  # must contain 'sample_id' and 'grouping'


expression_matrix$gene_id <- rownames(expression_matrix)
expression_long <- pivot_longer(expression_matrix, -gene_id, names_to = "sample_id", values_to = "TPM")
expression_annotated <- expression_long %>% left_join(metadata, by = "sample_id")
expression_annotated <- expression_annotated %>% mutate(log_expression = log2(TPM + 1))

# Boxplot to visualize gene expression per sample group
ggplot(expression_annotated, aes(x = grouping, y = log_expression, fill = grouping)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression by Grouping",
       x = "Grouping",
       y = "log2(expression + 1)")