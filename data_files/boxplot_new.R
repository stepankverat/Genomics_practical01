knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  error = TRUE
)

library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/home/verozvest/Documents/genomics_practical01/data_files")
expression_matrix <- read.csv("imprinted_tpm.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata.csv") 

expression_matrix$gene_id <- rownames(expression_matrix)
expression_long <- pivot_longer(expression_matrix, -gene_id, names_to = "sample", values_to = "TPM")
expression_annotated <- expression_long %>% left_join(metadata, by = "sample")
expression_annotated <- expression_annotated %>% mutate(log_expression = log2(TPM + 1))


ggplot(expression_annotated, aes(x = cross_spe, y = log_expression, fill = cross_spe)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Gene Expression by Grouping",
       x = "Grouping",
       y = "log2(expression + 1)")