setwd("/home/verozvest/Documents/Genomics_of_speciation/genomics_practical01/data_files")

read_counts <- read.csv("read_counts.csv", check.names = FALSE)
gene_lengths <- read.csv("gene_lengths.csv")

rownames(read_counts) <- read_counts$gene
read_counts <- read_counts[, -1]

gene_length_vector <- gene_lengths$length
names(gene_length_vector) <- gene_lengths$gene_id


calculate_tpm <- function(counts, gene_lengths) {
  rpk <- counts / (gene_lengths / 1000)
  scaling_factors <- colSums(rpk) / 1e6
  tpm <- sweep(rpk, 2, scaling_factors, "/")
  return(tpm)
}

tpm_values <- calculate_tpm(as.matrix(read_counts), gene_length_vector)

tpm_df <- as.data.frame(tpm_values)
tpm_df$gene <- rownames(tpm_df)
tpm_df <- tpm_df[, c("gene", colnames(tpm_df)[colnames(tpm_df) != "gene"])]

write.csv(tpm_df, "tpm_values.csv", row.names = FALSE)
