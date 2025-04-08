setwd("/home/verozvest/Documents/genomics_practical01/data_files")

sample_info = read.csv("sample_info.csv")
read_counts = read.csv("read_counts.csv", check.names = FALSE)
tpm_vals = read.csv("tpm_values.csv",check.names = FALSE)
imprinted_genes = read.csv("imprinted_genes.csv")
inter_samples <- sample_info[sample_info$cross_type == "inter", ]

imprinted_names <- imprinted_genes$gene.ID

inter_sample_names <- inter_samples$sample
inter_columns <- c("gene", inter_sample_names)
inter_read_counts <- read_counts[, colnames(read_counts) %in% inter_columns]

inter_tpm <- tpm_vals[, colnames(tpm_vals) %in% inter_columns]

imprinted_read_counts <- inter_read_counts[inter_read_counts$gene %in% imprinted_names, ]

imprinted_tpm <- inter_tpm[inter_tpm$gene %in% imprinted_names, ]

write.csv(imprinted_read_counts, "imprinted_inter_read_counts.csv", row.names = FALSE)
write.csv(inter_read_counts, "inter_read_counts.csv", row.names = FALSE)

write.csv(imprinted_tpm, "imprinted_tpm.csv", row.names = FALSE)
write.csv(inter_tpm, "inter_tpm.csv", row.names = FALSE)

group1_samples <- c("A9", "A21", "A33")
group1_columns <- c("gene", group1_samples)
group1_table <- imprinted_read_counts[, colnames(imprinted_read_counts) %in% group1_columns]
tpm_peru <- imprinted_tpm[, colnames(imprinted_tpm) %in% group1_columns]

group2_samples <- c("A5", "A17", "A29")
group2_columns <- c("gene", group2_samples)
group2_table <- imprinted_read_counts[, colnames(imprinted_read_counts) %in% group2_columns]
tpm_chille <- imprinted_tpm[, colnames(imprinted_tpm) %in% group2_columns]

write.csv(group1_table, "imprinted_P_mat_samples.csv", row.names = FALSE)
write.csv(group2_table, "imprinted_C_mat_samples.csv", row.names = FALSE)

write.csv(tpm_peru, "imprinted_P_tpm.csv", row.names = FALSE)
write.csv(tpm_chille, "imprinted_C_tpm.csv", row.names = FALSE)

metadata <- data.frame(
  sample_id = c("A5", "A9", "A17", "A21", "A29", "A33"),
  grouping = c("CP", "PC", "CP", "PC", "CP", "PC")
)

# Write metadata file
write.csv(metadata, "metadata.csv", row.names = FALSE)