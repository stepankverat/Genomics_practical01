setwd("/home/verozvest/Documents/Genomics_of_speciation/genomics_practical01/data_files")

sample_info = read.csv("sample_info.csv") #metadata
read_counts = read.csv("read_counts.csv", check.names = FALSE)
tpm_vals = read.csv("tpm_values.csv",check.names = FALSE)
imprinted_genes = read.csv("imprinted_genes.csv")

#imprintovane geny rozdelene podle druhu a pohlavi
mat_peru = read.csv("mat_peru.csv")
pat_peru = read.csv("pat_peru.csv")
mat_chils = read.csv("mat_chils.csv")
pat_chils = read.csv("pat_chil.csv")

#rozdeleni vzorku podle toho jestli jsou to hybridi nebo ne (cross_type)
inter_samples <- sample_info[sample_info$cross_type == "inter", ]
intra_samples <- sample_info[sample_info$cross_type == "intra", ]

#jmena impritovanych genu
imprinted_names <- imprinted_genes$gene.ID
mat_peru_names <- mat_peru$gene.ID
pat_peru_names <- pat_peru$gene.ID
mat_chils_names <- mat_chils$gene.ID
pat_chils_names <- pat_chils$gene.ID

#read_counts pro hybridy PC a CP
inter_sample_names <- inter_samples$sample
inter_columns <- c("gene", inter_sample_names)
inter_read_counts <- read_counts[, colnames(read_counts) %in% inter_columns]

#read counts pro CC a PP
intra_sample_names <- intra_samples$sample
intra_columns <- c("gene", intra_sample_names)
intra_read_counts <- read_counts[, colnames(read_counts) %in% intra_columns]

#tpm rozdeleno podle cross_type
inter_tpm <- tpm_vals[, colnames(tpm_vals) %in% inter_columns]
intra_tpm <- tpm_vals[, colnames(tpm_vals) %in% intra_columns]

#read_counts a tpm dale filtrovane pro imprinted geny
imprinted_read_counts_inter <- inter_read_counts[inter_read_counts$gene %in% imprinted_names, ]
imprinted_read_counts_intra <- intra_read_counts[intra_read_counts$gene %in% imprinted_names, ]

imprinted_tpm_inter <- inter_tpm[inter_tpm$gene %in% imprinted_names, ]
imprinted_tpm_intra <- intra_tpm[intra_tpm$gene %in% imprinted_names, ]
imprinted_tpm <- tpm_vals[tpm_vals$gene %in% imprinted_names, ]

write.csv(imprinted_tpm, "imprinted_tpm.csv", row.names = FALSE)

write.csv(imprinted_read_counts, "imprinted_inter_read_counts.csv", row.names = FALSE)
write.csv(inter_read_counts, "inter_read_counts.csv", row.names = FALSE)

write.csv(imprinted_tpm_inter, "imprinted_tpm_inter.csv", row.names = FALSE)
write.csv(inter_tpm, "inter_tpm.csv", row.names = FALSE)

#group1 = PC = matka peru, otec chille
group1_samples <- c("A9", "A21", "A33")
group1_columns <- c("gene", group1_samples)
group1_table <- imprinted_read_counts_inter[, colnames(imprinted_read_counts_inter) %in% group1_columns]
tpm_peru <- imprinted_tpm_inter[, colnames(imprinted_tpm_inter) %in% group1_columns]
tpm_peru_M <- tpm_peru[tpm_peru$gene %in% mat_per_names, ]
write.csv(tpm_peru_M, "tpm_peru_M.csv", row.names = FALSE)
tpm_chils_P  <- tpm_peru[tpm_peru$gene %in% pat_chils_names, ]
write.csv(tpm_chils_P, "tpm_chils_P.csv", row.names = FALSE)

#goup2 = CP = matka chille, otec peru
group2_samples <- c("A5", "A17", "A29")
group2_columns <- c("gene", group2_samples)
group2_table <- imprinted_read_counts_inter[, colnames(imprinted_read_counts_inter) %in% group2_columns]
tpm_chille <- imprinted_tpm_inter[, colnames(imprinted_tpm_inter) %in% group2_columns]
tpm_peru_P <- tpm_chille[tpm_chille$gene %in% pat_peru_names, ]
write.csv(tpm_peru_P, "tpm_peru_P.csv", row.names = FALSE)
tpm_chils_M  <- tpm_chille[tpm_chille$gene %in% mat_chils_names, ]
write.csv(tpm_chils_M, "tpm_chils_M.csv", row.names = FALSE)

#group3 = PP = oba rodice peru
group3_PP_samples <- c("A3", "A11", "A15","A23","A27","A35")
group3_PP_columns <- c("gene", group3_PP_samples)
group3_PP_table <- imprinted_read_counts_intra[, colnames(imprinted_read_counts_intra) %in% group3_PP_columns]
tpm_PP <- imprinted_tpm_intra[, colnames(imprinted_tpm_intra) %in% group3_PP_columns]

#group4 = CC = oba rodice chille
group4_CC_samples <- c("A1", "A10", "A13","A22","A25","A34")
group4_CC_columns <- c("gene", group4_CC_samples)
group4_CC_table <- imprinted_read_counts_intra[, colnames(imprinted_read_counts_intra) %in% group4_CC_columns]
tpm_CC <- imprinted_tpm_intra[, colnames(imprinted_tpm_intra) %in% group4_CC_columns]

write.csv(group1_table, "imprinted_P_mat_samples.csv", row.names = FALSE)
write.csv(group2_table, "imprinted_C_mat_samples.csv", row.names = FALSE)

write.csv(tpm_peru, "imprinted_P_tpm.csv", row.names = FALSE)
write.csv(tpm_chille, "imprinted_C_tpm.csv", row.names = FALSE)

write.csv(group3_PP_table, "imprinted_PP_samples.csv", row.names = FALSE)
write.csv(group4_CC_table, "imprinted_CC_samples.csv", row.names = FALSE)

write.csv(tpm_PP, "PP_tpm.csv", row.names = FALSE)
write.csv(tpm_CC, "CC_tpm.csv", row.names = FALSE)
