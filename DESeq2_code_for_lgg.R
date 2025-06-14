library(tidyverse)
library(DESeq2)


#read count matrix
count_matrix <- read_csv("D:/Thesis_project/Thesis/all_lgg_with_normal.csv")


sample_metadata <- read_csv("D:/Thesis_project/Thesis/merged_meta.csv")

#make sample_id rowname of sample_metadata
rownames(sample_metadata) <- sample_metadata$sample_id

view(sample_metadata)
view(count_matrix)



#show what did not match
not_matching <- colnames(count_matrix)[!colnames(count_matrix) %in% rownames(sample_metadata)]
not_matching #"...1" not marching is expected due to the first column being gene names


# Step 3: Remove the gene column from count matrix
gene_names <- count_matrix[[1]]  # Save gene names if needed later


view(gene_names)

#check if there any duplicates in gene_names
any(duplicated(gene_names))

#drop first column of count_matrix
count_matrix <- count_matrix[, -1]


#set rownames to gene names otherwise deseq2 won't show gene name
rownames(count_matrix) <- gene_names
view(count_matrix)


# making sure the row names in colData matches to column names in counts_data
all(colnames(count_matrix) %in% rownames(sample_metadata)) 
#check if they are in same order
all(colnames(count_matrix) == rownames(sample_metadata)) 



# Make sure metadata factors are correctly set
sample_metadata$Target <- factor(sample_metadata$Target)





# Convert to matrix before passing to DESeq2
count_matrix <- as.matrix(count_matrix)

# Check again that the row names are gene names
head(rownames(count_matrix))  # should now show actual gene names




# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                              colData = sample_metadata,
                              design = ~ Target)
# Run DESeq2
dds <- DESeq(dds)

# Get results based on Target
res <- results(dds)

view(res)

#seperate differentially expressed genes
res_df <- as.data.frame(res)
res_df <- res_df %>%
  rownames_to_column(var = "gene_name") %>%
  arrange(padj)
# Save results to CSV
write_csv(res_df, "DESeq2_results/DESeq2_results_lgg_14June.csv")

# Filter significant DEGs (adjust thresholds as needed)
significant_DEGs <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj)
# Save significant DEGs to CSV
write_csv(significant_DEGs, "DESeq2_results/DESeq2_significant_DEGs_lgg_14june.csv")


#show how many upregulated and downregulated genes
upregulated_count <- sum(significant_DEGs$log2FoldChange > 0)
downregulated_count <- sum(significant_DEGs$log2FoldChange < 0)
cat("Upregulated genes:", upregulated_count, "\n")
cat("Downregulated genes:", downregulated_count, "\n")






