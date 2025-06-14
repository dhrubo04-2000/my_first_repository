# Load required libraries
library(DESeq2)
library(readr)
library(dplyr)
library(tidyverse)

setwd('D:/Thesis_project/Thesis')
# Read the CSV
data1 <- read_csv("filtered_predictive_genes_lgg.csv")

# Step 2: Set gene names as row names
count_data <- data1 %>% column_to_rownames(var = "gene_name")


# Step 3: Create dummy metadata for samples
sample_ids <- colnames(count_data)
col_data <- data.frame(row.names = sample_ids,
                       condition = rep("A", length(sample_ids)))  # dummy condition



# Ensure counts are integers
count_data[] <- lapply(count_data, function(x) as.integer(round(x)))

# Step 4: Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~1)  # design not needed for blind VST


# Step 5: Apply VST
vst_data <- varianceStabilizingTransformation(dds, blind = TRUE)

vst_matrix <- assay(vst_data)
vst_df <- as.data.frame(t(vst_matrix))
vst_df <- tibble::rownames_to_column(vst_df, var = "sample_id")

view(vst_df)

write_csv(vst_df, "vst_filtered_predictive_genes_lgg.csv")










