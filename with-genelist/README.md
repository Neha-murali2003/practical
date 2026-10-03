# with-genelist
############################################################
# EXPERIMENT: ORA AND GSEA USING KEGG AND REACTOME
############################################################


############################################################
# STEP 1: INSTALL AND LOAD REQUIRED PACKAGES
############################################################

# Install packages if they are not already installed

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c(
  "clusterProfiler",
  "ReactomePA",
  "DOSE",
  "enrichplot",
  "org.Hs.eg.db",
  "ggplot2"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# Load libraries

library(clusterProfiler)
library(ReactomePA)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)


############################################################
# STEP 2: LOAD THE INPUT CSV FILE
############################################################

# Select the CSV file containing gene symbols and logFC

df <- read.csv(file.choose())

# Display the first few rows

head(df)

# Check column names

colnames(df)

# Check dimensions

dim(df)


############################################################
# STEP 3: DATA PREPROCESSING
############################################################

# Remove missing values

df <- na.omit(df)

# Ensure correct data types

df$gene_symbol <- as.character(df$gene_symbol)
df$logFC <- as.numeric(df$logFC)

# Remove duplicate gene symbols

df <- df[!duplicated(df$gene_symbol), ]

# View cleaned data

head(df)

# Summary of logFC values

summary(df$logFC)


############################################################
# STEP 4: CONVERT GENE SYMBOLS TO ENTREZ IDs
############################################################

# Convert gene symbols into Entrez Gene IDs

gene_conversion <- bitr(
  df$gene_symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Merge converted IDs with original logFC values

merged_data <- merge(
  gene_conversion,
  df,
  by.x = "SYMBOL",
  by.y = "gene_symbol"
)

# Remove duplicate Entrez IDs

merged_data <- merged_data[
  !duplicated(merged_data$ENTREZID),
]

# Display converted data

head(merged_data)

# Check number of successfully mapped genes

nrow(merged_data)


############################################################
# STEP 5: PREPARE INPUTS FOR ORA
############################################################

# Select genes with absolute logFC greater than 1.5

ora_genes <- merged_data$ENTREZID[
  abs(merged_data$logFC) > 1.5
]

# Define the background gene universe

universe_genes <- merged_data$ENTREZID

# Display number of selected genes

length(ora_genes)

# Display selected genes

ora_genes


############################################################
# STEP 6: PREPARE INPUT FOR GSEA
############################################################

# Create a named numeric vector using logFC values

geneList <- merged_data$logFC

# Assign Entrez IDs as names

names(geneList) <- merged_data$ENTREZID

# Sort genes in decreasing order of logFC

geneList <- sort(geneList, decreasing = TRUE)

# View ranked gene list

head(geneList)


############################################################
# STEP 7: KEGG ORA
############################################################

# Perform Over-Representation Analysis using KEGG

kk <- enrichKEGG(
  gene = ora_genes,
  organism = "hsa",
  keyType = "ncbi-geneid",
  universe = universe_genes,
  pvalueCutoff = 0.05
)

# Convert results into a data frame

kk_df <- as.data.frame(kk)

# Display results

head(kk_df)

# Convert Entrez IDs into readable gene symbols

kk_readable <- setReadable(
  kk,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

# Visualize KEGG ORA results

dotplot(kk, showCategory = 10) +
  ggtitle("KEGG ORA - Dot Plot")

barplot(kk, showCategory = 10) +
  ggtitle("KEGG ORA - Bar Plot")


############################################################
# STEP 8: KEGG GSEA
############################################################

# Perform GSEA using KEGG

kk2 <- gseKEGG(
  geneList = geneList,
  organism = "hsa",
  keyType = "ncbi-geneid",
  minGSSize = 10,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# Convert results into a data frame

kk2_df <- as.data.frame(kk2)

# Display results

head(kk2_df)

# Visualize KEGG GSEA results

dotplot(kk2, showCategory = 10) +
  ggtitle("KEGG GSEA - Dot
