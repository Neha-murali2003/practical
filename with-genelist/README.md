```r
# ============================================================
# EXPERIMENT: GO CLASSIFICATION, ENRICHMENT AND PATHWAY ANALYSIS
# ============================================================


# STEP 1: INSTALLATION OF REQUIRED PACKAGES

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required Bioconductor packages
BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "ReactomePA"
), ask = FALSE, update = FALSE)

# Install required CRAN packages
install.packages(c("ggplot2", "ggridges"))

# Load all required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ReactomePA)
library(ggplot2)
library(ggridges)


# STEP 2: IMPORT THE CSV FILE

# Select the gene list CSV file
data <- read.csv(file.choose())

# View the first few rows
head(data)

# Check column names
colnames(data)

# Check dataset dimensions
dim(data)


# STEP 3: PREPROCESSING AND ID CONVERSION

# Remove rows containing missing values
data <- na.omit(data)

# Convert gene symbols into Entrez IDs
gene_ids <- bitr(
  data$gene_symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Combine gene symbols, logFC values and Entrez IDs
gene_data <- merge(
  data,
  gene_ids,
  by.x = "gene_symbol",
  by.y = "SYMBOL"
)

# Remove duplicate Entrez IDs
gene_data <- gene_data[!duplicated(gene_data$ENTREZID), ]

# View converted gene IDs
head(gene_data)

# Check dimensions of converted dataset
dim(gene_data)


# STEP 4: GO CLASSIFICATION

# Perform GO classification using Biological Process
go_classification <- groupGO(
  gene = gene_data$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  level = 3,
  readable = TRUE
)

# View GO classification results
head(as.data.frame(go_classification))

# Convert results into a data frame
go_classification_df <- as.data.frame(go_classification)

# Check dimensions of the results
dim(go_classification_df)

#### BarPlot GO classification
barplot(go_classification)

# Save the plot
ggsave("GO_Classification.png", width = 10, height = 8)


#####3plots in one frame (if only reqiured)
# Convert all three GO classification results into data frames

bp_df <- as.data.frame(go_bp)
cc_df <- as.data.frame(go_cc)
mf_df <- as.data.frame(go_mf)

# Add ontology labels
bp_df$ONTOLOGY <- "BP"
cc_df$ONTOLOGY <- "CC"
mf_df$ONTOLOGY <- "MF"

# Combine all three data frames
go_combined <- rbind(bp_df, cc_df, mf_df)

# Select the top 8 classifications from each ontology
go_top <- do.call(
    rbind,
    lapply(
        split(go_combined, go_combined$ONTOLOGY),
        function(x) head(x[order(-x$Count), ], 8)
    )
)

# Plot the combined GO classification results
ggplot(
    go_top,
    aes(x = reorder(Description, Count), y = Count, fill = ONTOLOGY)
) +
    geom_col() +
    coord_flip() +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    labs(
        title = "Top GO Classifications Across Ontologies",
        x = NULL,
        y = "Gene Count",
        fill = "ONTOLOGY"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 8)
    )



#####CC ALONE AS IN GRAPH PLOT
# Load ggplot2
library(ggplot2)
# BP,CC,MF IN A PLOT FOR GO
# Convert GO classification results into a data frame
cc_df <- as.data.frame(go_cc)

# Plot Cellular Component classification
ggplot(cc_df, aes(x = reorder(Description, Count), y = Count)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(
        title = "GO Classification (Cellular Component)",
        x = NULL,
        y = "Count"
    ) +
    theme_minimal() +
    theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
    )





# STEP 5: GO OVER-REPRESENTATION ANALYSIS (ORA)

# Select genes with absolute logFC greater than or equal to 1
selected_genes <- gene_data$ENTREZID[
  abs(gene_data$logFC) >= 1
]

# View selected genes
selected_genes

# Perform GO enrichment analysis
go_enrichment <- enrichGO(
  gene = selected_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# View GO enrichment results
head(as.data.frame(go_enrichment))

# Convert results into a data frame
go_enrichment_df <- as.data.frame(go_enrichment)

# Check dimensions of the results
dim(go_enrichment_df)


# BAR PLOT OF GO ENRICHMENT

# Plot GO enrichment results
barplot(go_enrichment)

# Save the plot
ggsave("GO_ORA_Barplot.png", width = 10, height = 8)


# DOT Plot GO enrichment results
dotplot(go_enrichment) +
  theme(axis.text.y = element_text(size = 8))

# Save the plot
ggsave("GO_ORA_Dotplot.png",
       width = 12,
       height = 12,
       dpi = 300)

#GSEA
# Create a ranked gene list using logFC
gene_list <- gene_data$logFC

# Assign Entrez IDs as names
names(gene_list) <- gene_data$ENTREZID

# Sort genes in decreasing order
gene_list <- sort(gene_list, decreasing = TRUE)

# View the ranked gene list
head(gene_list)
# Perform GSEA for GO Biological Process
go_gsea <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose = FALSE
)

# View GSEA results
head(as.data.frame(go_gsea))

# Convert results into a data frame
go_gsea_df <- as.data.frame(go_gsea)

# Check dimensions
dim(go_gsea_df)
# Plot GSEA results using dot plot
dotplot(go_gsea)

# Save the plot
ggsave("GSEA_Dotplot.png",
       width = 12,
       height = 10,
       dpi = 300)


# Convert GSEA results into a data frame
gsea_df <- as.data.frame(go_gsea)

# Select the top 10 pathways based on adjusted p-value
top10_gsea <- head(gsea_df[order(gsea_df$p.adjust), ], 10)

# Create a simple bar plot
ggplot(top10_gsea, aes(x = reorder(Description, NES), y = NES)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 GSEA Enriched Pathways",
    x = "GO Terms",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal()

# Save the plot
ggsave("GSEA_Barplot.png", width = 12, height = 8, dpi = 300)


# Plot GSEA results using ridge plot
ridgeplot(go_gsea, showCategory = 10)

# Save the plot
ggsave("GSEA_Ridgeplot.png",
       width = 12,
       height = 8,
       dpi = 300)



# STEP 7: KEGG PATHWAY ENRICHMENT

# Perform KEGG enrichment analysis
kegg_results <- enrichKEGG(
  gene = selected_genes,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# View KEGG enrichment results
head(as.data.frame(kegg_results))

# Convert results into a data frame
kegg_df <- as.data.frame(kegg_results)

# Check dimensions
dim(kegg_df)

# Plot KEGG enrichment results
barplot(kegg_results)

# Save the plot
ggsave("KEGG_Barplot.png", width = 12, height = 8, dpi = 300)
# Plot KEGG enrichment results
dotplot(kegg_results)

# Save the plot
ggsave("KEGG_Dotplot.png", width = 12, height = 8, dpi = 300)

# Export KEGG enrichment results
write.csv(kegg_df,
          "KEGG_Results.csv",
          row.names = FALSE)


# STEP 8: REACTOME PATHWAY ENRICHMENT

# Perform Reactome pathway enrichment
reactome_results <- enrichPathway(
  gene = selected_genes,
  organism = "human",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  readable = TRUE
)

# View Reactome enrichment results
head(as.data.frame(reactome_results))

# Convert results into a data frame
reactome_df <- as.data.frame(reactome_results)

# Check dimensions
dim(reactome_df)

# Plot Reactome enrichment results
dotplot(reactome_results) +
  theme(axis.text.y = element_text(size = 9))

# Save the plot
ggsave("Reactome_Dotplot.png",
       width = 14,
       height = 10,
       dpi = 300)

# Plot Reactome enrichment results
barplot(reactome_results)

# Save the plot
ggsave("Reactome_Barplot.png",
       width = 12,
       height = 8,
       dpi = 300)
# Export Reactome enrichment results
write.csv(reactome_df,
          "Reactome_Results.csv",
          row.names = FALSE)

# ============================================================
# END OF CURRENT WORKFLOW
# ============================================================
```
