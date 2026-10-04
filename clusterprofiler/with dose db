# EXPERIMENT: GO CLASSIFICATION, ENRICHMENT ANALYSIS AND GSEA USING CLUSTERPROFILER

# --------------------------------------------------
# STEP 1: INSTALLATION AND LOADING OF PACKAGES
# --------------------------------------------------

# Install BiocManager if it is not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install the required Bioconductor packages for GO and pathway analysis
BiocManager::install(c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "DOSE",
  "enrichplot",
  "ReactomePA",
  "pathview"
), ask = FALSE, update = FALSE)

# Install ggplot2 for customized data visualization
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

# Load all the required packages into the R environment
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ReactomePA)
library(pathview)
library(ggplot2)


# --------------------------------------------------
# STEP 2: LOAD AND EXPLORE THE GENE LIST
# --------------------------------------------------

# Load the predefined ranked geneList dataset from the DOSE package
data(geneList, package = "DOSE")

# Display the first six genes in the dataset
head(geneList)

# Display the last six genes in the dataset
tail(geneList)

# Examine the structure and data type of the gene list
str(geneList)

# Calculate the total number of genes present in the list
length(geneList)

# Display the gene identifiers associated with the values
head(names(geneList))


# --------------------------------------------------
# STEP 3: GENE ONTOLOGY (GO) CLASSIFICATION
# --------------------------------------------------

# Perform GO classification using Biological Process ontology at level 3
ggo <- groupGO(
  gene = gene,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  level = 3,
  readable = TRUE
)

# Display the first few GO classification results
head(as.data.frame(ggo))

# Calculate the total number of GO categories identified
nrow(as.data.frame(ggo))


# --------------------------------------------------
# STEP 3.1: GO CLASSIFICATION BAR PLOT
# --------------------------------------------------

# Generate a bar plot displaying the top 10 GO categories
p_ggo <- barplot(
  ggo,
  drop = TRUE,
  showCategory = 10,
  title = "GO Classification - Biological Process"
)

# Wrap long GO term labels and adjust the font size and title alignment
p_ggo <- p_ggo +
  scale_y_discrete(
    labels = function(x) {
      vapply(x, function(y) {
        paste(strwrap(y, width = 30), collapse = "\n")
      }, character(1))
    }
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

# Display the GO classification bar plot
print(p_ggo)

# Save the bar plot as a high-resolution PNG image
ggsave(
  "GO_Classification_Barplot.png",
  plot = p_ggo,
  width = 12,
  height = 8,
  dpi = 300
)


# --------------------------------------------------
# STEP 3.2: GO CLASSIFICATION DOT PLOT
# --------------------------------------------------

# Convert GO classification results into a data frame
ggo_df <- as.data.frame(ggo)

# Display the available columns in the data frame
colnames(ggo_df)

# Select the top 10 GO categories based on decreasing gene count
ggo_top <- head(
  ggo_df[order(ggo_df$Count, decreasing = TRUE), ],
  10
)

# Create a dot plot showing gene counts across GO categories
p_ggo_dot <- ggplot(
  ggo_top,
  aes(
    x = Count,
    y = reorder(Description, Count),
    size = Count,
    color = Count
  )
) +
  geom_point(alpha = 0.8) +
  scale_y_discrete(
    labels = function(x) {
      vapply(x, function(y) {
        paste(strwrap(y, width = 30), collapse = "\n")
      }, character(1))
    }
  ) +
  scale_color_gradient(
    low = "blue",
    high = "red"
  ) +
  theme_minimal() +
  labs(
    title = "GO Classification - Biological Process",
    x = "Gene Count",
    y = "GO Terms",
    size = "Count",
    color = "Count"
  ) +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )

# Display the GO classification dot plot
print(p_ggo_dot)

# Save the dot plot as a high-resolution PNG image
ggsave(
  "GO_Classification_Dotplot.png",
  plot = p_ggo_dot,
  width = 12,
  height = 8,
  dpi = 300
)


# --------------------------------------------------
# STEP 4: GO OVER-REPRESENTATION ANALYSIS (ORA)
# --------------------------------------------------

# Perform GO enrichment analysis using the selected genes
ego <- enrichGO(
  gene = gene,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


# --------------------------------------------------
# STEP 4.1: GO ORA BAR PLOT
# --------------------------------------------------

# Generate a bar plot showing the top 10 enriched GO terms
p_ora_bar <- barplot(
  ego,
  showCategory = 10,
  title = "GO ORA - Top Enriched Terms"
)

# Display the ORA bar plot
print(p_ora_bar)

# Generate a bar plot without specifying the number of categories
p_ora_bar_all <- barplot(
  ego,
  title = "GO ORA - All Enriched Terms"
)

# Display the bar plot containing all available enriched categories
print(p_ora_bar_all)


# --------------------------------------------------
# STEP 4.2: GO ORA DOT PLOT
# --------------------------------------------------

# Create a dot plot showing the top 10 enriched GO terms
p_ora_dot <- dotplot(
  ego,
  showCategory = 10,
  title = "GO ORA - Top Enriched Terms"
)

# Display the ORA dot plot
print(p_ora_dot)

# Save the dot plot as a high-resolution PNG image
ggsave(
  "GO_ORA_Dotplot.png",
  plot = p_ora_dot,
  width = 12,
  height = 8,
  dpi = 300
)

# Export the complete GO ORA results into a CSV file
write.csv(
  as.data.frame(ego),
  "GO_ORA_Results.csv",
  row.names = FALSE
)


# --------------------------------------------------
# STEP 5: GENE SET ENRICHMENT ANALYSIS (GSEA)
# --------------------------------------------------

# Display the first few genes in the ranked gene list
head(geneList)

# Check the total number of genes in the ranked list
length(geneList)

# Check whether the gene list is arranged in decreasing order
head(sort(geneList, decreasing = TRUE))

# Sort the gene list in decreasing order if it is not already ordered
# geneList <- sort(geneList, decreasing = TRUE)


# --------------------------------------------------
# STEP 5.1: PERFORM GSEA
# --------------------------------------------------

# Perform Gene Set Enrichment Analysis using the complete ranked gene list
ego3 <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)


# --------------------------------------------------
# STEP 5.2: EXPLORE AND EXPORT GSEA RESULTS
# --------------------------------------------------

# Convert GSEA results into a data frame for inspection and further analysis
gsea_df <- as.data.frame(ego3)

# Display the first few GSEA results
head(gsea_df)

# Calculate the total number of enriched GO terms
nrow(gsea_df)

# Display the column names available in the GSEA results
colnames(gsea_df)

# Save the GSEA results as a CSV file
write.csv(
  gsea_df,
  "GO_GSEA_Results.csv",
  row.names = FALSE
)


# --------------------------------------------------
# STEP 5.3: GSEA DOT PLOT
# --------------------------------------------------

# Create a dot plot to visualize the enriched GO terms identified by GSEA
p_dot <- dotplot(
  ego3,
  title = "GO Gene Set Enrichment Analysis - Dot Plot"
)

# Display the GSEA dot plot
print(p_dot)

# Save the dot plot as a high-resolution PNG image
ggsave(
  "GSEA_Dotplot.png",
  plot = p_dot,
  width = 14,
  height = 9,
  dpi = 300
)


# --------------------------------------------------
# STEP 5.4: GSEA DOT PLOT WITH IMPROVED READABILITY
# --------------------------------------------------

# Display the top 10 GO terms and wrap long labels to reduce overcrowding
p_dot <- dotplot(
  ego3,
  showCategory = 10,
  label_format = 35,
  font.size = 12,
  title = "GO Gene Set Enrichment Analysis - Dot Plot"
)

# Display the improved GSEA dot plot
print(p_dot)

# Save the plot with increased dimensions for better readability
ggsave(
  "GSEA_Dotplot.png",
  plot = p_dot,
  width = 14,
  height = 10,
  dpi = 300
)


# --------------------------------------------------
# STEP 5.5: GSEA BAR PLOT
# --------------------------------------------------

# Convert the GSEA results into a data frame for customized plotting
gsea_df <- as.data.frame(ego3)

# Select the top 15 GO terms based on the lowest adjusted p-values
top15 <- head(
  gsea_df[order(gsea_df$p.adjust), ],
  15
)

# Create a horizontal bar plot using NES values and adjusted p-values
p_bar <- ggplot(
  top15,
  aes(
    x = reorder(Description, NES),
    y = NES,
    fill = p.adjust
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "GO Gene Set Enrichment Analysis - Bar Plot",
    x = "GO Terms",
    y = "Normalized Enrichment Score (NES)",
    fill = "Adjusted P-value"
  ) +
  theme_minimal()

# Display the GSEA bar plot
print(p_bar)

# Save the bar plot as a high-resolution PNG image
ggsave(
  "GSEA_Barplot.png",
  plot = p_bar,
  width = 12,
  height = 10,
  dpi = 300
)


# --------------------------------------------------
# STEP 5.6: GSEA RIDGE PLOT
# --------------------------------------------------

# Install ggridges if it is not already installed
install.packages("ggridges")

# Load the package required for creating ridge plots
library(ggridges)

# Create a ridge plot displaying the top 10 enriched GO terms
p_ridge <- ridgeplot(
  ego3,
  showCategory = 10
) +
  ggtitle("GO Gene Set Enrichment Analysis - Ridge Plot")

# Display the GSEA ridge plot
print(p_ridge)

# Save the ridge plot as a high-resolution PNG image
ggsave(
  "GSEA_Ridgeplot.png",
  plot = p_ridge,
  width = 12,
  height = 9,
  dpi = 300
)
