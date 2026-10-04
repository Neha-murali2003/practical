# Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c(
    "clusterProfiler",
    "org.Hs.eg.db",
    "DOSE",
    "enrichplot",
    "ReactomePA",
    "pathview"
), ask = FALSE, update = FALSE)

# Install ggplot2
if (!requireNamespace("ggplot2", quietly = TRUE))
    install.packages("ggplot2")

# Load all packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ReactomePA)
library(pathview)
library(ggplot2)

# Load geneList dataset from DOSE
data(geneList, package = "DOSE")

# Display the first few genes
head(geneList)

# Display the last few genes
tail(geneList)

# Check the structure
str(geneList)

# Find the total number of genes
length(geneList)

# View gene identifiers
head(names(geneList))

# 3. GO Classification

# Perform GO classification
ggo <- groupGO(
    gene = gene,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    level = 3,
    readable = TRUE
)

# Display results
head(as.data.frame(ggo))

# Check the number of GO categories
nrow(as.data.frame(ggo))

# 3.1 GO Classification Bar Plot

# Generate bar plot with fewer categories
p_ggo <- barplot(
    ggo,
    drop = TRUE,
    showCategory = 10,
    title = "GO Classification - Biological Process"
)

# Wrap long GO term labels
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

# Display plot
print(p_ggo)

# Save high-resolution plot
ggsave(
    "GO_Classification_Barplot.png",
    plot = p_ggo,
    width = 12,
    height = 8,
    dpi = 300
)


# 3.3 GO Classification Dot Plot

# Convert GO classification results into a data frame
ggo_df <- as.data.frame(ggo)

# View available columns
colnames(ggo_df)

# Select top 10 GO categories
ggo_top <- head(
    ggo_df[order(ggo_df$Count, decreasing = TRUE), ],
    10
)

# Create dot plot
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

# Display plot
print(p_ggo_dot)

# Save high-resolution plot
ggsave(
    "GO_Classification_Dotplot.png",
    plot = p_ggo_dot,
    width = 12,
    height = 8,
    dpi = 300
)



# 4.1 GO Over-Representation Analysis (ORA)

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

# 4.2 GO ORA Bar Plot

p_ora_bar <- barplot(
    ego,
    showCategory = 10,
    title = "GO ORA - Top Enriched Terms"
)

print(p_ora_bar)

#  OR Bar plot without specifying the number of categories
p_ora_bar_all <- barplot(
    ego,
    title = "GO ORA - All Enriched Terms"
)

print(p_ora_bar_all)

# 4.3 GO ORA Dot Plot

p_ora_dot <- dotplot(
    ego,
    showCategory = 10,
    title = "GO ORA - Top Enriched Terms"
)

# Display plot
print(p_ora_dot)

# Save high-resolution plot
ggsave(
    "GO_ORA_Dotplot.png",
    plot = p_ora_dot,
    width = 12,
    height = 8,
    dpi = 300
)
# Export GO ORA results
write.csv(
    as.data.frame(ego),
    "GO_ORA_Results.csv",
    row.names = FALSE
)

#GSEA
# 6.1 Check the ranked gene list
head(geneList)

# Check the total number of genes
length(geneList)

# Check whether the gene list is sorted
head(sort(geneList, decreasing = TRUE))

#IF NOT ODERED
## Sort gene list in decreasing order
#geneList <- sort(geneList, decreasing = TRUE)

# 6.2 Perform Gene Set Enrichment Analysis
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
# Convert GSEA results into a data frame
gsea_df <- as.data.frame(ego3)

# Display the first few results
head(gsea_df)

# Check the number of enriched terms
nrow(gsea_df)

# View column names
colnames(gsea_df)


# Save GSEA results as a CSV file
write.csv(
    gsea_df,
    "GO_GSEA_Results.csv",
    row.names = FALSE
)


#OR # STEP 6: Create a GSEA dot plot

p_dot <- dotplot(
  ego3,
  title = "GO Gene Set Enrichment Analysis - Dot Plot"
)

# Display the dot plot
print(p_dot)

# Save the dot plot
ggsave(
  "GSEA_Dotplot.png",
  plot = p_dot,
  width = 14,
  height = 9,
  dpi = 300
)
# -----------------------------------------------
# STEP 6: Create a clear GSEA dot plot
# -----------------------------------------------

# Display only the top 10 GO terms
p_dot <- dotplot(
  ego3,
  showCategory = 10,
  label_format = 35,
  font.size = 12,
  title = "GO Gene Set Enrichment Analysis - Dot Plot"
)

# Display the dot plot
print(p_dot)

# Save the plot with larger dimensions
ggsave(
  "GSEA_Dotplot.png",
  plot = p_dot,
  width = 14,
  height = 10,
  dpi = 300
)

# STEP 7: Create a GSEA bar plot


# Convert GSEA results into a data frame
gsea_df <- as.data.frame(ego3)

# Select top 15 terms
top15 <- head(
  gsea_df[order(gsea_df$p.adjust), ],
  15
)

# Create bar plot
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

# Display the plot
print(p_bar)

# Save the plot
ggsave(
  "GSEA_Barplot.png",
  plot = p_bar,
  width = 12,
  height = 10,
  dpi = 300
)


#STEP 8: Create a GSEA Ridge Plot
install.packages("ggridges")
# -----------------------------------------------
# STEP 8: Create a GSEA Ridge Plot
# -----------------------------------------------

# Load the required package
library(ggridges)

# Create the ridge plot
p_ridge <- ridgeplot(ego3, showCategory = 10) +
  ggtitle("GO Gene Set Enrichment Analysis - Ridge Plot")

# Display the plot
print(p_ridge)

# Save the plot
ggsave(
  "GSEA_Ridgeplot.png",
  plot = p_ridge,
  width = 12,
  height = 9,
  dpi = 300
)
