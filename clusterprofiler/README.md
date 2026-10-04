# clusterprofiler

# -----------------------------------------------
# STEP 1: Install required packages
# -----------------------------------------------

# Install BiocManager if it is not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor packages if required
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

if (!requireNamespace("DOSE", quietly = TRUE)) {
  BiocManager::install("DOSE")
}

if (!requireNamespace("enrichplot", quietly = TRUE)) {
  BiocManager::install("enrichplot")
}

# Install ggplot2 if required
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}


# -----------------------------------------------
# STEP 2: Load the required libraries
# -----------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggplot2)


# -----------------------------------------------
# STEP 3: Load the example gene list
# -----------------------------------------------

data(geneList, package = "DOSE")


# -----------------------------------------------
# STEP 4: Perform GO Gene Set Enrichment Analysis
# -----------------------------------------------

ego3 <- gseGO(
  geneList = geneList,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  minGSSize = 100,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE
)


# -----------------------------------------------
# STEP 5: View the GSEA results
# -----------------------------------------------

# Display the first few enriched GO terms
head(as.data.frame(ego3))

# Count the number of GO terms returned
nrow(as.data.frame(ego3))


# -----------------------------------------------
# STEP 6: Create a dot plot
# -----------------------------------------------

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
  width = 10,
  height = 8,
  dpi = 300
)


# -----------------------------------------------
# STEP 7: Create a less crowded bar plot
# Display only the top 15 GO terms
# -----------------------------------------------

# Select the 15 terms with the lowest adjusted p-values
top15 <- head(
  as.data.frame(ego3)[order(as.data.frame(ego3)$p.adjust), ],
  15
)

# Create the bar plot
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

# Display the bar plot
print(p_bar)

# Save the bar plot
ggsave(
  "GSEA_Barplot.png",
  plot = p_bar,
  width = 12,
  height = 10,
  dpi = 300
)


# Create a bar plot using the top 15 GO terms
ggplot(
  head(
    as.data.frame(ego3)[order(as.data.frame(ego3)$p.adjust), ],
    15
  ),
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

  #or
  ggsave(
 # "GSEA_Barplot.png",
 # width = 14,
 # height = 18,
  #dpi = 300
#)
# -----------------------------------------------
# STEP 8: Check the location of saved plots
# -----------------------------------------------

getwd()
