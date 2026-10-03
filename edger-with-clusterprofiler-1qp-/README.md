# edger-with-clusterprofiler-1qp-

############################################################
# QUESTION 2: DEG ANALYSIS AND PATHWAY ENRICHMENT
# USING edgeR AND clusterProfiler
############################################################


############################################################
# STEP 1: INSTALL AND LOAD PACKAGES
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages <- c(
  "edgeR",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "ggplot2"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}

library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)


############################################################
# STEP 2: IMPORT RAW RNA-SEQ COUNT MATRIX
############################################################

counts <- read.csv(
  file.choose(),
  row.names = 1,
  check.names = FALSE
)

head(counts)
dim(counts)


############################################################
# STEP 3: IMPORT SAMPLE METADATA
############################################################

metadata <- read.csv(
  file.choose(),
  row.names = 1
)

# Match metadata order to count matrix

metadata <- metadata[colnames(counts), , drop = FALSE]

# Define experimental groups

metadata$Condition <- factor(
  metadata$Condition,
  levels = c("Control", "Treatment")
)

table(metadata$Condition)


############################################################
# STEP 4: CREATE DGEList OBJECT
############################################################

counts <- as.matrix(counts)
storage.mode(counts) <- "integer"

y <- DGEList(
  counts = counts,
  group = metadata$Condition
)


############################################################
# STEP 5: FILTER LOWLY EXPRESSED GENES
############################################################

keep <- filterByExpr(y)

y <- y[keep, , keep.lib.sizes = FALSE]

nrow(y)


############################################################
# STEP 6: NORMALIZATION
############################################################

y <- calcNormFactors(y)


############################################################
# STEP 7: ESTIMATE DISPERSION
############################################################

y <- estimateDisp(y)

plotBCV(y)


############################################################
# STEP 8: DIFFERENTIAL EXPRESSION ANALYSIS
############################################################

et <- exactTest(
  y,
  pair = c("Control", "Treatment")
)

results <- topTags(
  et,
  n = Inf
)$table

head(results)


############################################################
# STEP 9: IDENTIFY SIGNIFICANT DEGs
############################################################

# Apply the exact thresholds in the question

deg <- subset(
  results,
  FDR <= 0.05 &
    (logFC <= -2 | logFC >= 2)
)

# Display significant genes

head(deg)


############################################################
# STEP 10: REPORT THE NUMBER OF DEGs
############################################################

total_degs <- nrow(deg)

cat(
  "Total significant DEGs:",
  total_degs,
  "\n"
)

cat(
  "Upregulated genes:",
  sum(deg$logFC >= 2),
  "\n"
)

cat(
  "Downregulated genes:",
  sum(deg$logFC <= -2),
  "\n"
)


############################################################
# STEP 11: SAVE DEG RESULTS
############################################################

write.csv(
  deg,
  "Significant_DEGs.csv",
  row.names = TRUE
)


############################################################
# STEP 12: CONVERT GENE SYMBOLS TO ENTREZ IDs
############################################################

# Use this section if the row names are human gene symbols.

gene_conversion <- bitr(
  rownames(deg),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Remove duplicate Entrez IDs

gene_conversion <- gene_conversion[
  !duplicated(gene_conversion$ENTREZID),
]

# Prepare gene list for pathway analysis

pathway_genes <- gene_conversion$ENTREZID

# Convert the background genes as well

background_conversion <- bitr(
  rownames(y),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

background_genes <- unique(
  background_conversion$ENTREZID
)


############################################################
# STEP 13: PATHWAY ENRICHMENT USING clusterProfiler
############################################################

kk <- enrichKEGG(
  gene = pathway_genes,
  organism = "hsa",
  keyType = "ncbi-geneid",
  universe = background_genes,
  pvalueCutoff = 0.05
)

# Convert enrichment results into a data frame

pathway_results <- as.data.frame(kk)

# Display pathway results

head(pathway_results)


############################################################
# STEP 14: IDENTIFY TOP 5 PATHWAYS
############################################################

# Select top five pathways by adjusted p-value

top5 <- head(
  pathway_results[
    order(pathway_results$p.adjust),
  ],
  5
)

print(top5)


############################################################
# STEP 15: CREATE BAR PLOT
############################################################

if (nrow(top5) > 0) {

  p <- barplot(
    kk,
    showCategory = 5,
    title = "Top 5 KEGG Pathways"
  )

  print(p)

  # Save bar plot

  ggsave(
    "Top5_KEGG_Pathways.png",
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )

} else {

  cat("No significant pathways found.")

}


############################################################
# STEP 16: SAVE PATHWAY RESULTS
############################################################

write.csv(
  pathway_results,
  "KEGG_Pathway_Results.csv",
  row.names = FALSE
)

write.csv(
  top5,
  "Top5_KEGG_Pathways.csv",
  row.names = FALSE
)


############################################################
# STEP 17: FINAL SUMMARY
############################################################

cat("Total significant DEGs:", total_degs, "\n")

cat("Total enriched pathways:", nrow(pathway_results), "\n")

cat("Top five pathways saved successfully.\n")

cat("Analysis completed!")
