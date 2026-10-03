# MICROARRAY

# ============================================================
# EXPERIMENT 4: MICROARRAY DATA ANALYSIS USING LIMMA
# ============================================================


# PART A: INSTALLATION AND LOADING OF PACKAGES

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("limma", quietly = TRUE))
    BiocManager::install("limma")

library(limma)


# Set working directory
getwd()

# Read target file
targets <- readTargets("Targets.txt")

# View target information
targets
head(targets)


# ============================================================
# PART B: IMPORTING RAW MICROARRAY DATA
# ============================================================

# 1. Import using Agilent median signal
RG <- read.maimages(
    targets$FileName,
    source = "agilent"
)

# Inspect raw data
dim(RG)
str(RG)


# 2. Import using Agilent mean signal
y <- read.maimages(
    targets,
    source = "agilent.mean"
)

dim(y)


# 3. Import selected columns for custom analysis
r_only <- read.maimages(
    targets,
    source = "agilent.mean",
    columns = list(
        E = "rMeanSignal",
        Eb = "rBGMedianSignal"
    )
)

dim(r_only)


# ============================================================
# PART C: WITHIN-ARRAY NORMALIZATION AND QUALITY CONTROL
# ============================================================

# Background correction and LOESS normalization
MA <- normalizeWithinArrays(
    RG,
    method = "loess",
    bc.method = "normexp",
    offset = 50
)

# Inspect dimensions and summary statistics
dim(RG)
summary(RG$R)
summary(RG$G)
dim(MA)


# Generate MA plot for sample 2
plotMD(
    RG,
    column = 2,
    main = "MA Plot Before Normalization"
)

# Generate MA plot using sample name
# Run only if sample2 is the exact column name
# plotMD(RG, column = which(colnames(RG) == "sample2"))


# Visualize MA plots across arrays
plotMA3by2(MA)


# ============================================================
# PART D: BETWEEN-ARRAY NORMALIZATION
# ============================================================

# Apply quantile normalization
MA.q <- normalizeBetweenArrays(
    MA,
    method = "quantile"
)


# ============================================================
# DENSITY PLOT COMPARISON
# BEFORE, WITHIN-ARRAY AND BETWEEN-ARRAY NORMALIZATION
# ============================================================

# Save all three plots in one image
png(
    "Density_Comparison.png",
    width = 2400,
    height = 900,
    res = 150
)

# Arrange plots horizontally
par(mfrow = c(1, 3), mar = c(5, 5, 5, 2))

# Before normalization
plotDensities(
    RG,
    main = "Before Normalization"
)

# After within-array normalization
plotDensities(
    MA,
    main = "After Within-array Normalization"
)

# After between-array normalization
plotDensities(
    MA.q,
    main = "After Between-array Normalization"
)

# Close and save image
dev.off()

# Reset plotting layout
par(mfrow = c(1, 1))


# ============================================================
# PART E: LINEAR MODELING AND EMPIRICAL BAYES ANALYSIS
# ============================================================

# Create design matrix
design <- modelMatrix(
    targets,
    ref = "normal"
)

# View design matrix
design


# Fit linear model
fit <- lmFit(
    MA.q,
    design
)


# Apply empirical Bayes moderation
eb <- eBayes(fit)


# Extract top differentially expressed probes
tb <- topTable(
    eb,
    adjust.method = "BH"
)

head(tb)


# Extract complete differential expression table
tt <- topTable(
    eb,
    number = Inf,
    adjust.method = "BH"
)

# View results
head(tt)

# Check dimensions
dim(tt)


# ============================================================
# PART F: FILTERING AND EXPORTING SIGNIFICANT PROBES
# ============================================================

# Filter based on p-value and log-fold change
select <- rownames(tt)[
    tt$P.Value <= 0.05 &
    (tt$logFC >= 0.5 | tt$logFC <= -0.5)
]

# Count selected probes
length(select)


# Extract normalized expression values
selectedProbes <- MA.q[select, ]


# Export selected probes
write.table(
    selectedProbes,
    "outputSelected.txt",
    row.names = TRUE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
)


# Export complete differential expression results
write.table(
    tt,
    "Differential_Expression_Results.txt",
    row.names = TRUE,
    col.names = TRUE,
    sep = "\t",
    quote = FALSE
)


# ============================================================
# PART G: SAVE ANALYSIS OBJECTS
# ============================================================

save(
    targets,
    RG,
    MA,
    MA.q,
    design,
    fit,
    eb,
    tt,
    selectedProbes,
    file = "Microarray_Analysis.RData"
)
# ============================================================
# COMBINED DENSITY PLOT
# Before Normalization, After Within-array Normalization,
# and After Between-array Normalization
# ============================================================

png("Combined_Density_Plot.png",
    width = 2400, height = 900, res = 150)

par(mfrow = c(1, 3),
    mar = c(5, 5, 5, 2))

# 1. Before normalization
plotDensities(RG,
              main = "Before Normalization",
              legend = FALSE)

# 2. After within-array normalization
plotDensities(MA,
              main = "After Within-array Normalization",
              legend = FALSE)

# 3. After between-array normalization
plotDensities(MA.q,
              main = "After Between-array Normalization",
              legend = FALSE)

dev.off()

# Reset plotting layout
par(mfrow = c(1, 1))

# ============================================================
# END OF EXPERIMENT
# ============================================================
