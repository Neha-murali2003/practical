# dplyr-https://datacarpentry.github.io/genomics-r-intro/05-dplyr.html
# ============================================================
# EXPERIMENT: DATA MANIPULATION USING DPLYR AND TIDYR
# ============================================================


# 1. INSTALL AND LOAD REQUIRED PACKAGES

# Install packages (run only once)

install.packages("dplyr")
install.packages("tidyr")
install.packages("readr")

# Load packages into the current R session

library(dplyr)
library(tidyr)
library(readr)


# ============================================================
# 2. IMPORTING THE DATASET
# ============================================================

# METHOD 1: Load the dataset from the local computer
# A file selection window will open.
# Select the required CSV file.

data <- read_csv(file.choose())

# Display the imported dataset

View(data)


# METHOD 2: Load the dataset directly from GitHub
# The dataset is imported directly using its URL.
# No manual downloading is required.

data_github <- read_csv(
  "https://raw.githubusercontent.com/datacarpentry/genomics-r-intro/main/episodes/data/combined_tidy_vcf.csv"
)

# Display the imported dataset

View(data_github)


# METHOD 3: Download and save the dataset locally
# The file is downloaded from GitHub and saved as variants.csv.

download.file(
  "https://raw.githubusercontent.com/datacarpentry/genomics-r-intro/main/episodes/data/combined_tidy_vcf.csv",
  destfile = "variants.csv"
)

# Read the downloaded CSV file into R

variants <- read_csv("variants.csv")

# Display the dataset

View(variants)


# NOTE:
# All three methods load the same dataset.
# For the remaining commands, use 'variants' as the
# dataset name if you have executed Method 3.
# If you use Method 1 or Method 2, replace 'variants'
# with 'data' or 'data_github', respectively.


# ============================================================
# 3. EXPLORING THE DATASET
# ============================================================

# Display the first six rows

head(variants)

# Display the last six rows

tail(variants)

# Display the dimensions (number of rows and columns)

dim(variants)

# Display column names

names(variants)

# Display the structure of the dataset

str(variants)

# Display a compact overview of the dataset

glimpse(variants)

# Generate summary statistics

summary(variants)


# ============================================================
# 4. SELECTING COLUMNS USING select()
# ============================================================

# Select specific columns

select(variants, sample_id, REF, ALT, DP)

# Select columns ending with the letter B

select(variants, ends_with("B"))

# Select columns containing the letter i

select(variants, contains("i"))

# Select columns containing i along with POS

select(variants, contains("i"), POS)

# Select columns and exclude Indiv and FILTER

select(variants, contains("i"), POS, -Indiv, -FILTER)


# ============================================================
# 5. FILTERING ROWS USING filter()
# ============================================================

# Select rows belonging to a specific sample

filter(variants, sample_id == "SRR2584863")

# Select rows where depth is greater than 20

filter(variants, DP > 20)

# Select rows where quality is greater than or equal to 100

filter(variants, QUAL >= 100)

# Select rows where depth is not equal to 10

filter(variants, DP != 10)

# AND operator: both conditions must be TRUE

filter(variants, DP > 20 & QUAL > 100)

# OR operator: at least one condition must be TRUE

filter(variants, DP > 20 | QUAL > 100)

# NOT operator: exclude INDEL variants

filter(variants, !INDEL)

# Filter a specific sample using multiple conditions

filter(
  variants,
  sample_id == "SRR2584863",
  MQ >= 50 | QUAL >= 100
)

# Select rows belonging to either of two samples

filter(
  variants,
  sample_id %in% c("SRR2584863", "SRR2584866")
)

# Select rows where DP contains missing values

filter(variants, is.na(DP))

# Select rows where DP is not missing

filter(variants, !is.na(DP))


# ============================================================
# 6. COMBINING filter() AND select() USING PIPE
# ============================================================

# First filter rows, then select required columns

variants %>%
  filter(sample_id == "SRR2584863") %>%
  select(REF, ALT, DP)


# ============================================================
# 7. FILTERING VARIANTS WITHIN A GENOMIC REGION
# ============================================================

# Select variants between positions 1,000,000 and 2,000,000
# Retain variants with quality greater than 200
# Exclude INDEL variants

filter(
  variants,
  POS >= 1e6 &
  POS <= 2e6 &
  QUAL > 200 &
  !INDEL
)


# ============================================================
# 8. SELECTING ROWS USING slice()
# ============================================================

# Create a dataset containing variants from one sample

SRR2584863_variants <- variants %>%
  filter(sample_id == "SRR2584863")

# Select the first six rows

SRR2584863_variants %>%
  slice(1:6)


# ============================================================
# 9. CREATING NEW COLUMNS USING mutate()
# ============================================================

# Calculate polymorphism probability

variants %>%
  mutate(
    POLPROB = 1 - (10 ^ -(QUAL / 10))
  )

# Display selected columns along with the new variable

variants %>%
  mutate(
    POLPROB = 1 - (10 ^ -(QUAL / 10))
  ) %>%
  select(sample_id, POS, QUAL, POLPROB)


# ============================================================
# 10. GROUPING DATA USING group_by()
# ============================================================

# Group variants according to sample ID

variants %>%
  group_by(sample_id)


# ============================================================
# 11. COUNTING OBSERVATIONS USING tally() AND count()
# ============================================================

# Count the number of variants in each sample

variants %>%
  group_by(sample_id) %>%
  tally()

# Count observations using count()

variants %>%
  count(sample_id)

# Count INDEL and non-INDEL variants

variants %>%
  count(INDEL)


# ============================================================
# 12. CALCULATING SUMMARY STATISTICS USING summarize()
# ============================================================

# Calculate mean, median, minimum and maximum depth
# separately for each sample

variants %>%
  group_by(sample_id) %>%
  summarize(
    mean_DP = mean(DP),
    median_DP = median(DP),
    min_DP = min(DP),
    max_DP = max(DP)
  )

# Calculate mean depth while ignoring missing values

variants %>%
  group_by(sample_id) %>%
  summarize(
    mean_DP = mean(DP, na.rm = TRUE)
  )


# ============================================================
# 13. USING group_by() WITH mutate()
# ============================================================

# Calculate mean depth for each sample
# Retain all original rows

variants %>%
  group_by(sample_id) %>%
  mutate(
    mean_DP = mean(DP, na.rm = TRUE)
  )


# ============================================================
# 14. REMOVING GROUPING USING ungroup()
# ============================================================

# Group by sample, calculate mean depth,
# and then remove the grouping structure

variants %>%
  group_by(sample_id) %>%
  mutate(
    mean_DP = mean(DP, na.rm = TRUE)
  ) %>%
  ungroup()


# ============================================================
# 15. CONVERTING LONG FORMAT TO WIDE FORMAT
# USING pivot_wider()
# ============================================================

# Group the data by sample ID and chromosome
# Calculate the mean depth for each group

variants_summary <- variants %>%
  group_by(sample_id, CHROM) %>%
  summarize(
    mean_DP = mean(DP, na.rm = TRUE),
    .groups = "drop"
  )

# Display the summarized dataset

variants_summary

# Convert long format into wide format
# Each sample ID becomes a separate column
# Mean depth values are placed in the corresponding columns

variants_wide <- variants_summary %>%
  pivot_wider(
    names_from = sample_id,
    values_from = mean_DP
  )

# Display the wide-format dataset

variants_wide


# ============================================================
# 16. CONVERTING WIDE FORMAT TO LONG FORMAT
# USING pivot_longer()
# ============================================================

# Convert the wide-format dataset back into long format
# Keep the CHROM column unchanged
# Store former column names in sample_id
# Store the corresponding values in mean_DP

variants_long <- variants_wide %>%
  pivot_longer(
    cols = -CHROM,
    names_to = "sample_id",
    values_to = "mean_DP"
  )

# Display the long-format dataset

variants_long


# ============================================================
# 17. EXAMPLE: USING pivot_longer() DIRECTLY
# ============================================================

# Convert all columns except CHROM into long format

variants_wide %>%
  pivot_longer(
    cols = -CHROM,
    names_to = "sample_id",
    values_to = "mean_DP"
  )


# ============================================================
# 18. EXAMPLE: USING pivot_wider() DIRECTLY
# ============================================================

# Convert the long-format dataset back into wide format

variants_long %>%
  pivot_wider(
    names_from = sample_id,
    values_from = mean_DP
  )


# ============================================================
# 19. COMPLETE DATA MANIPULATION PIPELINE
# ============================================================

# Filter variants from a specific sample
# Retain variants with quality greater than 100
# Exclude INDEL variants
# Calculate polymorphism probability
# Select the required columns

variants %>%
  filter(
    sample_id == "SRR2584863",
    QUAL > 100,
    !INDEL
  ) %>%
  mutate(
    POLPROB = 1 - (10 ^ -(QUAL / 10))
  ) %>%
  select(
    sample_id,
    CHROM,
    POS,
    REF,
    ALT,
    QUAL,
    DP,
    POLPROB
  )


# ============================================================
# END OF EXPERIMENT
# ============================================================
