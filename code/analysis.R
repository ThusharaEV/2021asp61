# Dengue Virus Bioinformatics Project
# ICA 3 - Bioinformatics

# Project information


project_name <- "Dengue Virus Sequence Analysis"

research_topic <- "Dengue Virus Serotype Analysis"

print(project_name)

print(research_topic)

# Display project initialization message

print("Project setup completed successfully!")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("Biostrings")


# Load Biostrings package
library(Biostrings)
library(dplyr)
library(ggplot2)

# Import four dengue serotype FASTA files
denv1 <- readDNAStringSet("data/DENV1.fasta")
denv2 <- readDNAStringSet("data/DENV2.fasta")
denv3 <- readDNAStringSet("data/DENV3.fasta")
denv4 <- readDNAStringSet("data/DENV4.fasta")

# Add serotype names
names(denv1) <- paste0("DENV1_", seq_along(denv1))
names(denv2) <- paste0("DENV2_", seq_along(denv2))
names(denv3) <- paste0("DENV3_", seq_along(denv3))
names(denv4) <- paste0("DENV4_", seq_along(denv4))

# Combine all sequences
all_sequences <- c(denv1, denv2, denv3, denv4)

# Check total number of sequences
length(all_sequences)

# Display sequence lengths
width(all_sequences)

# Display sequence names
names(all_sequences)

# Check number of sequences per serotype
cat("DENV-1:", length(denv1), "\n")
cat("DENV-2:", length(denv2), "\n")
cat("DENV-3:", length(denv3), "\n")
cat("DENV-4:", length(denv4), "\n")
cat("Total:", length(all_sequences), "\n")

# Calculate sequence lengths
seq_lengths <- width(all_sequences)

# Calculate GC content
gc_values <- letterFrequency(
  all_sequences,
  letters = c("G", "C")
)

gc_content <- rowSums(gc_values) / seq_lengths * 100

# Create analysis data frame
results <- data.frame(
  Sequence = names(all_sequences),
  Length = seq_lengths,
  GC_Content = gc_content
)

# Display first results
head(results)

# Summary statistics
summary(results)

# Identify the serotype of each sequence
serotype <- c(
  rep("DENV-1", length(denv1)),
  rep("DENV-2", length(denv2)),
  rep("DENV-3", length(denv3)),
  rep("DENV-4", length(denv4))
)

# Count A, T, G and C nucleotides
nuc_counts <- letterFrequency(
  all_sequences,
  letters = c("A", "T", "G", "C")
)

# Calculate percentages using A, T, G and C bases
valid_bases <- rowSums(nuc_counts)

nuc_percent <- sweep(
  nuc_counts, 1, valid_bases, "/"
) * 100

# Create the final analysis table
results <- data.frame(
  Sequence = names(all_sequences),
  Serotype = serotype,
  Length = width(all_sequences),
  A = nuc_percent[, "A"],
  T = nuc_percent[, "T"],
  G = nuc_percent[, "G"],
  C = nuc_percent[, "C"],
  GC_Content = nuc_percent[, "G"] + nuc_percent[, "C"]
)

# Display results
head(results)

# Summary by serotype
composition_summary <- results %>%
  group_by(Serotype) %>%
  summarise(
    Mean_Length = mean(Length),
    Mean_GC = mean(GC_Content),
    Mean_A = mean(A),
    Mean_T = mean(T),
    Mean_G = mean(G),
    Mean_C = mean(C),
    .groups = "drop"
  )

print(composition_summary)


p1 <- ggplot(results, aes(x = Serotype,
                          y = GC_Content,
                          fill = Serotype)) +
  geom_boxplot() +
  labs(
    title = "GC Content Across Dengue Virus Serotypes",
    x = "Dengue Serotype",
    y = "GC Content (%)"
  ) +
  theme_minimal()

print(p1)

ggsave("results/figures/gc_content.png",
       plot = p1, width = 8, height = 5)


p2 <- ggplot(results, aes(x = Serotype,
                          y = Length,
                          fill = Serotype)) +
  geom_boxplot() +
  labs(
    title = "Sequence Length Across Dengue Virus Serotypes",
    x = "Dengue Serotype",
    y = "Sequence Length (bases)"
  ) +
  theme_minimal()

print(p2)

ggsave("results/figures/sequence_length.png",
       plot = p2, width = 8, height = 5)

dir.create("results/figures",
           recursive = TRUE,
           showWarnings = FALSE)

table(results$Serotype)

summary(results$Length)

# Save analysis results

write.csv(
  results,
  "results/sequence_results.csv",
  row.names = FALSE
)

write.csv(
  composition_summary,
  "results/composition_summary.csv",
  row.names = FALSE
)

# Check number of sequences per serotype
table(results$Serotype)

# Check sequence length distribution
summary(results$Length)

# Check ambiguous N bases
n_counts <- letterFrequency(
  all_sequences,
  letters = "N"
)

summary(n_counts)

n_quality <- data.frame(
  Sequence = names(all_sequences),
  N_Count = as.numeric(n_counts)
)

n_quality <- n_quality %>%
  arrange(desc(N_Count))

print(head(n_quality, 10))

# Calculate N percentage for each sequence

n_quality$Length <- as.numeric(
  width(all_sequences)[match(
    n_quality$Sequence,
    names(all_sequences)
  )]
)

n_quality$N_Percentage <- (
  n_quality$N_Count / n_quality$Length
) * 100

print(head(n_quality, 10))

# Identify sequences with more than 5% ambiguous bases

n_quality %>%
  filter(N_Percentage > 5)

# Remove sequences with more than 5% ambiguous bases

excluded_sequences <- n_quality %>%
  filter(N_Percentage > 5)

results_clean <- results %>%
  filter(!Sequence %in% excluded_sequences$Sequence)

# Check cleaned dataset
print(table(results_clean$Serotype))

# Check remaining sequence lengths
summary(results_clean$Length)

# Save cleaned dataset
write.csv(
  results_clean,
  "results/cleaned_sequence_results.csv",
  row.names = FALSE
)

# Sequence length summary by serotype

results_clean %>%
  group_by(Serotype) %>%
  summarise(
    Number_of_Sequences = n(),
    Minimum_Length = min(Length),
    Maximum_Length = max(Length),
    Mean_Length = mean(Length),
    .groups = "drop"
  )


# Select near-complete sequences

results_near_complete <- results_clean %>%
  filter(Length >= 9000)

# Count sequences by serotype

results_near_complete %>%
  group_by(Serotype) %>%
  summarise(
    Number_of_Sequences = n(),
    Mean_Length = mean(Length),
    Mean_GC = mean(GC_Content),
    .groups = "drop"
  )

# Save near-complete dataset

write.csv(
  results_near_complete,
  "results/near_complete_sequences.csv",
  row.names = FALSE
)


p3 <- ggplot(results_near_complete,
             aes(x = Serotype, y = GC_Content, fill = Serotype)) +
  geom_boxplot() +
  labs(
    title = "GC Content of Near-Complete Dengue Virus Sequences",
    x = "Dengue Serotype",
    y = "GC Content (%)"
  ) +
  theme_minimal()

print(p3)

ggsave(
  "results/figures/near_complete_gc_content.png",
  plot = p3, width = 8, height = 5
)


# P4: Sequence Length Comparison of Near-Complete Sequences

p4 <- ggplot(
  results_near_complete,
  aes(x = Serotype, y = Length, fill = Serotype)
) +
  geom_boxplot() +
  labs(
    title = "Sequence Length Comparison of Near-Complete Dengue Genomes",
    x = "Dengue Serotype",
    y = "Sequence Length (bases)"
  ) +
  theme_minimal()

print(p4)

ggsave(
  "results/figures/near_complete_sequence_length.png",
  plot = p4,
  width = 8,
  height = 5
)


# Save summary statistics

near_complete_summary <- results_near_complete %>%
  group_by(Serotype) %>%
  summarise(
    Number_of_Sequences = n(),
    Mean_Length = mean(Length),
    Mean_GC = mean(GC_Content),
    Mean_A = mean(A),
    Mean_T = mean(T),
    Mean_G = mean(G),
    Mean_C = mean(C),
    .groups = "drop"
  )

write.csv(
  near_complete_summary,
  "results/near_complete_summary.csv",
  row.names = FALSE
)

print(near_complete_summary)
