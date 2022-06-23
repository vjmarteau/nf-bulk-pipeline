#!/usr/bin/env Rscript
'
Usage:
  Reformat_pseudobulk.R --gene_counts=<gene_counts> --samplesheet=<samplesheet> --prefix=<prefix> [options]

Mandatory arguments:
  --gene_counts=<gene_counts>   nfcore/rnaseq gene counts in TSV format
  --samplesheet=<samplesheet>   Samplesheet from nfcore/rnaseq pipeline
  --prefix=<prefix>             Prefix for output filenames

Optional arguments:
  --results_dir=<dir>         Output directory [default: ./]
' -> doc

# load required packages
library(docopt)
arguments <- docopt(doc, version = "0.1")
print(arguments)

library(conflicted)
library(tidyverse)

conflict_prefer("rename", "dplyr")
conflict_prefer("select", "dplyr")

# Load parameters
data <- read_tsv(arguments$gene_counts)
metadata <- read_csv(arguments$samplesheet)
prefix <- arguments$prefix
results_dir <- arguments$results_dir

# Data wrangling
metadata <- metadata |>
  mutate(label = internal_id) |>
  mutate(condition = factor(group)) |>
  mutate(condition = fct_relevel(condition, c("Ctrl", "A8", "A9", "A8A9"))) |>
  mutate(patient = factor(gsub("-.*", "", internal_id))) |>
  mutate(order = parse_number(sample)) |>
  arrange(order) |>
  select(sample, label, condition, patient) |>
  column_to_rownames(var = "sample")

count_mat <- data |>
  mutate(gene_id = gsub("\\.[0-9]*$", "", gene_id)) |>
  column_to_rownames(var = "gene_id") |>
  select(-gene_name) |>
  mutate_all(as.integer) |>
  select(rownames(metadata))

key <- data |>
  select(gene_id, gene_name) |>
  mutate(gene_id = gsub("\\.[0-9]*$", "", gene_id))

lapply(c("metadata", "count_mat", "key"), function(x) saveRDS(get(x), file.path(results_dir, paste0(prefix, "_", x, ".rds"))))