#!/usr/bin/env Rscript
'
Usage:
  Reformat_data.R --gene_counts=<gene_counts> --samplesheet=<samplesheet> --prefix=<prefix> [options]

Mandatory arguments:
  --gene_counts=<gene_counts>   nfcore/rnaseq gene counts in TSV format
  --samplesheet=<samplesheet>   Samplesheet from nfcore/rnaseq pipeline
  --prefix=<prefix>             Prefix for output filenames

Optional arguments:
  --resDir=<resDir>         Output directory [default: ./]
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
metadata <- read_csv(arguments$samplesheet)
data <- read_tsv(arguments$gene_counts)
prefix <- arguments$prefix
resDir <- arguments$resDir

# Data wrangling
metadata <- metadata |>
  mutate(label = internal_id, treatment = group) |>
  mutate(patient = gsub("-.*", "", internal_id)) |>
  mutate(order = parse_number(sample)) |>
  arrange(order) |>
  select(sample, patient, treatment, label)

count_mat <- data |>
  mutate(gene_id = gsub("\\.[0-9]*$", "", gene_id)) |>
  select(-gene_name) |>
  mutate_at(vars(-("gene_id")), ceiling) |>
  select(c("gene_id", metadata |> pull(1) )) |>
  rename_with(.col = "gene_id", ~"ensembl")

gene_cnvan_key <- data |>
  select(gene_id, gene_name) |>
  mutate(gene_id = gsub("\\.[0-9]*$", "", gene_id)) |>
  rename_with(.col = c("gene_id", "gene_name"), ~c("ensembl", "symbol"))

lapply(c("metadata", "count_mat", "gene_cnvan_key"), function(x) write_tsv(get(x), file.path(resDir, paste0(prefix, "_", x, ".tsv"))))